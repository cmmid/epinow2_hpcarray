
suppressPackageStartupMessages({
  require(data.table)
  require(qs)
})

.debug <- "ZAF"
.args <- if (interactive()) sprintf(c(
  "~/Dropbox/Covid_LMIC/All_Africa_paper/r0_fitting/%s/result.rds",
  "~/Dropbox/covidm_reports/hpc_inputs/%s/params_set.rds",
  "~/Dropbox/covidm_reports/hpc_inputs/%s/contact_matrices.rds",
  "~/Dropbox/covidm_reports/hpc_inputs/covidm_fit_yu.qs",
  "~/Dropbox/Covid_LMIC/All_Africa_paper/r0_fitting/fits/alt_%s.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

R0ref <- readRDS(.args[1])
pars <- readRDS(.args[2])[[1]]
cmatrix <- readRDS(.args[3])
if (length(cmatrix) > 4) {
  cmatrix <- lapply(Reduce(function(l, r) {
    mapply("+", l, r, SIMPLIFY = FALSE)
  }, cmatrix), function(cm) cm/length(cmatrix))
}
yu_fits <- qread(.args[4])[order(ll)]
yu_fits[, eqs := (1:.N)/.N ]
#' using the median yu fits
medyu <- yu_fits[which.max(eqs > 0.5)]
yref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("y_",colnames(medyu))]))
uref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("u_",colnames(medyu))]))
ys <- rep(yref[1, ], each = 2)
us <- rep(uref[1, ], each = 2)

pars$pop <- lapply(
  pars$pop,
  function(x){
    x$matrices <- cmatrix
    x$y <- ys
    x$u <- us
    return(x)
  }
)

#' TODO handle exceptions
preR0 <- R0ref[era == "pre", c(lo.lo, lo, med, hi, hi.hi)]
postR0 <- R0ref[era == "post", c(lo.lo, lo, med, hi, hi.hi)]

#' TODO for each scenario
#' get the single parameter combination that minimizes
#' the difference between resulting shifts from preR0 => postR0

load("NGM.rda")

scens <- data.table(
  expand.grid(
#    home = c("none", "some"),
    work = c("small", "large"),
    school = c("small", "large"),
    other = c("small", "large"),
#    symptrans = c("none","some"),
    # qtile = 1:length(preR0),
    stringsAsFactors = FALSE
  )
)[!(work == "small" & school == "small" & other == "small")]
#' none = 0 reduction
#' small < large

fillargs <- function(ps) {
  names(ps)[1:2] <- c("large", "symptfact")
  return(modifyList(
    list(smallfact = 0), as.list(ps)
  ))
}

R0_recalc <- function(ps, umul, poppars, ls_cats) {
  arglist <- fillargs(ps)
  cmr <- rep(0, 4)
#  cmr[1] <- arglist$large*arglist$home
  cmr[1+which(ls_cats == "large")] <- arglist$large
  cmr[1+which(ls_cats == "small")] <- arglist$large*arglist$smallfact
  cm_ngm(
    poppars, umul, cmr, arglist$large*arglist$symptfact
  )$R0
}

cost_fun <- function(ps, ls_cats) {
  arglist <- fillargs(ps)
  numlarge <- sum(ls_cats == "large")
  return(
    10-(
      numlarge*log(1-arglist$large) +
      (3-numlarge)*log(1-arglist$large*arglist$smallfact) +
      log(1-arglist$large*arglist$symptfact)
    )
  )
}

optimTar <- function(ps, umul, poppars, ls_cats, tarR) {
  sse <- 0
  for (i in 1:length(umul)) {
    sse <- sse + (R0_recalc(ps, umul[i], poppars, ls_cats) - tarR[i])^2
  }
  # cost = closer to 100% reduction is increasingly expensive
  sse*cost_fun(ps, ls_cats)
}

factbaseline <- 0.5
refbaseline <- 0.5

prg <- txtProgressBar(max = scens[,.N], style = 3)

fitfun <- function(poppars, preR, postR, scenarios) {
  
  ured <- preR / cm_ngm(poppars)$R0
  
  scenarios[,{
    optpars <- c(large = refbaseline, symptfact = factbaseline)
    if ("small" %in% c(work, school, other)) {
      optpars <- c(optpars, smallfact = factbaseline)
    }
    optup <- rep(.99, length(optpars))
    optlow <- rep(0, length(optpars))
    # if (length(optpars) == 1) {
    #   opars <- optimize(f = optimizeTar, interval = c(0, 1),
    #     umul = ured,#[qtile],
    #     poppars = poppars, ls_cats = force(c(work, school, other)),
    #     tarR = postR
    #   )$minimum
    # } else {
      opars <- optim(
        optpars, optimTar, umul = ured,
        poppars = poppars, ls_cats = force(c(work, school, other)),
        tarR = postR, lower = optlow, upper = optup,
        method = if (length(optpars) == 1) "Brent" else "L-BFGS-B"
      )$par
      res <- within(fillargs(opars),{
        smallred <- large*smallfact
        largered <- large
        sympred <- large*symptfact
      })
      #names(opars)[1] <- "large"
    # }
    
    rcalc <- sapply(ured, function(u)
      R0_recalc(opars, u, poppars = poppars, ls_cats = force(c(work, school, other)))
    )
    names(rcalc) <- sprintf("R.%s", c("ll","lo","md","hi","hh"))
    setTxtProgressBar(prg, .GRP)
    c(res[c("largered","smallred","sympred")], as.list(rcalc))
  }, by=.(work, school, other)]

}

ret <- fitfun(pars, preR0, postR0, scens)

saveRDS(ret, tail(.args, 1))
