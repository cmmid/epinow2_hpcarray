
suppressPackageStartupMessages({
  require(data.table)
  require(qs)
})

.debug <- "ZAF"
.args <- if (interactive()) sprintf(c(
  "./output_data/%s/result.rds",
  "~/Dropbox/covidm_reports/hpc_inputs/%s/params_set.rds",
  "~/Dropbox/covidm_reports/hpc_inputs/%s/contact_matrices.rds",
  "~/Dropbox/covidm_reports/hpc_inputs/covidm_fit_yu.qs",
  "%s.rds"
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

cm_calc_R0_extended <- function(
  params,
  umultiplier = 1,
  cm_reds = c(0, 0, 0, 0),
  fIsred = rep(0, 16)
) {
  
  infected_states <- c("E","Ia","Ip","Is")
  infected_states_entry <- c(1,0,0,0)
  ages <- c(1:length(params$pop[[1]]$size))
  pops <- sapply(params$pop, "[[", "name")
  
  duration <- function(distributed_times, tstep=params$time_step){
    #calculates the mean of the distribution
    sum(distributed_times * seq(0, by=tstep, length.out = length(distributed_times)))
  }
  
  #reduced transmission matrix
  #transmission matrix T times inverse of auxiliary matrix E
  transmission_reduced <- matrix(
    0,
    sum(infected_states_entry)*length(ages)*length(pops),
    length(infected_states)*length(ages)*length(pops)
  )
  
  #reduced transition matrix
  #negative of inversed transition matrix Sigma times auxilliary matrix E
  transition_reduced <- matrix(
    0,
    length(infected_states)*length(ages)*length(pops),
    sum(infected_states_entry)*length(ages)*length(pops)
  )
  
  for(p1 in 1:length(params$pop)){
    cm = Reduce('+', mapply(
      function(c, m, reduction) c * m * (1-reduction),
      params$pop[[p1]]$contact,
      params$pop[[p1]]$matrices,
      cm_reds,
      SIMPLIFY = F)
    )
    for(a1 in 1:length(params$pop[[p1]]$size)){
      for(p2 in 1:length(params$pop)){
        for(a2 in 1:length(params$pop[[p2]]$size)){
          for(s in 1:length(infected_states)){
            sj <- ti <- (p1-1)*length(params$pop[[p1]]$size)+a1
            si <- tj <- (p2-1)*length(params$pop[[p2]]$size)*length(infected_states)+(a2-1)*length(infected_states)+s
            
            trates <- c(
              "E" = 0,
              "Ia" =  params$pop[[p1]]$u[a1] * umultiplier *
                ifelse(
                  params$pop[[p1]]$size[a1] != 0,
                  params$pop[[p2]]$size[a2]/params$pop[[p1]]$size[a1],
                  0
                ) *
                #adjust beta if population size is scaled down
                # as probability with which people are contacted will be scale down if there
                # are less people of that age
                ifelse(
                  !is.null(params$pop[[p2]]$size_original),
                  params$pop[[p2]]$size[a2]/params$pop[[p2]]$size_original[a2],
                  1
                ) *
                cm[a1,a2] * 
                params$pop[[p2]]$fIa[a2] *
                params$travel[p1,p2] *
                ifelse(p1==p2,1,params$pop[[p2]]$tau[a2]),
              "Ip" = params$pop[[p1]]$u[a1] * umultiplier *
                ifelse(
                  params$pop[[p1]]$size[a1] != 0,
                  params$pop[[p2]]$size[a2]/params$pop[[p1]]$size[a1],
                  0
                ) *
                #adjust beta if population size is scaled down
                ifelse(
                  !is.null(params$pop[[p2]]$size_original),
                  params$pop[[p2]]$size[a2]/params$pop[[p2]]$size_original[a2],
                  1
                ) *
                cm[a1,a2] * 
                params$pop[[p2]]$fIp[a2] *
                params$travel[p1,p2] *
                ifelse(p1==p2,1,params$pop[[p2]]$tau[a2]),
              "Is" = params$pop[[p1]]$u[a1] * umultiplier *
                ifelse(
                  params$pop[[p1]]$size[a1] != 0,
                  params$pop[[p2]]$size[a2]/params$pop[[p1]]$size[a1],
                  0
                ) *
                #adjust beta if population size is scaled down
                ifelse(
                  !is.null(params$pop[[p2]]$size_original),
                  params$pop[[p2]]$size[a2]/params$pop[[p2]]$size_original[a2],
                  1
                ) *
                cm[a1,a2] * 
                params$pop[[p2]]$fIs[a2] * (1-fIsred)[a2] *
                params$travel[p1,p2] *
                ifelse(p1==p2,1,params$pop[[p2]]$tau[a2])
            )
            
            #only applicable within one age
            if(p1==p2 & a1==a2){
              srates <- c(
                "E" = duration(params$pop[[p1]]$dE),
                "Ia" = (1-params$pop[[p1]]$y[a1])*duration(params$pop[[p1]]$dIa),
                "Ip" = params$pop[[p1]]$y[a1]*duration(params$pop[[p1]]$dIp),
                "Is" = params$pop[[p1]]$y[a1]*duration(params$pop[[p1]]$dIs)
              )  
            } else ( srates <- rep(0, 4) )
            
            transmission_reduced[ti,tj] <- trates[s]
            transition_reduced[si,sj] <- srates[s]
          }
        }
      }
    }
  }
  k <- transmission_reduced %*% transition_reduced
  #if(is.complex(eigen(k)$values)){
  #  warning("Eigenvalue is complex")
  #}
  return(max(Re(eigen(k)$values)))
}

scens <- data.table(
  expand.grid(
    home = c("none", "some"),
    work = c("small", "large"),
    school = c("small", "large"),
    other = c("small", "large"),
    symptrans = c("none","some"),
    # qtile = 1:length(preR0),
    stringsAsFactors = FALSE
  )
)[!(work == "small" & school == "small" & other == "small")]
#' none = 0 reduction
#' small < large

R0_recalc <- function(ps, umul, poppars, ls_cats) {
  names(ps)[1] <- "large"
  arglist <- modifyList(
    list(smallfact = 0, homefact = 0, symptfact = 0),
    as.list(ps)
  )
  cmr <- rep(0, 4)
  cmr[1] <- arglist$large*arglist$home
  cmr[1+which(ls_cats == "large")] <- arglist$large
  cmr[1+which(ls_cats == "small")] <- arglist$large*arglist$smallfact
  cm_calc_R0_extended(
    poppars, umul, cmr,
    fIsred = rep(arglist$large*arglist$symptfact, 16)
  )
}

cost_fun <- function(ps, ls_cats) {
  names(ps)[1] <- "large"
  arglist <- modifyList(
    list(smallfact = 0, homefact = 0, symptfact = 0),
    as.list(ps)
  )
  numlarge <- sum(ls_cats == "large")
  return(
    10-(numlarge*log(1-arglist$large) + (3-numlarge)*log(1-arglist$smallfact) + log(1-arglist$symptfact) + log(1-arglist$homefact))
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
  
  ured <- preR / cm_calc_R0_extended(poppars)
  
  scenarios[,{
    optpars <- c(large = refbaseline)
    if ("small" %in% c(work, school, other)) {
      optpars <- c(optpars, smallfact = factbaseline)
      if (home == "none") {
          if (symptrans != "none") optpars <- c(optpars, symptfact = factbaseline)
      } else {
        optpars <- c(optpars, homefact = factbaseline)
        if (symptrans != "none") optpars <- c(optpars, symptfact = factbaseline)
      }
    } else {
      if (home == "none") {
        if (symptrans != "none") optpars <- c(optpars, symptfact = factbaseline)
      } else {
        optpars <- c(optpars, homefact = factbaseline)
        if (symptrans != "none") optpars <- c(optpars, symptfact = factbaseline)
      }
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
      names(opars)[1] <- "large"
    # }
    
    res <- within(modifyList(
      list(smallfact = 0, homefact = 0, symptfact = 0),
      as.list(opars)
    ), {
      smallred <- large*smallfact
      homered <- large*homefact
      largered <- large
      sympred <- large*symptfact
    })
    
    rcalc <- sapply(ured, function(u)
      R0_recalc(opars, u, poppars = poppars, ls_cats = force(c(work, school, other)))
    )
    names(rcalc) <- sprintf("R.%s", c("ll","lo","md","hi","hh"))
    setTxtProgressBar(prg, .GRP)
    c(res[c("largered","smallred","homered","sympred")], as.list(rcalc))
  }, by=.(home, work, school, other, symptrans)]

}

saveRDS(fitfun(pars, preR0, postR0, scens), tail(.args, 1))
