suppressPackageStartupMessages({
  require(EpiNow2)
  require(data.table)
  require(qs)
})

.debug <- "ZAF"
.args <- if (interactive()) sprintf(c(
  "cases.rds", "rt_bounds.rds",
  "~/Dropbox/covidm_reports/hpc_inputs/%s/contact_matrices.rds",
  "~/Dropbox/covidm_reports/hpc_inputs/%s/params_set.rds",
  "~/Dropbox/covidm_reports/hpc_inputs/covidm_fit_yu.qs",
  "%s",
  "~/Dropbox/Covid_LMIC/All_Africa_paper/r0_fitting/%s/result.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

tariso <- .args[6]

case.dt <- readRDS(.args[1])[iso3 == tariso][, .(date, confirm = cases )]
fill.case <- case.dt[
  case.dt[, .(date = seq(min(date),max(date),by="day"))],
  on=.(date),
  .(date, confirm = fifelse(is.na(confirm), 0, confirm))
  ]

lims.dt <- readRDS(.args[2])[iso3 == tariso]
lims.dt[, transitionstart := transitionstart - 5 ]
#lims.dt[, transitionstart := as.Date("2020-03-23") ]
lims.dt[, transitionend := transitionstart + 14 ]
lims.dt[, interventionend := transitionend + 30 ]

contact_matrices <- readRDS(.args[3])

refcm <- if (is.null(names(contact_matrices))) { 
  contact_matrices <- lapply(Reduce(function(l, r) {
    mapply("+", l, r, SIMPLIFY = FALSE)
  }, contact_matrices), function(cm) cm/length(contact_matrices))
} else { contact_matrices }
names(refcm) <- gsub("cm_","",names(refcm))

params <- readRDS(.args[4])[[1]]

yu_fits <- qread(.args[5])[order(ll)]
yu_fits[, eqs := (1:.N)/.N ]
#' using the median yu fits
medyu <- yu_fits[which.max(eqs > 0.5)]
yref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("y_",colnames(medyu))]))
uref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("u_",colnames(medyu))]))
ys <- rep(yref[1, ], each = 2)
us <- rep(uref[1, ], each = 2)

params$pop <- lapply(
  params$pop,
  function(x){
    x$matrices <- refcm
    x$y <- ys
    x$u <- us
    return(x)
  }
)

load("NGM.rda")

# Set up example generation time
generation_time <- as.list(EpiNow2::covid_generation_times[1,
  .(mean, mean_sd, sd, sd_sd, max=30)
])

generation_time$mean <- unname(cm_generation_time(params))

# Set delays between infection and case report
# (any number of delays can be specifed here)
incubation_period <- as.list(EpiNow2::covid_incubation_period[1,
  .(mean, mean_sd, sd, sd_sd, max=30)
])

reporting_delay <- list(
  mean = log(5), mean_sd = log(2),
  sd = log(2), sd_sd = log(1.5),
  max = 30
)

# additional time to include for algorithm
est.window <- 30
smps <- 2e3
crs <- 4

early_reported_cases <- case.dt[date <= (lims.dt$interventionend + est.window)]
with(lims.dt,{
  early_reported_cases[, era := fifelse(del > 7, "pre", "short")  ]
  early_reported_cases[date >= transitionstart, era := "window" ]
  early_reported_cases[date > transitionend, era := "post" ]
  early_reported_cases[date >= interventionend, era := "censor" ]
})

# Add breakpoints
early_reported_cases[, breakpoint := era %in% c("window", "censor") ]
rebreak <- copy(early_reported_cases)
rebreak[cumsum(confirm) <= 10, era := "introduction"]
rebreak[era == "introduction", breakpoint := TRUE ]

re.est <- estimate_infections(
  rebreak, family = "negbin",
  generation_time = generation_time,
  delays = list(incubation_period, reporting_delay),
  samples = smps, cores = crs,
  chains = crs,
  estimate_breakpoints = TRUE,
  fixed = TRUE, horizon = 0,
  verbose = TRUE, return_fit = TRUE
)$samples[variable == "R", .(value), by=.(sample, date)]

results <- re.est[
  between(date, lims.dt$transitionstart-1, lims.dt$transitionend+1)
  ][, {
    qs <- quantile(value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
    as.list(qs)
  }, keyby = .(date)][
    rebreak[,
            .(date, ccases = cumsum(confirm), era)
            ][
              between(date, lims.dt$transitionstart-1, lims.dt$transitionend+1)
              ],
    on=.(date)
    ]

saveRDS(results, tail(.args, 1))
