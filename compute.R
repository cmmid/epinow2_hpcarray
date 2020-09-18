require(EpiNow2)
require(data.table)

.debug <- "SWE"
.args <- if (interactive()) sprintf(c(
  "cases.rds", "rt_bounds.rds", "%s", "2", "res/%s/result.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

iso <- .args[3]

bounds <- as.list(readRDS(.args[2])[
  iso3 == iso
])

# Set up example generation time
generation_time <- as.list(EpiNow2::covid_generation_times[1,
  .(mean, mean_sd, sd, sd_sd, max=30)
])

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
crs <- as.integer(.args[4])
smps <- 1e4

if (length(bounds$del) && !is.na(bounds$del) && (bounds$del > 0)) {
  
  reported_cases <- readRDS(.args[1])[
    (iso3 == iso) & (date <= bounds$interventionend + est.window)
  ][, .(date, confirm = cases )]
  
  reported_cases[, era := fifelse(bounds$del > 7, "pre", "short")  ]
  reported_cases[date >= bounds$transitionstart, era := "window" ]
  reported_cases[date > bounds$transitionend, era := "post" ]
  reported_cases[date >= bounds$interventionend, era := "censor" ]
  
  
  # Add breakpoints
  reported_cases[,
    breakpoint := era %in% c("window", "censor")
  ]
  
  # Run model with breakpoints but otherwise static Rt
  # This formulation may increase the apparent effect of
  # the breakpoint but needs to be tested using
  # model fit criteria (i.e LFO).
  fbkp <- estimate_infections(
    reported_cases, family = "negbin",
    generation_time = generation_time,
    delays = list(incubation_period, reporting_delay),
    samples = smps, warmup = smps*0.1, cores = crs,
    chains = crs, estimate_breakpoints = TRUE, fixed = TRUE, horizon = 0,
    verbose = FALSE, return_fit = TRUE
  )$samples[variable == "R", .(value), by=.(sample, date)]
  
  results <- fbkp[
    between(date, bounds$transitionstart-1, bounds$transitionend+1)
  ][, {
    qs <- quantile(value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
    as.list(qs)
  }, keyby = .(date)][
    reported_cases[,
      .(date, ccases = cumsum(confirm), era)
    ][between(date, bounds$transitionstart-1, bounds$transitionend+1)],
    on=.(date)
  ]
  
} else {
  
  reported_cases <- readRDS(.args[1])[
    (iso3 == iso)
  ][1:(2*est.window)][, .(date, confirm = cases, breakpoint = FALSE )]
 
  reported_cases[, era := "censor" ]
  reported_cases[1:est.window, era := "post" ]

   
  fbkp <- estimate_infections(
    reported_cases, family = "negbin",
    generation_time = generation_time,
    delays = list(incubation_period, reporting_delay),
    samples = smps, warmup = smps*0.1, cores = crs,
    chains = crs, fixed = TRUE, horizon = 0,
    verbose = FALSE, return_fit = TRUE
  )$samples[variable == "R", .(value), by=.(sample, date)]

  ref <- reported_cases[,
    .(date, ccases = cumsum(confirm), era)
  ][est.window]
  
  results <- fbkp[
    date == ref$date
  ][,{
    qs <- quantile(value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
    as.list(qs)
  }, keyby = .(date)][
    ref,
    on=.(date)
  ]
  
}

saveRDS(results, tail(.args, 1))
