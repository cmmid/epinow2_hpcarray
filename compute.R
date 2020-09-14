require(EpiNow2)
require(data.table)

.debug <- "BDI"
.args <- if (interactive()) sprintf(c(
  "cases.rds", "rt_bounds.rds", "%s", "res/%s/result.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

iso <- .args[3]

bounds <- as.list(readRDS(.args[2])[
  iso3 == iso
][, .SD ])

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
crs <- 4
smps <- 1e4

if (length(bounds$del) & (bounds$del > 0)) {
  
  reported_cases <- readRDS(.args[1])[
    (iso3 == iso) & (date <= bounds$interventionend + est.window)
  ][, .(date, confirm = cases )]
  
  # Add breakpoints
  reported_cases[,
    breakpoint := between(
     date,
     bounds$transitionstart, bounds$transitionend
    ) | (date >= bounds$interventionend)
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
    verbose = TRUE, return_fit = TRUE
  )$samples[variable == "R", .(value), by=.(sample, date)]
  
  results <- fbkp[
    (abs(date - (bounds$transitionstart-1)) < 0.5) |
    (abs(date - (bounds$transitionend+1)) < 0.5)
  ][,{
    qs <- quantile(value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
    as.list(qs)
  }, keyby = .(date)]
  
  results[order(date), era := c(ifelse(bounds$del > 7,"pre","weak") ,"post")]
} else {
  reported_cases <- readRDS(.args[1])[
    (iso3 == iso)
  ][1:est.window][, .(date, confirm = cases, breakpoint = FALSE )]
  
  fbkp <- estimate_infections(
    reported_cases, family = "negbin",
    generation_time = generation_time,
    delays = list(incubation_period, reporting_delay),
    samples = smps, warmup = smps*0.1, cores = crs,
    chains = crs, fixed = TRUE, horizon = 0,
    verbose = TRUE, return_fit = TRUE
  )$samples[variable == "R", .(value), by=.(sample, date)]

  results <- fbkp[
    date == min(date)
  ][,{
    qs <- quantile(value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
    as.list(qs)
  }, keyby = .(date)]
  
  results[order(date), era := "post" ]
  
}

saveRDS(results, tail(.args, 1))
