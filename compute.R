require(EpiNow2)
require(data.table)

.args <- if (interactive()) c(
  "interventions.rds", "cases.rds", "KEN", "./KEN/result.rds"
) else commandArgs(trailingOnly = TRUE)

iso <- .args[3]
intervention.window <- as.list(readRDS(.args[1])[
  iso3 == iso
][, .(start, end, periodend) ])
# additional time to include for algorithm
est.window <- 30
reported_cases <- readRDS(.args[2])[
  (iso3 == iso) & (date <= intervention.window$periodend + est.window)
][, .(date, confirm = cases )]


# Add breakpoints
reported_cases[,
  breakpoint := between(
    date,
    intervention.window$start, intervention.window$end
  ) | (date >= intervention.window$periodend)
]
# Set up example generation time
generation_time <- as.list(
  EpiNow2::covid_generation_times[
    1, .(mean, mean_sd, sd, sd_sd, max=30)
  ]
)

# Set delays between infection and case report
# (any number of delays can be specifed here)
incubation_period <- as.list(
  EpiNow2::covid_incubation_period[
    1, .(mean, mean_sd, sd, sd_sd, max=30)
  ]
)

reporting_delay <- list(
  mean = log(5), mean_sd = log(2),
  sd = log(2), sd_sd = log(1.5),
  max = 30
)

# Run model with breakpoints but otherwise static Rt
# This formulation may increase the apparent effect of
# the breakpoint but needs to be tested using
# model fit criteria (i.e LFO).
fbkp <- estimate_infections(
  reported_cases, family = "negbin",
  generation_time = generation_time,
  delays = list(incubation_period, reporting_delay),
  samples = 1000, warmup = 200, cores = 4,
  chains = 4, estimate_breakpoints = TRUE, fixed = TRUE, horizon = 0,
  verbose = TRUE, return_fit = TRUE
)$samples[variable == "R", .(value), by=.(sample, date)]

# pickout pre-intervention + post-intervention dates
results <- fbkp[
  date %in% c(intervention.window$start-1, intervention.window$end+1), {
    qs <- quantile(value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
    as.list(qs)
  },
  keyby = .(date)
]

results[order(date), era := c("pre","post")]

saveRDS(results, tail(.args, 1))
