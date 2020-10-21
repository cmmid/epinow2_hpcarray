suppressPackageStartupMessages({
  require(EpiNow2)
  require(data.table)
  require(ggplot2)
  require(patchwork)
})

.debug <- "ZAF"
.args <- if (interactive()) sprintf(c(
  "cases.rds",
  "rt_bounds.rds",
  "%s",
  "~/Dropbox/Covid_LMIC/All_Africa_paper/r0_fitting/%s/alt_result.rds"
),.debug) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

case.dt <- readRDS(.args[1])[iso3 == tariso][, .(date, confirm = cases )]
lims.dt <- readRDS(.args[2])[iso3 == tariso]
lims.dt[, transitionstart := transitionstart - 3 ]
#lims.dt[, transitionstart := as.Date("2020-03-23") ]
lims.dt[, transitionend := transitionstart + 14 ]
lims.dt[, interventionend := transitionend + 30 ]

#' @examples
#' ggplot(case.dt) + aes(date, confirm) + geom_line() +
#' geom_col(aes(x=transitionstart-shft, y=1000), data = lims.dt, alpha = 0.2) +
#' geom_col(aes(x=transitionend-shft, y=1000), data = lims.dt, alpha = 0.2) +
#' theme_minimal() + scale_y_log10() 

fill.case <- case.dt[
  case.dt[, .(date = seq(min(date),max(date),by="day"))],
  on=.(date),
  .(date, confirm = fifelse(is.na(confirm), 0, confirm))
]

generation_time <- as.list(
  EpiNow2::covid_generation_times[
    1, .(mean, mean_sd, sd, sd_sd, max=30)
])

generation_time$mean <- 6.790176

# Set delays between infection and case report
# (any number of delays can be specifed here)
incubation_period <- as.list(
  EpiNow2::covid_incubation_period[
    1, .(mean, mean_sd, sd, sd_sd, max=30)
])

reporting_delay <- list(
  mean = log(6), mean_sd = log(2),
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
  
  # Run model with breakpoints but otherwise static Rt
  # This formulation may increase the apparent effect of
  # the breakpoint but needs to be tested using
  # model fit criteria (i.e LFO).
# fbkp <- estimate_infections(
#   early_reported_cases, family = "negbin",
#   generation_time = generation_time,
#   delays = list(incubation_period, reporting_delay),
#   samples = smps, cores = crs,
#   chains = crs*2,
#   estimate_breakpoints = TRUE,
#   fixed = TRUE, horizon = 0,
#   verbose = TRUE, return_fit = FALSE
# )$samples[variable == "R", .(value), by=.(sample, date)]
# 
# nobrk <- estimate_infections(
#   early_reported_cases[,.(date, confirm)], family = "negbin",
#   generation_time = generation_time,
#   delays = list(incubation_period, reporting_delay),
#   samples = smps, cores = crs,
#   chains = crs*2,
#   horizon = 0, verbose = TRUE, return_fit = FALSE
# )$samples[variable == "R", .(value), by=.(sample, date)]
# 
# 
# results <- cmbn(fbkp, nobrk)

#' @examples 
#' xscale <- scale_x_date(
#'   date_breaks="week", date_labels = c("%b %d"),
#'   minor_breaks = NULL
#' )
#' 
#' p1 <- function(dt) ggplot(dt) + aes(date) +
#'   geom_ribbon(aes(ymin=lo.lo, ymax=hi.hi, fill=analysis), alpha = 0.2) +
#'   geom_ribbon(aes(ymin=lo, ymax=hi, fill=analysis), alpha = 0.2) +
#'   geom_line(aes(y=med, color=analysis)) +
#'   theme_minimal() +
#'   theme(
#' #    axis.text.x = element_blank(),
#' #    axis.title.x = element_blank()
#'   ) +
#'   xscale +
#'   scale_y_continuous("Rt")
#' 
#' p2 <- function(dt) ggplot(dt[date <= (lims.dt$interventionend + est.window)]) + aes(date, confirm) +
#'   geom_line() +
#'   theme_minimal() +
#'   theme(
#'     axis.text.x = element_text(angle = 90, vjust = .5)
#'   ) +
#'   xscale +
#'   scale_y_log10("Confirmed Cases")
#'
#' cmbn <- function(brks, non) {
#'   res <- rbind(
#'     brks[, analysis := "breakpoint"],
#'     non[, analysis := "continuous"]
#'   )
#'   res[, {
#'     qs <- quantile(value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
#'     names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
#'     as.list(qs)
#'   }, keyby = .(analysis, date)]
#' }
#'
#' (p1(results) / p2(fill.case)) & coord_cartesian(
#'   xlim = range(re.est$date)
#' )

rebreak <- copy(early_reported_cases)

# rebreak[era %in% c("window","introduction"), breakpoint := FALSE]
# rebreak[era %in% c("window", "pre"), era := "post"]
# rebreak[era == "introduction", era := "pre"]
# 
# cntr <- as.Date(lims.dt$transitionstart) - shft + 4
# 
# rebreak[
#   between(date, cntr-4, cntr+3),
#   era := "window"
# ]
# rebreak[date < cntr - 4, era := "pre" ]
# rebreak[era == "window", breakpoint := TRUE]
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

# 7-centre'd moving average estimate_week_eff
#' @examples 
#' re.cont <- estimate_infections(
#'   rebreak[, .(date, confirm = round(frollmean(confirm, 7, align = "center")))][!is.na(confirm)],
#'   family = "negbin",
#'   generation_time = generation_time,
#'   delays = list(incubation_period, reporting_delay),
#'   samples = smps, cores = crs,
#'   chains = crs*2,
#'   fixed = FALSE, horizon = 0,
#'   verbose = TRUE, return_fit = FALSE,
#'   estimate_week_eff = FALSE
#' )$samples[variable == "R", .(value), by=.(sample, date)]
#' 
#' moreres <- cmbn(re.est, re.cont)
#' smooth.cases <- case.dt[, .(date, confirm = round(frollmean(confirm, 7, align = "center")))][!is.na(confirm)]
#' 
#' (p1(moreres) / p2(smooth.cases)) & coord_cartesian(
#'   xlim = range(re.est$date)
#' )
# 
# (p1(results) + p1(moreres) + plot_layout(guides = "collect")) &
#   coord_cartesian(ylim = range(
#     rbind(results, moreres)[, c(lo.lo, hi.hi)]
# ))
