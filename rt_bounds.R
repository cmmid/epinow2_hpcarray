
suppressPackageStartupMessages({
  require(data.table)
})

.debug <- "."
.args <- if (interactive()) sprintf(c(
  "%s/cases.rds", "%s/interventions.rds", "%s/rt_bounds.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

case.dt <- readRDS(.args[1])
int.dt <- readRDS(.args[2])

jt <- case.dt[
  (cases + deaths) != 0
][, .(firstcases = date[1]), keyby=.(continent, iso3)][int.dt, on=.(iso3)][
  !is.na(continent),
  .(
    del = start - firstcases,
    from = firstcases,
    transitionstart = start,
    transitionend = end,
    interventionend = periodend
  ),
  keyby=.(continent, iso3)
]

saveRDS(jt, tail(.args, 1))

#' @examples 
#' require(ggplot2)
#' ggplot(jt) +
#'  aes(del, fill = ifelse(del <= 0, "none", ifelse(del > 7, ">week","<week"))) +
#'  facet_grid(. ~ continent) +
#'  geom_histogram(binwidth=1) +
#'  scale_fill_discrete("time between first cases\nandintervention") +
#'  theme_minimal()
#' 
#' jt[!is.na(continent), .(del=start - firstcases, continent)][,
#'   .N, keyby=.(continent, duration = ifelse(del <= 0, "none", ifelse(del > 7, ">week","<week")))
#' ]