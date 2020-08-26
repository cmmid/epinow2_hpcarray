require(data.table)

.args <- if (interactive()) c(
  "interventions.rds"
) else commandArgs(trailingOnly = TRUE)

# TODO: move over intervention start timing analysis
# dayeff: is the estimate of the day the intervention becomes effective

res <- data.table(
  iso3 = "KEN",
  dayeff = as.Date("2020-03-31")
)
# days to censor around intervention for estimation
window <- 14
# days of intervention period to consider
est.window <- 30

res[,
  c("start", "end", "periodend") := as.list(
    dayeff + c(-window/2, window/2, window/2+est.window)
  )
]

saveRDS(res, tail(.args, 1))