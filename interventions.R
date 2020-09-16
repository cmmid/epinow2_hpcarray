require(data.table)
require(lubridate)

.args <- if (interactive()) c(
  "rawinterventions.csv", "interventions.rds"
) else commandArgs(trailingOnly = TRUE)

# TODO: move over intervention start timing analysis
# dayeff: is the estimate of the day the intervention becomes effective

ref <- fread(.args[1])[
  admin_level == "national" &
  measure_stage == "new" &
  enforcement %in% c("Required", "Monitored", "Not applicable") &
  !(who_category %in% c(
    "International travel measures", "Other measures", "Environmental measures"
  ))
]

first.days <- ref[,.(
  date_start = mean(as.Date(date_start, "%m-%d-%y"))
),by=.(iso3=iso, who_category, who_subcategory, who_measure)]

res <- first.days[,
  .(dayeff=mean(date_start)),
  keyby=iso3
]

# presumed delay of intervention impact
offset <- 0
# days to censor around intervention for estimation
window <- 7
# days of intervention period to consider
est.window <- 30

res[, start := round_date(dayeff + offset, "day") ]
res[, end := start + window ] # ceiling date
res[, periodend := end + est.window ]

saveRDS(res, tail(.args, 1))