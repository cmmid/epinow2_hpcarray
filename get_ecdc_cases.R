require(data.table)

.args <- if (interactive()) c(
  "cases.rds"
) else commandArgs(trailingOnly = TRUE)

target <- tail(.args, 1)
ecdcurl <- "https://opendata.ecdc.europa.eu/covid19/casedistribution/csv"

res <- fread(ecdcurl)
res[, date := as.Date(dateRep, format = "%d/%m/%Y") ]

final <- res[countryterritoryCode != "",
  .(cases, deaths),
  keyby=.(continent = continentExp, iso3=countryterritoryCode, date)
]

if (final[, any(cases < 0)]) {
  warning(sprintf(
    "the following countries have negative cases:\n%s",
    final[, .(neg=any(cases < 0)), by=iso3][neg==TRUE, paste(iso3, collapse = "\n")]
  ))
  excludeisos <- final[cases < 0, unique(iso3)]
  final <- final[!(iso3 %in% excludeisos)]
}

record <- if (file.exists(target)) {
  ref <- readRDS(target)
  (final[,.N] != ref[,.N]) || !all(ref == final)
} else TRUE

if (record) {
  saveRDS(final, target)
} else warning(sprintf("current %s is identical", target))