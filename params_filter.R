
suppressPackageStartupMessages({
  require(data.table)
})

.debug <- "output_data"
.args <- if (interactive()) sprintf(c(
  "%s", "~/Dropbox/covidm_reports/hpc_inputs", "%s/eligible.csv"
), .debug) else commandArgs(trailingOnly = TRUE)

tars <- grep(
  "^\\w{3}",
  list.files(.args[1], "result\\.rds", recursive = TRUE),
  value = T
)

inputisos <- list.dirs(.args[2], full.names = F)[-1]

#' which of tars have pre / post estimates?

subtars <- data.table(iso=intersect(inputisos, gsub("/result\\.rds", "", tars[sapply(tars, function(fn) {
  res <- readRDS(sprintf("%s/%s", .args[1], fn))
  res[, era[1] == "pre" & era[.N] == "post"] & res[, med[1] > med[.N]]
})])))

fwrite(subtars, tail(.args, 1), col.names = F)
