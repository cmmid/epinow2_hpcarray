require(data.table)

.args <- if (interactive()) c(
  "cases.rds", "iso3.csv"
) else commandArgs(trailingOnly = TRUE)

target <- tail(.args, 1)
res <- readRDS(.args[1])[, .(iso3=unique(iso3))]

record <- if (file.exists(target)) {
  ref <- fread(target)
  (res[,.N] != ref[,.N]) || !all(res == ref)
} else TRUE

if (record) {
  fwrite(res, tail(.args, 1), col.names = FALSE)
} else warning(sprintf("current %s is identical", target))
