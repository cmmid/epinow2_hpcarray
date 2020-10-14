
.args <- if (interactive()) c(
  "~/Dropbox/Covid_LMIC/All_Africa_paper/r0_fitting"
) else commandArgs(trailingOnly = TRUE)

fls <- grep("[[:upper:]]{3}/res",
  list.files(.args[1], "result\\.rds$", recursive = TRUE),
  value = TRUE
)

cups <- fls[sapply(fls, function(fl) readRDS(
  sprintf("%s/%s", .args[1], fl)
)[era == "window"][, !all(sign(diff(med))<1) ]
)]
