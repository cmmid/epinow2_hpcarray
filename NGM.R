
.args <- if (interactive()) c(
  "NGM.rds"
) else commandArgs(trailingOnly = TRUE)

ngm <- function(
  pops,
  
) {
  
}

saveRDS(ngm, tail(.args, 1))