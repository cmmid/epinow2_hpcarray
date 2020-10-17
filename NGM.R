
.args <- if (interactive()) c(
  "NGM.rds"
) else commandArgs(trailingOnly = TRUE)

ngm <- function(
  pops,
  u_multiplier = 1,
  cm_reductions = c(0, 0, 0, 0),
  fIs_reductions = rep(0, pops[[1]]$size)
) {
  
  
  
   
}

saveRDS(ngm, tail(.args, 1))