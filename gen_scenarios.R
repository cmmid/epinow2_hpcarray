suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("output_data", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/fits/%s.rds", "%s/%s/result.rds", "~/Dropbox/covidm_reports/hpc_inputs/ZAF/timing.rds", "%s/scens/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

#' parameter fitting results
fit.dt <- readRDS(.args[1])
#' intervention start date
ts.dt <- readRDS(.args[2])
start_date <- ts.dt[era == "window", round(median(date))]

day0 <- readRDS(.args[3])$day0date
#' target output
outfile <- tail(.args, 1)

translated.dt <- fit.dt[,.(
  self_iso = sympred,
  school = fifelse(school == "large", largered, smallred),
  home = 0,
  work = fifelse(work == "large", largered, smallred),
  other = fifelse(other == "large", largered, smallred),
  start_day = start_date - day0,
  end_day = Inf,
  travel = 0,
  population = -1, coverage = 1,
  scen_type = "mitigated"
)]

dg <- function(...) data.table(expand.grid(..., stringsAsFactors = FALSE))

ref <- function(
  age_split = NA,
  self_iso = 0,
  population = -1, coverage = 1, # coverage generally only applies to multi-pop models
  school = 0, home = 0, work = 0, other = 0, # recall == these are reductions
  travel = 0,
  start_day = NA_integer_, end_day = Inf,
  trigger_type = NA_character_,
  trigger_value = NA_real_,
  ...
) do.call(dg, c(as.list(environment()), list(...)))

scen_counter <- 1
tagscenario <- function(dt, sc) { 
  dt[, scen_id := (1:.N) + sc ]
  return(dt[.N, scen_id])
}

#' TODO: ignore? going to just iterate over rows where scen_id = 1, and if it's empty, nothing to do
unmitigated <- ref(scen_type = "unmitigated")[, scen_id := 1 ]

all.dt <- rbind(
  unmitigated, translated.dt, fill = TRUE
)[, scen_id := 1L:.N ]

saveRDS(all.dt, outfile)
