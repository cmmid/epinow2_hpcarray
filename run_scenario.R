suppressPackageStartupMessages({
  require(data.table)
  require(qs)
})

.debug <- "ZAF"
.args <- if (interactive()) sprintf(c(
  "~/Dropbox/Covid_LMIC/All_Africa_paper/r0_fitting/scens/alt_%s.rds",
  "~/Dropbox/covidm_reports/hpc_inputs/%s/params_set.rds",
  "~/Dropbox/covidm_reports/hpc_inputs/covidm_fit_yu.qs",
  "~/Dropbox/Covid_LMIC/All_Africa_paper/r0_fitting/%s/alt_result.rds",
  "~/Dropbox/covidm_reports/hpc_inputs/%s/contact_matrices.rds",
  "intros/intros.rds",
  "intros/urban.rds",
  "../covidm",
  "~/Dropbox/Covid_LMIC/All_Africa_paper/r0_fitting/%s/alt_projection.qs"
), .debug) else commandArgs(trailingOnly = TRUE)

load("NGM.rda")

# load covidm
cm_path = tail(.args, 2)[1]
cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T
cm_version = 1

suppressPackageStartupMessages({
  source(paste0(cm_path, "/R/covidm.R"))
})

# identify country / scenario
scenario <- readRDS(.args[1])
tarfile <- tail(.args, 1)

yu_fits <- qread(.args[3])[order(ll)]
yu_fits[, eqs := (1:.N)/.N ]
#' using the median yu fits
medyu <- yu_fits[which.max(eqs > 0.5)]
yref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("y_",colnames(medyu))]))
uref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("u_",colnames(medyu))]))
ys <- rep(yref[1, ], each = 2)
us <- rep(uref[1, ], each = 2)
#if (.debug == "NGA") {
#  ushft <- max(us)/us
#  us <- 1/(1:length(us)+1)
#  ys <- ys / ushft
#}

params <- readRDS(.args[2])[[1]]

tariso <- gsub(".+([[:upper:]]{3})\\.rds$","\\1", .args[1])

intros <- readRDS(.args[6])[iso3 == tariso, Reduce(c, mapply(rep, intro.day, intros, SIMPLIFY = FALSE))]
urbfrac <- readRDS(.args[7])[iso3 == tariso, value / 100]
#urbfrac <- 1

#urbfrac <- if (.debug == "NGA") .51 else 1

params$pop[[1]]$seed_times <- intros
params$pop[[1]]$size <- round(params$pop[[1]]$size*urbfrac)
# params$pop[[1]]$dist_seed_ages <- c(rep(0,4), rep(1, 6), rep(0, 6))

contact_matrices <- readRDS(.args[5])

refcm <- if (is.null(names(contact_matrices))) { 
  contact_matrices <- lapply(Reduce(function(l, r) {
    mapply("+", l, r, SIMPLIFY = FALSE)
  }, contact_matrices), function(cm) cm/length(contact_matrices))
} else { contact_matrices }
names(refcm) <- gsub("cm_","",names(refcm))

params$pop <- lapply(
  params$pop,
  function(x){
    x$matrices <- refcm
    x$y <- ys
    x$u <- us
    return(x)
  }
)

run_options <- melt(
  readRDS(.args[4])[era == "pre"],
  measure.vars = c("lo.lo","lo","med","hi","hi.hi"),
  value.name = "r0"
)[, model_seed := 1234L ]

params_back <- params

allbind <- data.table()

scenario[, waning_dur := NA_integer_ ]
#scenario[!(scen_type == "unmitigated"), waning_dur := 90 ]

prg <- txtProgressBar(max = scenario[,.N]*run_options[,.N], style = 3)
prgind <- 0

for (scenario_index in 1:max(scenario$scen_id)) {
  #' sub
  iv_data <- scenario[scen_id == scenario_index][order(trigger_type)]
  for(i in 1:nrow(run_options)){
    
    params <- params_back
    
    #adjust r0 to that in current sample
    target_R0 <- run_options[i, r0]
    uf <- target_R0 / cm_ngm(params)$R0
    params$pop <- lapply(
      params$pop,
      function(x){
        x$u <- x$u * uf
        return(x)
      }
    )
    
    if (iv_data[!(scen_type == "unmitigated"), .N]) {
      iv = cm_iv_build(params)
      
      # generic interventions
      for (j in 1:nrow(iv_data[population == -1])) {
        #pars <- as.list(iv_data[population == -1][j])
        
        with(as.list(iv_data[population == -1][j]), {
          contact <- c(home, work, school, other)
          if (is.na(waning_dur)) {
            cm_iv_set(iv,
                      as.Date(params$date0) + start_day,
                      as.Date(params$date0) + ifelse(is.finite(end_day),end_day,params_back$time1),
                      fIs = rep(1-self_iso, 16),
                      contact = 1-contact # TODO: manage splits
            )
          } else {
            for (inc in (start_day:ifelse(is.finite(end_day),end_day-1,params_back$time1-1))-start_day) {
              waning_mult <- exp(-(1/waning_dur)*inc)
              cm_iv_set(iv,
                        as.Date(params$date0) + start_day + inc,
                        as.Date(params$date0) + start_day + inc + 1,
                        fIs = rep(1-self_iso*waning_mult, 16),
                        contact = 1-contact*waning_mult # TODO: manage splits
              )
            }
          }
        })
      }
      
      params = cm_iv_apply(params, iv)
      rm(iv)
      
    }
    #run the model
    sim <- cm_simulate(
      params, 1,
      model_seed = run_options[i, model_seed]
    )$dynamics
    
    result <- sim[,
      .(value = sum(value)),
      keyby = .(run, t, group, compartment)
    ][, run := i ][, scen_id := scenario_index ]
    
    rm(params)
    rm(sim)
    
    allbind <- rbind(allbind, result)
    prgind <- prgind + 1
    setTxtProgressBar(prg, prgind)
  }
  
}

#' @examples 
#' require(ggplot2)
#' ggplot(allbind[compartment == "cases"]) + aes(t, value) + facet_grid(group ~ ., scales = "free_y") + geom_line()

qsave(allbind, tarfile)
