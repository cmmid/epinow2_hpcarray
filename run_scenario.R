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
    return(x)
  }
)

params$pop <- lapply(
  params$pop,
  function(x){
    x$y <- ys
    x$u <- us
    return(x)
  }
)

params_back <- params

run_options <- melt(
  readRDS(.args[4])[era == "pre"],
  measure.vars = c("lo.lo","lo","med","hi","hi.hi"),
  value.name = "r0"
)[, model_seed := 1234L ]

cm_calc_R0_extended <- function(
  params
){
  
  infected_states <- c("E","Ia","Ip","Is")
  infected_states_entry <- c(1,0,0,0)
  ages <- c(1:length(params$pop[[1]]$size))
  pops <- sapply(params$pop, "[[", "name")
  
  duration <- function(distributed_times, tstep=params$time_step){
    #calculates the mean of the distribution
    sum(distributed_times * seq(0, by=tstep, length.out = length(distributed_times)))
  }
  
  #reduced transmission matrix
  #transmission matrix T times inverse of auxiliary matrix E
  transmission_reduced <- matrix(
    0,
    sum(infected_states_entry)*length(ages)*length(pops),
    length(infected_states)*length(ages)*length(pops)
  )
  
  #reduced transition matrix
  #negative of inversed transition matrix Sigma times auxilliary matrix E
  transition_reduced <- matrix(
    0,
    length(infected_states)*length(ages)*length(pops),
    sum(infected_states_entry)*length(ages)*length(pops)
  )
  
  for(p1 in 1:length(params$pop)){
    cm = Reduce('+', mapply(function(c, m) c * m, params$pop[[p1]]$contact, params$pop[[p1]]$matrices, SIMPLIFY = F))
    for(a1 in 1:length(params$pop[[p1]]$size)){
      for(p2 in 1:length(params$pop)){
        for(a2 in 1:length(params$pop[[p2]]$size)){
          for(s in 1:length(infected_states)){
            sj <- ti <- (p1-1)*length(params$pop[[p1]]$size)+a1
            si <- tj <- (p2-1)*length(params$pop[[p2]]$size)*length(infected_states)+(a2-1)*length(infected_states)+s
            
            trates <- c(
              "E" = 0,
              "Ia" =  params$pop[[p1]]$u[a1] *
                ifelse(
                  params$pop[[p1]]$size[a1] != 0,
                  params$pop[[p2]]$size[a2]/params$pop[[p1]]$size[a1],
                  0
                ) *
                #adjust beta if population size is scaled down
                # as probability with which people are contacted will be scale down if there
                # are less people of that age
                ifelse(
                  !is.null(params$pop[[p2]]$size_original),
                  params$pop[[p2]]$size[a2]/params$pop[[p2]]$size_original[a2],
                  1
                ) *
                cm[a1,a2] * 
                params$pop[[p2]]$fIa[a2] *
                params$travel[p1,p2] *
                ifelse(p1==p2,1,params$pop[[p2]]$tau[a2]),
              "Ip" = params$pop[[p1]]$u[a1] *
                ifelse(
                  params$pop[[p1]]$size[a1] != 0,
                  params$pop[[p2]]$size[a2]/params$pop[[p1]]$size[a1],
                  0
                ) *
                #adjust beta if population size is scaled down
                ifelse(
                  !is.null(params$pop[[p2]]$size_original),
                  params$pop[[p2]]$size[a2]/params$pop[[p2]]$size_original[a2],
                  1
                ) *
                cm[a1,a2] * 
                params$pop[[p2]]$fIp[a2] *
                params$travel[p1,p2] *
                ifelse(p1==p2,1,params$pop[[p2]]$tau[a2]),
              "Is" = params$pop[[p1]]$u[a1] *
                ifelse(
                  params$pop[[p1]]$size[a1] != 0,
                  params$pop[[p2]]$size[a2]/params$pop[[p1]]$size[a1],
                  0
                ) *
                #adjust beta if population size is scaled down
                ifelse(
                  !is.null(params$pop[[p2]]$size_original),
                  params$pop[[p2]]$size[a2]/params$pop[[p2]]$size_original[a2],
                  1
                ) *
                cm[a1,a2] * 
                params$pop[[p2]]$fIs[a2] *
                params$travel[p1,p2] *
                ifelse(p1==p2,1,params$pop[[p2]]$tau[a2])
            )
            
            #only applicable within one age
            if(p1==p2 & a1==a2){
              srates <- c(
                "E" = duration(params$pop[[p1]]$dE),
                "Ia" = (1-params$pop[[p1]]$y[a1])*duration(params$pop[[p1]]$dIa),
                "Ip" = params$pop[[p1]]$y[a1]*duration(params$pop[[p1]]$dIp),
                "Is" = params$pop[[p1]]$y[a1]*duration(params$pop[[p1]]$dIs)
              )  
            } else ( srates <- rep(0, 4) )
            
            transmission_reduced[ti,tj] <- trates[s]
            transition_reduced[si,sj] <- srates[s]
          }
        }
      }
    }
  }
  k <- transmission_reduced %*% transition_reduced
  #if(is.complex(eigen(k)$values)){
  #  warning("Eigenvalue is complex")
  #}
  return(max(Re(eigen(k)$values)))
}

allbind <- data.table()

for (scenario_index in 1:max(scenario$scen_id)) {
  #' sub
  iv_data <- scenario[scen_id == scenario_index][order(trigger_type)]
  for(i in 1:nrow(run_options)){
    
    params <- params_back
    
    #adjust r0 to that in current sample
    target_R0 <- run_options[i, r0]
    uf <- target_R0 / cm_calc_R0_extended(params)
    # params$pop <- lapply(
    #   params$pop,
    #   function(x){
    #     x$u <- x$u * uf
    #     return(x)
    #   }
    # )
    
    if (iv_data[,.N]) {
      iv = cm_iv_build(params)
      
      # generic interventions
      for (j in 1:nrow(iv_data[population == -1])) {
        #pars <- as.list(iv_data[population == -1][j])
        
        with(as.list(iv_data[population == -1][j]), {
          contact <- c(home, work, school, other)
          cm_iv_set(iv,
                    as.Date(params$date0) + start_day,
                    as.Date(params$date0) + ifelse(is.finite(end_day),end_day,1e3),
                    fIs = rep(1-self_iso, 16),
                    contact = 1-contact # TODO: manage splits
          )
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
  }
  
}

qsave(allbind, tarfile)
