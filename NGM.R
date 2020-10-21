
require(Matrix)
require(RSpectra)

.args <- if (interactive()) c(
  "NGM.rds"
) else commandArgs(trailingOnly = TRUE)

dur <- function(dX, ts=0.25) sum(dX*seq(0,by=ts,length.out=length(dX)))

#' assumes a single population; doesn't work for multi-pop model w/
#' travel matrix
ngm <- function(
  pop,
  u_multiplier = 1,
  cm_reductions = c(0, 0, 0, 0),
  fIs_reductions = 0
) {
  
  stay_vec <- function(dX) {
    moveX <- dX[which(dX != 0)]
    #' at index i, how much probability do we have left?
    ref <- c(1,(1-cumsum(moveX))[-length(moveX)])
    #' p(move at i) = dX[i]/remaining probability
    pmax(1 - moveX/ref, 0)
  }
  
  lblr <- function(fmt, vec) sprintf(fmt, (1:length(vec))-1)
  
  redfIs <- pop$fIs*(1 - fIs_reductions)
  redCon <- pop$contact*(1 - cm_reductions)
  
  umod <- u_multiplier*pop$u
  
  #' make the transition matrices
  #' E -> Ia or Ip
  #' Ip -> Is
  #' dX[i] is the probability of having the duration of X being i-1
  #' so if someone is in X[i], dX[i] is the probability that they *exit* X
  #' (except for dX[i==max(i)])
  #' otherwise they go to i+1 with probability dX[i]
  E_stay <- stay_vec(pop$dE)
  Ia_stay <- stay_vec(pop$dIa)
  Ip_stay <- stay_vec(pop$dIp)
  Is_stay <- stay_vec(pop$dIs)
  
  len_allstates <- sum(
    length(E_stay), length(Ia_stay), length(Ip_stay), length(Is_stay)
  )
  
  Enms <- lblr("E%03i", E_stay)
  Ianms <- lblr("Ia%03i", Ia_stay)
  Ipnms <- lblr("Ip%03i", Ip_stay)
  Isnms <- lblr("Is%03i", Is_stay)
  
  nms <- c(
    Enms, Ianms, Ipnms, Isnms
  )
  
  Tij <- Matrix(
    0,
    nrow = len_allstates, ncol = len_allstates,
    dimnames = list(to=nms, from=nms),
    sparse = TRUE
  )
  
  #' all the states move to their next time step with the stay probability
  diag(
    Tij[Enms[-1], Enms[-length(E_stay)]]
  ) <- E_stay[-length(E_stay)]
  diag(
    Tij[Ianms[-1], Ianms[-length(Ia_stay)]]
  ) <- Ia_stay[-length(Ia_stay)]
  diag(
    Tij[Ipnms[-1], Ipnms[-length(Ip_stay)]]
  ) <- Ip_stay[-length(Ip_stay)]
  diag(
    Tij[Isnms[-1], Isnms[-length(Is_stay)]]
  ) <- Is_stay[-length(Is_stay)]
  
  #' presymptomatic moves to step 0 of symptomatic with it's move prob
  Tij[Isnms[1], Ipnms] <- 1-Ip_stay
  
  #' the rest of transitions are age specific - i.e.
  #' based on the fraction symptomatic for an age,
  #' Tij[Ianms[1], Enms] <- (1-E_stay)*(1-clinfrac)
  #' Tij[Ipnms[1], Enms] <- (1-E_stay)*clinfrac
  
  #' TODO age specific transmissions
  age_tij <- mapply(function(tij, symp_by_age) {
    tij[Ianms[1], Enms] <- (1-E_stay)*(1-symp_by_age)
    tij[Ipnms[1], Enms] <- (1-E_stay)*symp_by_age
    tij
  }, tij=list(Tij), pop$y, SIMPLIFY = FALSE)
  
  #' transitions
  #' make the contact matrices
  cm = Reduce(
    '+',
    mapply(
      function(c, m) c * m,
      redCon,
      pop$matrices,
      SIMPLIFY = F
    )
  )
  
  age_cats <- length(pop$size)
  
  #' alls age reference
  Mij <- Matrix(
    0,
    nrow = len_allstates*age_cats, ncol = len_allstates*age_cats,
    dimnames = list(
      to=outer(nms, sprintf("a%02i",1:16), paste, sep=""),
      from=outer(nms, sprintf("a%02i",1:16), paste, sep="")
    ),
    sparse = TRUE
  )
  
  bigTij <- Mij #' copy constructor
  
  for (age in 1:age_cats) {
    slc <- paste(nms, sprintf("a%02i",age), sep="")
    bigTij[slc,slc] <- age_tij[[age]]
  }
  
  toages <- paste(nms[1], sprintf("a%02i", 1:age_cats), sep="")
  
  for (fromage in 1:age_cats) {
    scont <- cm[, fromage] * umod
    Mij[toages, paste(Ianms, sprintf("a%02i", fromage), sep="")] <-
      scont * pop$fIa[fromage]
    Mij[toages, paste(Ipnms, sprintf("a%02i", fromage), sep="")] <-
      scont * pop$fIp[fromage]
    Mij[toages, paste(Isnms, sprintf("a%02i", fromage), sep="")] <-
      scont * redfIs[fromage]
  }
  
  totM <- Mij + bigTij
  RSpectra::eigs(totM, 10)
  
}

# sig <- function(sf) {
#   res <- matrix(0, 4, 4)
#   res[,1] <- c(0, 1-sf, sf, 0)
#   res[,3] <- c(0, 0, 0, 1)
#   res
# }
# sig2 <- function(sf) {
#   res <- matrix(0, 4, 4)
#   diag(res) <- -1
#   res[2:3,1] <- c(1-sf, sf)
#   res[4,3] <- 1
#   res
# }
# 
# trs <- function(foi) {
#   res <- matrix(0, 4, 4)
#   res[1,] <- c(0, foi, 2*foi, 2*foi)
#   res
# }

#' @examples 
#' thing <- ngm(params_back$pop[[1]])
#' thing2 <- ngm(params_back$pop[[1]], 1/thing$values)
#' thing2$values

saveRDS(ngm, tail(.args, 1))
