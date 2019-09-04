#' @title Prepare metacommunity parameters
#' @description This function is used to prepare the parameters in a useful format and avoid lengthy code for every simulation. It will output all parameters in a list ready for the main_sim_fx()
#' @param N number of patches for the simulation
#' @param D number of environmental variables simulated
#' @param R number of species for the simulation
#' @param breadth niche breadth. As reference 0.8 as narrow and 2 as broad with a quadratic response.
#' @param alpha dispersal parameter associated to an exponential distribution
#' @param interex Species sensitivity to interactions. For our simulations, sensitivity is set the same for all species, and both for colonization and extinction
#' @param initocup initial occupancy. For our simulations set to 0.8
#'

prep_pars <- function(N = NULL,
                      D = NULL,
                      R = NULL,
                      nicheOpt = NULL,
                      breadth = NULL,
                      alpha = NULL,
                      interx_col = NULL,
                      interx_ext = NULL,
                      makeRDS = FALSE,
                      whereToSave = NULL,
                      objName = NULL){

  # Number of patches, environmental variables and species
  N <- N
  D <- D
  R <- R

  u_c <- nicheOpt

  # Effect of environment on colonization
  if(is.null(nicheOpt) == TRUE){
    u_c <- matrix(nrow = D, ncol = R)
    u_c[1,] <- seq(0.1, 0.9, length = R)
  }

  s_c <- matrix(breadth, nrow = D, ncol = R)

  # Effect of the environment on extinction
  u_e <- matrix(nrow = D, ncol = R)
  u_e[1,] <- rep(0.5, R)
  s_e <- matrix(Inf, nrow = D, ncol = R)

  # Mean dispersal
  alpha <- alpha

  # Immigration
  m <- 0.001

  # Colonication function
  c_0 <- rep(0.4, R) # Colonization at 0 interactions
  c_max <- rep(1, R) # Colonization at max interactions


  # Extinction function
  e_0 <- rep(0.025, R) # Extinction at 0 interactions
  e_min <- rep(0, R) # Extinction at max interactions

  # Sensitivity to interactions
  d_c <- interx_col
  d_e <- interx_ext

  # Interaction matrix
  A <- matrix(0, ncol = R, nrow = R)
  d <- as.matrix(dist(c(1:R), upper = TRUE, diag = TRUE))
  A[d <= 1] <- -1
  diag(A) <- 0


  pars <- list(N = N, D = D, R = R,
               u_c = u_c, u_e = u_e, s_c = s_c, s_e = s_e, alpha = alpha, m = m,
               c_0 = c_0, e_0 = e_0, c_max = c_max, e_min = e_min, d_c = d_c, d_e = d_e, A = A)

  if(makeRDS == TRUE){
    nameFile <- paste0(whereToSave, objName, "-params.RDS")
    saveRDS(pars, file = nameFile)
  }

  return(pars)
}
