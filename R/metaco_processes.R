#' @title Metacommunity process: immigration
#' @description This function will calculate the propagule pressure. The effects of immigration are given as a weighted average of the occurrence probability of species i in neighborhood of z.
#' @param Y matrix of species occurrence
#' @param K connectivity matrix
#' @param m independent constant to account for immigration from outside the simulated metacommunity
#' @keywords immigration metacommunity
#'

I_f <- function(Y = "Species occurrence matrix",
                K = "Patch connectivity matrix",
                m = "outside immigration"){
  N <- nrow(Y)
  R <- ncol(Y)
  I <- (1-m) * (K %*% Y) / (K %*% matrix(1, nrow = N, ncol = R)) + m
  return(I)
}


#' @title Metacommunity process: patch connectivity
#' @description This function will calculate the connectivity matrix for a given set of patches and a dispersal parameter
#' @param XY corresponds to the patch coordinates
#' @param alpha dispersal parameter associated to an exponential distribution used for dispersal
#' @keywords patch


get_K <- function(XY, alpha){
  N <- nrow(XY)
  distMat <- as.matrix(dist(XY, method = "euclidean", upper = T, diag = T))
  ConMat <- exp((-1/alpha)*distMat)
  diag(ConMat) <- 0
  return(ConMat)
}


#' @title Metacommunity process: Environmental filtering with gaussian response
#' @description The effect of the environment on the local performance and colonization. Probability of establishing a viable local population given environmental conditions in the patch.
#' @param E matrix of environmental covariates for each patch
#' @param u_c considered as the niche optima for each species, for each environmental variable.
#' @param s_c understood as niche breadth for each species for each environmental variable.
#' @keywords environmental filtering gaussian response

S_f_gaussian <- function(E, u_c, s_c) {
  R <- ncol(u_c)
  N <- nrow(E)
  D <- ncol(E)
  S <- matrix(1, nr = N, nc = R)
  for(i in 1:D){
    optima <- matrix(u_c[i,],nrow = N,ncol = R,byrow = TRUE)
    breadth <- matrix(s_c[i,],nrow = N,ncol = R,byrow = TRUE)^2
    S <- S*exp(-(E[,i]-optima)^2 / breadth)
  }
  return(S)
}

#' @title Metacommunity process: Environmental filtering with quadratic response
#' @description The effect of the environment on the local performance and colonization. Probability of establishing a viable local population given environmental conditions in the patch.
#' @param E matrix of environmental covariates for each patch
#' @param u_c considered as the niche optima for each species, for each environmental variable.
#' @param s_c understood as niche breadth for each species for each environmental variable.
#' @keywords environmental filtering quadratic response
#'

S_f_quadratic <- function(E, u_c, s_c) {
  R <- ncol(u_c)
  N <- nrow(E)
  D <- ncol(E)
  S <- matrix(1, nr = N, nc = R)
  for(i in 1:D){
    optima <- matrix(u_c[i,],nrow = N,ncol = R,byrow = TRUE)
    breadth <- matrix(s_c[i,],nrow = N,ncol = R,byrow = TRUE)
    S <- S * ((-1 / (breadth/2)^2) * (E[,i] - optima)^2 + 1)
    S <- ifelse(S < 0, 0 , S)

  }
  return(S)
}

#' @title Metacommunity process: ecological interactions
#' @description Sum of ecological interactions. Computes the sum of ecological interactions for every location and every species
#' @param A species interaction matrix
#' @param Y species occurrence matrix
#' @keywords ecological interactions
#'

sum_interactions <- function(A, Y){
  t(A %*% t(Y))
}

#' @title Metacommunity process: colonization
#' @description The total effect of ecological interactions on the colonization probabilty
#' @param v sum of interactions
#' @param d_c sensitivity to interactions
#' @param c_0 colonization at zero interactions
#' @param c_max colonization at max interactions
#'

C_f <- function(v, d_c, c_0, c_max){

  c_max * (1 + ((1 / c_0) - 1) * exp(-v * d_c))^-1
  # Is this the same as:
  # c_max / (1 + ((1 / c_0) - 1) * exp(-v * d_c))
}



#' @title Metacomminbut process: extinction
#' @description Responses of the extinction probability to the local environment. Extinction rate should be minimal at environmental optimum
#' @param E Matrix of environmental covariates at each patch
#' @param u_e extinction rate at environmental optimum
#' @param s_e extintion rate at minimum
#'

M_f <- function(E, u_e, s_e) {
  R <- ncol(u_e)
  N <- nrow(E)
  D <- ncol(E)
  M <- matrix(1, nrow = N, ncol = R)
  for(i in 1:D){
    optima <- matrix(u_e[i,],nrow = N,ncol = R,byrow = TRUE)
    minima <- matrix(s_e[i,],nrow = N,ncol = R,byrow = TRUE)^2
    M <- M*(1 - exp(-(E[,i] - optima)^2 / minima))
  }
  return(M)
}


#' @title Metacommunity process: interactions and extinction
#' @description Effect of ecological interactions on species extinction at the patch level. In this case, extinction probability must be larger when interactions are negative and smaller when the interactions are positive.
#' @param v sum of ecological interactions
#' @param d_e sensitivity to interactions
#' @param e_0 extinction at zero interactions
#' @param e_min extinction at max interactions
#'

E_f <- function(v, d_e, e_0, e_min) {
  N <- nrow(v)
  R <- ncol(v)

  e_min_mat <- matrix(e_min, nrow = N, ncol = R, byrow=TRUE)

  e_min_mat + (1 / (1-e_min_mat) + (1/(e_0-e_min_mat) - 1 / (1 - e_min_mat)) * exp(d_e * v))^-1

}


