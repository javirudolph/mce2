#' @title Metacommunity process: immigration
#' @description This function will calculate the propagule pressure. The effects of immigration are given as a weighted average of the occurrence probability of species i in neighborhood of z.
#' @param Y matrix of species occurrence
#' @param K connectivity matrix
#' @param m independent constant to account for immigration from outside the simulated metacommunity
#' @keywords immigration metacommunity
#' @examples
#' I_f(Y, K, m)
#'

I_f <- function(Y = "Species occurrence matrix",
                K = "Patch connectivity matrix",
                m = "outside immigration"){
  N <- nrow(Y)
  R <- ncol(Y)
  I <- (1-m) * (K %*% Y) / (K %*% matrix(1, nrow = N, ncol = R)) + m
  return(I)
}
