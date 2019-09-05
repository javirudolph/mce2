#' @title Metacommunity Simulation one set of parameters
#' @description This function will use the set of given parameters and run the main metacommunity simulation for a set amount of time steps. The output is a list of matrices, where each one of the matrices is an iteration of the metacommunity with the same parameters. It gives the option to save RDS files to the specified directory so these RDS files can be accessed later by other HMSC functions.
#' @param XY coordinates for each of the sites or patches
#' @param E matrix of environmental variables measured at each site
#' @param pars list of parameters given by the output of the prep_pars() function
#' @param nsteps numeric, number of time steps before getting the "snapshot" of a metacommunity
#' @param occupancy numeric value between 0-1, to set as the initial conditions occupancy of sites by the different species.
#' @param niter number of iterations for simulating the metacommunity with the same parameters.
#' @param envResp type of response to the environment: "gaussian" or "quadratic"
#' @param makeRDS should the function make an RDS file from the output and save it? Default is FALSE
#' @param whereToSave file path for the RDS file to be saved in
#' @param objName should the RDS file be saved, what should it be called?
#'
#'
metacom_sim4HMSC <- function(XY, E, pars, nsteps,
                             occupancy, niter,
                             envResp = "quadratic",
                             makeRDS = FALSE,
                             whereToSave = NULL,
                             objName = NULL){

  N <- pars$N
  D <- pars$D
  R <- pars$R

  Y0 <- ifelse(matrix(runif(N * D), nrow = N, ncol = R) < occupancy, 1, 0)

  res <- vector("list", length = niter)
  for(i in 1:niter){
    run <- mainfx(XY, E, pars, Y0, nsteps, envResp = envResp)
    res[[i]] <- run[[nsteps]]
  }

  if(makeRDS == TRUE){
    nameFile <- paste0(whereToSave, objName, "-metacomSim.RDS")
    saveRDS(res, file = nameFile)
  }

  return(res)

}


#' @title Metacommunity Simulation for multiple sets of parameters
#' @description This function will use the set of given parameters and run the main metacommunity simulation for a set amount of time steps. The output is a list of matrices, where each one of the matrices is an iteration of the metacommunity with the same parameters. It gives the option to save RDS files to the specified directory so these RDS files can be accessed later by other HMSC functions. It is different from the simulation with one set of parameters in that you can assign different dispersal and competition parameters to the species within your simulation.
#' @param XY coordinates for each of the sites or patches
#' @param E matrix of environmental variables measured at each site
#' @param pars in this case, the parameters are a list of prep_pars() outputs or the output of prep_multiparam()
#' @param nsteps numeric, number of time steps before getting the "snapshot" of a metacommunity
#' @param occupancy numeric value between 0-1, to set as the initial conditions occupancy of sites by the different species.
#' @param niter number of iterations for simulating the metacommunity with the same parameters.
#' @param envResp type of response to the environment: "gaussian" or "quadratic"
#' @param makeRDS should the function make an RDS file from the output and save it? Default is FALSE
#' @param whereToSave file path for the RDS file to be saved in
#' @param objName should the RDS file be saved, what should it be called?
#'
#'

# This function is for figure 3, in which we have separate groups of species with different dispersal or interaction parameters.
metacom_sim4HMSC_multParams <- function(XY, E, pars, nsteps,
                                        occupancy, niter,
                                        envResp = "quadratic",
                                        makeRDS = FALSE,
                                        whereToSave = NULL,
                                        objName = NULL){


  res <- vector("list", length = niter)
  for(i in 1:niter){
    bindruns <- NULL
    for(k in 1:length(pars)){
      subpars <- pars[[k]]

      N <- subpars$N
      D <- subpars$D
      R <- subpars$R

      Y0 <- ifelse(matrix(runif(N * D), nrow = N, ncol = R) < occupancy, 1, 0)
      run <- mainfx(XY, E, subpars, Y0, nsteps, envResp = envResp)
      lastrun <- run[[nsteps]]
      bindruns <- cbind(bindruns, lastrun)
    }
    res[[i]] <- bindruns
  }

  if(makeRDS == TRUE){
    nameFile <- paste0(whereToSave, objName, "-metacomSim.RDS")
    saveRDS(res, file = nameFile)
  }

  return(res)

}
