#' @title Format metacommunity data to "HMSC" classes
#' @description This function will use the simulated metacommunity matrices and format the data into "HMSC" classes
#' @param metacomData list of simulated metacommunity matrices. Output from metacom_sim4HMSC() or metacom_sim4HMSC_multiparam()
#' @param numClusters number of clusters based on number of cores available to work with DoParallel
#' @param E matrix of environmental variables associated to each site
#' @param MEMsel distance-based Moran' Eigenvector Map's eigenvector maps from a geographic distance matrix (xy matrix). This comes from adespatial::dbmem()
#' @param hmscPars list of parameters to be included for iterations, burn and thin components of the as.HMSCdata()
#' @param makeRDS should the function make an RDS file from the output and save it? Default is FALSE
#' @param whereToSave file path for the RDS file to be saved in
#' @param objName should the RDS file be saved, what should it be called?
#'
metacom_as_HMSCdata <- function(metacomData, numClusters, E, MEMsel,
                                hmscPars = NULL,
                                makeRDS = FALSE,
                                whereToSave = NULL,
                                objName = NULL){

  N <- nrow(E)

  run <- metacomData
  nrun <- length(run)

  clusters <- makeCluster(numClusters)
  registerDoParallel(clusters)

  if(is.null(hmscPars) == TRUE){
    hmscniter <- 100000
    hmscnburn <- 10000
    hmscthin <- 20
  } else{
    hmscniter <- hmscPars$niter
    hmscnburn <- hmscPars$nburn
    hmscthin <- hmscPars$thin
  }

  ### Estimate models
  model <- foreach(j = 1:nrun) %dopar% {
    library(HMSC)
    formData <- as.HMSCdata(Y = run[[j]], X = cbind(scale(E),scale(E)^2, MEMsel),
                            Random = as.factor(1:N),
                            scaleX = TRUE, interceptX = TRUE)

    hmsc(formData, family = "probit",
         niter = hmscniter, nburn = hmscnburn, thin = hmscthin)
  }

  ### Stop clusters
  stopCluster(clusters)

  if(makeRDS == TRUE){
    nameFile <- paste0(whereToSave, objName, "-model.RDS")
    saveRDS(model, file = nameFile)
  }

  return(model)
}
