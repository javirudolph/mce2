#' @title Variation partinioning: species
#' @desctiption use of the VariPart function from HMSC, gets results at the level of species
#' @param HMSCmodel output from the metacom_as_HMSCdata() or an object from as.HMSCdata()
#' @param MEMsel distance-based Moran' Eigenvector Map's eigenvector maps from a geographic distance matrix (xy matrix). This comes from adespatial::dbmem()
#' @param numClusters number of clusters based on number of cores available to work with DoParallel
#' @param makeRDS should the function make an RDS file from the output and save it? Default is FALSE
#' @param whereToSave file path for the RDS file to be saved in
#' @param objName should the RDS file be saved, what should it be called?
#'
get_VPresults <- function(HMSCmodel, MEMsel, numClusters,
                          makeRDS = FALSE,
                          whereToSave = NULL,
                          objName = NULL){

  model <- HMSCmodel
  nmodel <- length(model)

  clusters <- makeCluster(numClusters)
  registerDoParallel(clusters)

  ### Estimate models
  vpRes <- foreach(j = 1:nmodel) %dopar% {
    library(HMSC)
    variPart(model[[j]], groupX = c(rep("env",3),rep("spa",length(MEMsel))),
             type = "III", R2adjust = TRUE)
  }

  ### Stop clusters
  stopCluster(clusters)

  if(makeRDS == TRUE){
    nameFile <- paste0(whereToSave, objName, "-vpspp.RDS")
    saveRDS(vpRes, file = nameFile)
  }

  return(vpRes)

}

# This function is necessary when running infor on sites.
get_VPresults_SITE <- function(HMSCmodel, MEMsel, numClusters,
                               makeRDS = FALSE,
                               whereToSave = NULL,
                               objName = NULL){

  model <- HMSCmodel
  nmodel <- length(model)

  clusters <- makeCluster(numClusters)
  registerDoParallel(clusters)

  ### Estimate models
  vpRes <- foreach(j = 1:nmodel) %dopar% {
    library(HMSC)
    variPart(model[[j]], groupX = c(rep("env",3),rep("spa",length(MEMsel))),
             indSite = TRUE,
             type = "III", R2adjust = TRUE)
  }

  ### Stop clusters
  stopCluster(clusters)

  if(makeRDS == TRUE){
    nameFile <- paste0(whereToSave, objName, "-vpsites.RDS")
    saveRDS(vpRes, file = nameFile)
  }

  return(vpRes)
}
