#' @title Create the species correlation matrices
#' @description This mostly follows the code provided from the HMSC vignette, therefore it is not in a tidyverse framework and requires the corrplot package
#' @param folderpath path to stored RDS files from the get_VP_results()
#' @param scenario the specific name of the file to use. Should be able to access the HMSC model.
#' @param iteration Which iteration of the model should it use? Default is 1
#' @param corTitle Title for the correlation plot


interaction_plot <- function(folderpath, scenario, iteration = NULL, corTitle = NULL){

  if(is.null(iteration) == TRUE){
    iteration <- 1
  }


  modelfile <- readRDS(paste0(folderpath, scenario, "-model.RDS"))
  assoMat <- corRandomEff(modelfile[[iteration]])
  siteMean <- apply(assoMat[ , , , 1], 1:2, mean)

  siteDrawCol <- matrix(NA, nrow = nrow(siteMean), ncol = ncol(siteMean))
  siteDrawCol[which(siteMean > 0.4, arr.ind=TRUE)]<-"red"
  siteDrawCol[which(siteMean < -0.4, arr.ind=TRUE)]<-"blue"

  # Build matrix of "significance" for corrplot
  siteDraw <- siteDrawCol
  siteDraw[which(!is.na(siteDraw), arr.ind = TRUE)] <- 0
  siteDraw[which(is.na(siteDraw), arr.ind = TRUE)] <- 1
  siteDraw <- matrix(as.numeric(siteDraw), nrow = nrow(siteMean), ncol = ncol(siteMean))

  Colour <- colorRampPalette(c("blue", "white", "red"))(200)
  corrplot(siteMean, method = "color", col = Colour, type = "lower",
           diag = FALSE, p.mat = siteDraw, tl.srt = 45, title = corTitle)

}
