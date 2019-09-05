#' @title Organize and format site level data
#' @description Manuscript specific function. This function takes the simulated metacommunity to get the species richness for each site. It then binds that information to the Variation Partinioning output from the HMSC::variPart(). This function also takes the output of HMSC::variPart() and organizes each fraction: into the env, spa and codist components
#' @param folderpath path to stored RDS files from the get_VP_results()
#' @param scenario the specific name of the file to use. Should be able to access the both the metacomSim, the model and the VP results with this name.
#'
#'
get_sites_data <- function(folderpath, scenario){
  richness <- readRDS(paste0(folderpath, scenario, "-metacomSim.RDS")) %>%
    set_names(imap(., ~ paste0("iter", .y))) %>%
    map(., rowSums) %>%
    bind_rows() %>%
    rownames_to_column(var = "sites") %>%
    gather(., key = "iteration", value = "richness", -sites) %>%
    mutate(identifier = paste0("site", sites, "_", iteration)) %>%
    dplyr::select(., -c(sites, iteration))


  vp <- readRDS(paste0(folderpath, scenario, "-vpsites.RDS"))

  overlap1 <- map(vp, "overlap1")
  nspp <- as.numeric(dim(overlap1[[1]])[2])
  overlap2 <- map(vp, "overlap2")
  overlap3 <- map(vp, "overlap3")

  vpALL <- vector("list", length = 5)
  for(i in 1:5){
    workingVP1 <- overlap1[[i]]
    workingVP2 <- overlap2[[i]]
    workingVP3 <- overlap3[[i]]

    c <- rowSums(workingVP1[,,1])/nspp
    b <- rowSums(workingVP1[,,2])/nspp
    a <- rowSums(workingVP1[,,3])/nspp

    e <- rowSums(workingVP2[,,1])/nspp
    f <- rowSums(workingVP2[,,2])/nspp
    d <- rowSums(workingVP2[,,3])/nspp

    g <- rowSums(workingVP3)/nspp

    env <- a + f + 1/2 * d + 1/2 * g
    env <- ifelse(env < 0, 0, env)
    spa <- b + e + 1/2 * d + 1/2 * g
    spa <- ifelse(spa < 0, 0, spa)
    random <- c
    codist <- ifelse(random < 0, 0, random)
    r2 <- env + spa + codist
    iteration <- factor(paste0("iter", i), levels = paste0("iter", 1:5))

    cleanData <- cbind.data.frame(env, spa, codist, r2, iteration)
    cleanData$site <- paste0(row.names(cleanData))

    vpALL[[i]] <- cleanData
  }

  vpALL %>%
    bind_rows() %>%
    mutate(identifier = paste(site, iteration, sep = "_"),
           scenario = scenario) %>%
    left_join(., richness)
}

#' @title Organize and Format VP results at the species level
#' @description this function organizes species level data and gets it ready to plot. It takes the RDS file that is the output of the HMSC::variPart() function and reorganizes it into the env, codist and spa components.It also calculates the prevalence information for each species.
#' organizes each fraction: into the env, spa and codist components
#' @param folderpath path to stored RDS files from the get_VP_results()
#' @param scenario the specific name of the file to use. Should be able to access the both the metacomSim, the model and the VP results with this name.
#'

get_species_data <- function(folderpath, scenario){

  prevalence <- readRDS(paste0(folderpath, scenario, "-metacomSim.RDS")) %>%
    set_names(imap(., ~ paste0("iter_", .y))) %>%
    map(., colSums) %>%
    bind_cols() %>%
    rownames_to_column(var = "species") %>%
    gather(., key = "iteration", value = "prevalence", -species) %>%
    mutate(identifier = paste0("spp", species, "_", iteration)) %>%
    dplyr::select(., -c(species, iteration))




  readRDS(paste0(folderpath, scenario, "-vpspp.RDS")) %>%
    set_names(imap(., ~ paste0("iter_", .y))) -> VPdata

  fullData <- list()
  for(i in 1:length(VPdata)){
    fullData[[i]] <- VPdata[[i]] %>%
      map(as_tibble) %>%
      bind_cols() %>%
      rownames_to_column() %>%
      set_names(c("species", "c", "b", "a", "e", "f", "d", "g")) %>%
      transmute(species = species,
                env = a + f + 0.5 * d + 0.5 * g,
                env = ifelse(env < 0, 0, env),
                spa = b + e + 0.5 * d + 0.5 * g,
                spa = ifelse(spa < 0, 0, spa),
                codist = c,
                codist = ifelse(codist < 0, 0, codist),
                r2 = env + spa + codist,
                iteration = names(VPdata[i]))

  }

  fullData %>%
    bind_rows() %>%
    mutate(identifier = paste0("spp", species, "_", iteration),
           scenario = scenario) %>%
    left_join(., prevalence)

}
