
#' @title Get parameters for Figure 2
#' @description  Manuscript specific. This function is for Figure 2 only.It takes the RDS params file and organizes it from a list to a dataframe, so that it can be joined to the Variation Partition results at the species level. That way we can know what niche values each species has. The output of this function can be joined to the VP results at the species level and then plot.
#' @param folderpath path to stored RDS files from the get_VP_results()
#' @param scenario the specific name of the file to use. Should be able to access the both the metacomSim, the model and the VP results with this name.
#'

get_fig2_params <- function(folderpath, scenario){
  pars <- readRDS(paste0(folderpath, scenario, "-params.RDS"))

  with(pars, {enframe(u_c[1,], name = "species", value = "nicheOptima") %>%
      left_join(., enframe(s_c[1,], name = "species", value = "nicheBreadth")) %>%
      left_join(., enframe(c_0, name = "species", value = "colonizationProb")) %>%
      mutate(., dispersal = alpha,
             species = as.character(species),
             speciesChr = as.character(paste0("spp_", species)),
             interCol = d_c ,
             interExt = d_e)})
}

#' @title Get parameters for Figure 3
#' @description Manuscript specific, figure 3 only.Figure 3 has three sets of parameters because we divided simulations into groups.Half of the species in each simulation has interactions, and the other half doesn't.Or, a third of the species is assigned a different dispersal level.
#' @param folderpath path to stored RDS files from the get_VP_results()
#' @param scenario the specific name of the file to use. Should be able to access the both the metacomSim, the model and the VP results with this name.
#'

get_fig3_params <- function(folderpath, scenario){
  parsList <- readRDS(paste0(folderpath, scenario, "-params.RDS"))
  fullPars <-  data.frame()
  for(i in 1:length(parsList)){

    nspp <- nrow(fullPars)

    i_pars <- with(parsList[[i]], {
      enframe(u_c[1,], name = "species",value =  "nicheOpt") %>%
        left_join(., enframe(s_c[1,], name = "species", value =  "nicheBreadth")) %>%
        left_join(., enframe(c_0, name = "species", value =  "colProb")) %>%
        mutate(dispersal = as.factor(alpha),
               species = as.character(species + nspp),
               intercol = as.factor(d_c),
               interext = d_e)})


    fullPars <- rbind(fullPars, i_pars)

  }

  return(fullPars)
}

#' @title Base Plots: Species
#' @description Base species plot
#' @param data formatted data from get_spp_data()
#' @param plotMain main title for the figure
#' @param colorVar variable for the color variation. Color by species or by niche optima
#' @param colorLegend should the color legend be included, what is the title?
#'
base_spp_plot <- function(data, plotMain = NULL, colorVar = NULL, colorLegend = "none"){
data %>%
  ggtern(aes(x = env, z = spa, y = codist, size = r2)) +
  geom_point(aes_string(color = colorVar), alpha = 0.8) +
  scale_T_continuous(limits=c(0.0,1.0),
                     breaks=seq(0.0,1.0,by=0.1),
                     labels=seq(0.0,1.0,by=0.1)) +
  scale_L_continuous(limits=c(0.0,1),
                     breaks=seq(0,1,by=0.1),
                     labels=seq(0,1,by=0.1)) +
  scale_R_continuous(limits=c(0.0,1.0),
                     breaks=seq(0,1,by=0.1),
                     labels=seq(0,1,by=0.1)) +
  labs(title = plotMain,
       x = "E",
       xarrow = "Environment",
       y = "C",
       yarrow = "Co-Distribution",
       z = "S",
       zarrow = "Spatial Autocorrelation") +
  theme_light() +
  theme_showarrows() +
  #scale_colour_brewer(palette = "Set1") +
  #scale_colour_brewer(palette = "Spectral") +
  #scale_color_viridis_d() +
  scale_size_area(limits = c(0,1), breaks = seq(0,1,0.2)) +
  guides(color = guide_legend(colorLegend, order = 2),
         size = guide_legend(title = expression(R^2), order = 1)) +
  theme(panel.grid = element_line(color = "darkgrey"),
        axis.text = element_text(size =5),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 12, margin = margin(t = 10, b = -20)),
        tern.axis.arrow = element_line(size = 1))
}

#' @title Base Plots: Sites
#' @description Base species plot
#' @param data formatted data from get_sites_data()
#' @param plotMain main title for the figure
#' @param colorVar variable for the color variation. Color by species or by niche optima
#' @param colorLegend should the color legend be included, what is the title?
#'

base_sites_plot <- function(data, plotMain = NULL, colorVar = NULL, colorLegend = "none"){
  data %>%
    ggtern(aes(x = env, z = spa, y = codist, size = r2)) +
    geom_point(aes_string(color = colorVar), alpha = 0.6) +
    scale_T_continuous(limits=c(0,1.0),
                       breaks=seq(0,1,by=0.1),
                       labels=seq(0,1,by=0.1)) +
    scale_L_continuous(limits=c(0.0,1),
                       breaks=seq(0,1,by=0.1),
                       labels=seq(0,1,by=0.1)) +
    scale_R_continuous(limits=c(0.0,1.0),
                       breaks=seq(0,1,by=0.1),
                       labels=seq(0,1,by=0.1)) +
    labs(title = plotMain,
         x = "E",
         xarrow = "Environment",
         y = "C",
         yarrow = "Co-Distribution",
         z = "S",
         zarrow = "Spatial Autocorrelation") +
    theme_light() +
    theme_showarrows() +
    #scale_colour_brewer(palette = "Set1") +
    #scale_colour_brewer(palette = "Spectral") +
    #scale_color_viridis_d() +
    scale_size_area(limits = c(0, 0.003), breaks = seq(0, 0.003, 0.0005)) +
    guides(color = guide_colorbar(colorLegend, order = 2),
           size = guide_legend(title = expression(R^2), order = 1)) +
    theme(panel.grid = element_line(color = "darkgrey"),
          axis.text = element_text(size =5),
          axis.title = element_text(size = 8),
          plot.title = element_text(size = 12, margin = margin(t = 10, b = -20)),
          tern.axis.arrow = element_line(size = 1))
}



