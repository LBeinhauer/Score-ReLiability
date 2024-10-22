## Score ReLiability

# E1.5

# Clean & extract data


# library loading and installing as necessary


# relevant libraries required for this script
packages <- c("magrittr", "dplyr", "here", "patchwork", "ggplot2")

# check, whether library already installed or not - install and load as needed:
apply(as.matrix(packages), MARGIN = 1, FUN = function(x) {
  
  pkg_avail <- nzchar(system.file(package = x))   # check if library is installed on system
  
  if(pkg_avail){
    require(x, character.only = TRUE)             # load the library, if already installed
    
  }else{
    install.packages(x)                           # install the library, if missing
    require(x, character.only = TRUE)             # load after installation
  }
})

B_alpha_rma_df <- read.csv(here("Data/Processed/Reliability_analysis.csv"))

ES_rma_df <- read.csv(here("Data/Processed/Aggregates_ES_analysis.csv"))

agg_L <- readRDS(here("Data/Processed/Aggregates_simple.csv"))

MASC_names <- substr(list.files(here("Data/Extracted (Project) Data")), 1, nchar(list.files(here("Data/Extracted (Project) Data")))-4) 

# generate plots for meta-analyses of corrected AND uncorrected Cohen's d (slide 7)
plots2 <- lapply((1:length(agg_L)), FUN = function(idx){
  
  # name of MASC
  name <- MASC_names[idx]
  
  # get aggregates per sample
  y <- agg_L[[idx]]
  
  # get rma-results separated
  rma_raw <- ES_rma_df[ES_rma_df$MASC == name & ES_rma_df$corr == 0,]
  rma_corr <- ES_rma_df[ES_rma_df$MASC == name & ES_rma_df$corr == 1,]
  
  # get idea of minimum and maximum observations in MASC, to draw the respective distribution
  minmaxsequence_raw <- seq(from = min(y$d_raw), to = max(y$d_raw), length.out = 100)
  minmaxsequence_corr <- seq(from = min(y$d_corr), to = max(y$d_corr), length.out = 100)
  
  # Freedman-Diacons rule for binwidth:
  bw_raw <- (2 * IQR(y$d_raw) * (length(y$d_raw)^-(1/3)))
  bw_corr <- (2 * IQR(y$d_corr) * (length(y$d_corr)^-(1/3)))
  
  # make plot
  p <- ggplot() +
    scale_shape_identity() +
    # geom_point(aes(x = y$d_raw, y = 0), colour = "darkgrey", shape = 108, size = 3) +
    # geom_point(aes(x = y$d_corr, y = 0), colour = "black", shape = 108, size = 3) +
    # geom_histogram(aes(x = y$d_raw), colour = NULL, fill = "darkgrey", binwidth = bw_raw, alpha = .5) +
    # geom_histogram(aes(x = y$d_corr), colour = NULL, fill = "black", binwidth = bw_raw, alpha = .5) +
    # geom_point(aes(x = rma_raw$mu, y = 0), colour = "darkgrey", shape = 108, size = 8) +
    # geom_point(aes(x = rma_corr$mu, y = 0), colour = "black", shape = 108, size = 8) +
    geom_line(aes(x = minmaxsequence_raw,
                  y = dnorm(minmaxsequence_raw, mean = rma_raw$mu, sd = sqrt(rma_raw$tau^2))),
              linewidth = 1, colour = "darkgrey") +
    geom_line(aes(x = minmaxsequence_corr,
                  y = dnorm(minmaxsequence_corr, mean = rma_corr$mu, sd = sqrt(rma_corr$tau^2))),
              linewidth = 1, colour = "black") +
    theme(axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", colour = "transparent"), 
          panel.grid.major.y = element_line(colour = "transparent"),
          panel.grid.major.x = element_line(colour = "grey"),
          panel.grid.minor = element_line(colour = "transparent"),
          axis.ticks.x = element_line(colour = "grey"),
          strip.background = element_rect(fill = "transparent")) +
    labs(subtitle = name)
  
  # this is actually an artefact, we don't need this, heterogeneity is never truly zero in this selection
  p
  
})

# combine plots using pathwork-infrastructure
combined_plot2 <- Reduce('+', plots2) + plot_layout(ncol = 2) &
  plot_annotation(theme = theme(plot.background = element_rect(fill = "transparent", colour = "transparent")))
combined_plot2

# save plot
ggsave(filename = here("Graphics/densities_four.png"),
       plot = last_plot(),
       width = 11, height = 5)



forest_plot_rel <- function(rma_df, aggregates, ci.lvl = .975){
  rma_raw <- rma_df[rma_df$corr == 0,]
  rma_corr <- rma_df[rma_df$corr == 1,]
  
  aggregates %>% 
    mutate(cil_raw = d_raw - (1-(1-ci.lvl)/2) * (SE_d),
           ciu_raw = d_raw + (1-(1-ci.lvl)/2) * (SE_d),
           cil_corr = d_corr - (1-(1-ci.lvl)/2) * (SE_d_corr),
           ciu_corr = d_corr + (1-(1-ci.lvl)/2) * (SE_d_corr)) %>% 
    arrange(desc(d_raw)) %>% 
    ggplot() +
    # point estimate of raw d
    geom_point(aes(x = d_raw, y = 1:nrow(aggregates)), colour = "darkgrey") +
    # CI of point estimate of raw d
    geom_segment(aes(x = cil_raw, y = 1:nrow(aggregates), xend = ciu_raw, yend = 1:nrow(aggregates)), colour = "darkgrey") +
    geom_segment(aes(x = cil_raw, xend = cil_raw, y = (1:nrow(aggregates))+.3, yend = (1:nrow(aggregates))-.3), colour = "darkgrey") +
    geom_segment(aes(x = ciu_raw, xend = ciu_raw, y =( 1:nrow(aggregates))+.3, yend = (1:nrow(aggregates))-.3), colour = "darkgrey") +
    # point estiamte of corrected d
    geom_point(aes(x = d_corr, y = 1:nrow(aggregates)), colour = "black") +
    # CI of point estimate of corrected d
    geom_segment(aes(x = cil_corr, y = 1:nrow(aggregates), xend = ciu_corr, yend = 1:nrow(aggregates)), colour = "black") +
    geom_segment(aes(x = cil_corr, xend = cil_corr, y = (1:nrow(aggregates))+.3, yend = (1:nrow(aggregates))-.3), colour = "black") +
    geom_segment(aes(x = ciu_corr, xend = ciu_corr, y = (1:nrow(aggregates))+.3, yend = (1:nrow(aggregates))-.3), colour = "black") +
    # solid black line separating RMA-estimates from point estimates
    geom_abline(slope = 0, intercept = -1, colour = "black") +
    # adding diamond showing rma-estimate of raw d
    geom_polygon(data = data.frame(x = c(rma_raw$ci.ub, rma_raw$mu, rma_raw$ci.lb, rma_raw$mu),
                                   y = c(-3, -3-.7, -3, -3+.7)),
                 aes(x = x, y = y),
                 fill = "darkgrey", colour = "darkgrey") +
    # adding diamond showing rma-estimate of corrected d
    geom_polygon(data = data.frame(x = c(rma_corr$ci.ub, rma_corr$mu, rma_corr$ci.lb, rma_corr$mu),
                                   y = c(-3, -3-.7, -3, -3+.7)),
                 aes(x = x, y = y),
                 colour = "black", fill = "black") +
    # defining theme of plot (transparent background etc.)
    theme(legend.position = "bottom", 
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", colour = "transparent"), 
          panel.grid.major.y = element_line(colour = "transparent"),
          panel.grid.major.x = element_line(colour = "grey"),
          panel.grid.minor = element_line(colour = "transparent"),
          axis.ticks = element_line(colour = "grey"),
          strip.background = element_rect(fill = "transparent"),
          strip.text = element_text(size = 12)) +
    scale_y_continuous(breaks = c(-3, 1:nrow(aggregates)), labels = c("RMA-estimate", aggregates$source)) +
    labs(x = "Cohen's d",
         y = "Country") 
}

forest_plot_rel(ES_rma_df[c(1,2),], agg_L[[1]])



plots3 <- lapply((1:length(agg_L)), FUN = function(idx){
  
  # name of MASC
  name <- MASC_names[idx]
  
  # get aggregates per sample
  y <- agg_L[[idx]]
  
  # get rma-results separated
  rma_raw <- ES_rma_df[ES_rma_df$MASC == name & ES_rma_df$corr == 0,]
  rma_corr <- ES_rma_df[ES_rma_df$MASC == name & ES_rma_df$corr == 1,]
  
  forest_plot_rel(ES_rma_df[ES_rma_df$MASC == name,], y)
  
})

# combine plots using pathwork-infrastructure
combined_plot3 <- Reduce('+', plots3) + plot_layout(ncol = 2) &
  plot_annotation(theme = theme(plot.background = element_rect(fill = "transparent", colour = "transparent")))
combined_plot3

