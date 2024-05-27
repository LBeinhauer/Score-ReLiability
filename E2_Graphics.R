## Score ReLiability

# E2

# Graphics



# library loading and installing as necessary


# relevant libraries required for this script
packages <- c("tidyverse", "here", "metafor", "ggpubr", "gridExtra")

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


##### PHOTOGRAPHS

# Read data-files containing estimates of SMDs and Standard Errors
pos_SMDs.df <- read.csv(file = here("Data/Processed/SMDs_positive.csv"))
neg_SMDs.df <- read.csv(file = here("Data/Processed/SMDs_negative.csv"))


# fit a meta-analytic model to estimates of uncorrected Cohen's d, for positiveness-ratings
pos_SMDs_rma_raw <- metafor::rma(yi = d_raw, 
                                 sei = SE_d_raw, 
                                 data = pos_SMDs.df,
                                 measure = "GEN",
                                 method = "REML")

# fit a meta-analytic model to estimates of corrected Cohen's d, for positiveness-ratings
pos_SMDs_rma_cor <- metafor::rma(yi = d_cor, 
                                 sei = SE_d_cor, 
                                 data = pos_SMDs.df,
                                 measure = "GEN",
                                 method = "REML")


# fit a meta-analytic model to estimates of uncorrected Cohen's d, for negativeness-ratings
neg_SMDs_rma_raw <- metafor::rma(yi = d_raw, 
                                 sei = SE_d_raw, 
                                 data = neg_SMDs.df,
                                 measure = "GEN",
                                 method = "REML")


# fit a meta-analytic model to estimates of corrected Cohen's d, for negativeness-ratings
neg_SMDs_rma_cor <- metafor::rma(yi = d_cor, 
                                 sei = SE_d_cor, 
                                 data = neg_SMDs.df,
                                 measure = "GEN",
                                 method = "REML")



##### STATE EMOTIONS

# fit a meta-analytic model to estimates of uncorrected Cohen's d, for positiveness-ratings
pos_state_SMDs_rma_raw <- metafor::rma(yi = d_raw, 
                                       sei = SE_d_raw, 
                                       data = pos_state_SMDs.df,
                                       measure = "GEN",
                                       method = "REML")

# fit a meta-analytic model to estimates of corrected Cohen's d, for positiveness-ratings
pos_state_SMDs_rma_cor <- metafor::rma(yi = d_cor, 
                                       sei = SE_d_cor, 
                                       data = pos_state_SMDs.df,
                                       measure = "GEN",
                                       method = "REML")




# fit a meta-analytic model to estimates of uncorrected Cohen's d, for negativeness-ratings
neg_state_SMDs_rma_raw <- metafor::rma(yi = d_raw, 
                                       sei = SE_d_raw, 
                                       data = neg_state_SMDs.df,
                                       measure = "GEN",
                                       method = "REML")


# fit a meta-analytic model to estimates of corrected Cohen's d, for negativeness-ratings
neg_state_SMDs_rma_cor <- metafor::rma(yi = d_cor, 
                                       sei = SE_d_cor, 
                                       data = neg_state_SMDs.df,
                                       measure = "GEN",
                                       method = "REML")




##### ANTICIPATED EMOTIONS


# fit a meta-analytic model to estimates of uncorrected Cohen's d, for positiveness-ratings
pos_antc_SMDs_rma_raw <- metafor::rma(yi = d_raw, 
                                      sei = SE_d_raw, 
                                      data = pos_antc_SMDs.df,
                                      measure = "GEN",
                                      method = "REML")

# fit a meta-analytic model to estimates of corrected Cohen's d, for positiveness-ratings
pos_antc_SMDs_rma_cor <- metafor::rma(yi = d_cor, 
                                      sei = SE_d_cor, 
                                      data = pos_antc_SMDs.df,
                                      measure = "GEN",
                                      method = "REML")





# fit a meta-analytic model to estimates of uncorrected Cohen's d, for negativeness-ratings
neg_antc_SMDs_rma_raw <- metafor::rma(yi = d_raw, 
                                      sei = SE_d_raw, 
                                      data = neg_antc_SMDs.df,
                                      measure = "GEN",
                                      method = "REML")


# fit a meta-analytic model to estimates of corrected Cohen's d, for negativeness-ratings
neg_antc_SMDs_rma_cor <- metafor::rma(yi = d_cor, 
                                      sei = SE_d_cor, 
                                      data = neg_antc_SMDs.df,
                                      measure = "GEN",
                                      method = "REML")




forest_plot_rel <- function(rma.fit_raw, rma.fit_cor, rma.data, ci.lvl = .975){
  rma.data %>% 
    mutate(cil_raw = d_raw - (1-(1-ci.lvl)/2) * (SE_d_raw),
           ciu_raw = d_raw + (1-(1-ci.lvl)/2) * (SE_d_raw),
           cil_cor = d_cor - (1-(1-ci.lvl)/2) * (SE_d_cor),
           ciu_cor = d_cor + (1-(1-ci.lvl)/2) * (SE_d_cor)) %>% 
    arrange(desc(d_raw)) %>% 
    ggplot() +
    # point estimate of raw d
    geom_point(aes(x = d_raw, y = 1:nrow(rma.data)), colour = "darkgrey") +
    # CI of point estimate of raw d
    geom_segment(aes(x = cil_raw, y = 1:nrow(rma.data), xend = ciu_raw, yend = 1:nrow(rma.data)), colour = "darkgrey") +
    geom_segment(aes(x = cil_raw, xend = cil_raw, y = (1:nrow(rma.data))+.3, yend = (1:nrow(rma.data))-.3), colour = "darkgrey") +
    geom_segment(aes(x = ciu_raw, xend = ciu_raw, y =( 1:nrow(rma.data))+.3, yend = (1:nrow(rma.data))-.3), colour = "darkgrey") +
    # point estiamte of corrected d
    geom_point(aes(x = d_cor, y = 1:nrow(rma.data)), colour = "black") +
    # CI of point estimate of corrected d
    geom_segment(aes(x = cil_cor, y = 1:nrow(rma.data), xend = ciu_cor, yend = 1:nrow(rma.data)), colour = "black") +
    geom_segment(aes(x = cil_cor, xend = cil_cor, y = (1:nrow(rma.data))+.3, yend = (1:nrow(rma.data))-.3), colour = "black") +
    geom_segment(aes(x = ciu_cor, xend = ciu_cor, y = (1:nrow(rma.data))+.3, yend = (1:nrow(rma.data))-.3), colour = "black") +
    # solid black line separating RMA-estimates from point estimates
    geom_abline(slope = 0, intercept = -1, colour = "black") +
    # adding diamond showing rma-estimate of raw d
    geom_polygon(data = data.frame(x = c(rma.fit_raw$ci.ub, rma.fit_raw$b[1], rma.fit_raw$ci.lb, rma.fit_raw$b[1]),
                                   y = c(-3, -3-.7, -3, -3+.7)),
                 aes(x = x, y = y),
                 fill = "darkgrey", colour = "darkgrey") +
    # adding diamond showing rma-estimate of corrected d
    geom_polygon(data = data.frame(x = c(rma.fit_cor$ci.ub, rma.fit_cor$b[1], rma.fit_cor$ci.lb, rma.fit_cor$b[1]),
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
    scale_y_continuous(breaks = c(-3, 1:nrow(rma.data)), labels = c("RMA-estimate", rma.data$country)) +
    labs(x = "Cohen's d",
         y = "Country") 
}

forest_plot_rel(pos_SMDs_rma_raw, pos_SMDs_rma_cor, pos_SMDs.df)

ggsave(filename = "C:/Users/Lukas/Downloads/Forest_PCAR002_pos.png",
       plot = last_plot(),
       width = 8,
       height = 4.5)

forest_plot_rel(neg_SMDs_rma_raw, neg_SMDs_rma_cor, neg_SMDs.df)



forest_plot_rel(pos_state_SMDs_rma_raw, pos_state_SMDs_rma_cor, pos_state_SMDs.df)
forest_plot_rel(neg_state_SMDs_rma_raw, neg_state_SMDs_rma_cor, neg_state_SMDs.df)


forest_plot_rel(pos_antc_SMDs_rma_raw, pos_antc_SMDs_rma_cor, pos_antc_SMDs.df)
forest_plot_rel(neg_antc_SMDs_rma_raw, neg_antc_SMDs_rma_cor, neg_antc_SMDs.df)




underlying_density <- function(rma.fit_raw, rma.fit_cor){
  
  x_raw <- seq(rma.fit_raw$b[1] - 3*sqrt(rma.fit_raw$tau2), rma.fit_raw$b[1] + 3*sqrt(rma.fit_raw$tau2), length.out = 1000)
  density_values_raw <- dnorm(x_raw, mean = rma.fit_raw$b[1], sd = sqrt(rma.fit_raw$tau2)) 
  df_raw <- data.frame(x_raw, density_values_raw)
  
  x_cor <- seq(rma.fit_cor$b[1] - 3*sqrt(rma.fit_cor$tau2), rma.fit_cor$b[1] + 3*sqrt(rma.fit_cor$tau2), length.out = 1000)
  density_values_cor <- dnorm(x_cor, mean = rma.fit_cor$b[1], sd = sqrt(rma.fit_cor$tau2)) 
  df_cor <- data.frame(x_cor, density_values_cor)
  
  df <- data.frame(df_raw, df_cor)
  
  
  ggplot(df) +
    geom_line(aes(x = x_raw, y = density_values_raw), colour = "darkgrey") +
    geom_area(aes(x = x_raw, y = density_values_raw), fill = "darkgrey") +
    geom_line(aes(x = x_cor, y = density_values_cor), colour = "black", alpha = .7) +
    geom_area(aes(x = x_cor, y = density_values_cor), fill = "black", alpha = .7) +
    labs(x = "Cohen's d", 
         y = "Density") +
    theme(legend.position = "none", 
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", colour = "transparent"), 
          panel.grid.major.y = element_line(colour = "transparent"),
          panel.grid.major.x = element_line(colour = "black"),
          panel.grid.minor = element_line(colour = "transparent"),
          axis.ticks = element_line(colour = "black"))
  
  
}



##### PHOTOGRPAPHS

underlying_density(pos_SMDs_rma_raw, pos_SMDs_rma_cor)

ggsave(filename = "C:/Users/Lukas/Downloads/Density_PCAR002_pos.png",
       plot = last_plot(),
       width = 8,
       height = 4.5)

underlying_density(neg_SMDs_rma_raw, neg_SMDs_rma_cor)


##### STATE EMOTIONS

underlying_density(pos_state_SMDs_rma_raw, pos_state_SMDs_rma_cor)
underlying_density(neg_state_SMDs_rma_raw, neg_state_SMDs_rma_cor)


##### ANTICIPATED EMOTIONS

# throws error as no heterogeneity
# underlying_density(pos_antc_SMDs_rma_raw, pos_antc_SMDs_rma_cor)
underlying_density(neg_antc_SMDs_rma_raw, neg_antc_SMDs_rma_cor)




fig1.p <- gridExtra::arrangeGrob(forest_plot_rel(pos_SMDs_rma_raw, pos_SMDs_rma_cor, pos_SMDs.df) +
                                   labs(subtitle = "a)"),
                                 underlying_density(pos_SMDs_rma_raw, pos_SMDs_rma_cor) +
                                   labs(subtitle = "b)"),
                                 layout_matrix = matrix(c(rep(1, 90), rep(2, 54)), byrow = FALSE, ncol = 16, nrow = 9))

ggsave(filename = "C:/Users/Lukas/Downloads/ScoreReLiability_fig1_pos.png",
       plot = fig1.p,
       width = 8,
       height = 4.5)
