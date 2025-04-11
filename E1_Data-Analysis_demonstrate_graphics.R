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

source(here("ReLiability_Function-library.R"))


B_alpha_rma_df <- read.csv(here("Data/Processed/Reliability_analysis.csv"))

ES_rma_df <- read.csv(here("Data/Processed/Aggregates_ES_analysis.csv"))

agg_L <- readRDS(here("Data/Processed/Aggregates_simple.csv"))

MASC_names <- unique(ES_rma_df$MASC)


# generate plots for meta-analyses of corrected AND uncorrected Cohen's d (slide 7)
plots2 <- lapply((1:length(agg_L)), FUN = function(idx){
  
  # name of MASC
  name <- MASC_names[idx]
  
  # get aggregates per sample
  y <- na.omit(agg_L[[idx]])
  
  # get rma-results separated
  rma_raw <- ES_rma_df[ES_rma_df$MASC == name & ES_rma_df$corr == 0,]
  rma_corr <- ES_rma_df[ES_rma_df$MASC == name & ES_rma_df$corr == 1,]
  
  # # get idea of minimum and maximum observations in MASC, to draw the respective distribution
  minmaxsequence_raw <- seq(from = min(y$d_raw), to = max(y$d_raw), length.out = 100)
  minmaxsequence_corr <- seq(from = min(y$d_corr), to = max(y$d_corr), length.out = 100)
  # 
  # # Freedman-Diacons rule for binwidth:
  # bw_raw <- (2 * IQR(y$d_raw) * (length(y$d_raw)^-(1/3)))
  # bw_corr <- (2 * IQR(y$d_corr) * (length(y$d_corr)^-(1/3)))
  
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
combined_plot2 <- Reduce('+', plots2) + plot_layout(ncol = 5) &
  plot_annotation(theme = theme(plot.background = element_rect(fill = "transparent", colour = "transparent")))
combined_plot2

# save plot
ggsave(filename = here("Graphics/densities_twentytwo.png"),
       plot = last_plot(),
       width = 11, height = 5)







########################################
# Analysis using informative data-sets #
########################################


# Identify informative phenomena

## Criterion 3: non-null effects


nn_eff_idx <- which(ES_rma_df$pval[ES_rma_df$corr == 0] <= .05)

# generate plots for meta-analyses of corrected AND uncorrected Cohen's d (slide 7)
plots3 <- lapply((1:length(agg_L))[nn_eff_idx], FUN = function(idx){
  
  # name of MASC
  name <- MASC_names[idx]
  
  # get aggregates per sample
  y <- na.omit(agg_L[[idx]])
  
  # get rma-results separated
  rma_raw <- ES_rma_df[ES_rma_df$MASC == name & ES_rma_df$corr == 0,]
  rma_corr <- ES_rma_df[ES_rma_df$MASC == name & ES_rma_df$corr == 1,]
  
  # # get idea of minimum and maximum observations in MASC, to draw the respective distribution
  minmaxsequence_raw <- seq(from = min(y$d_raw), to = max(y$d_raw), length.out = 100)
  minmaxsequence_corr <- seq(from = min(y$d_corr), to = max(y$d_corr), length.out = 100)
  # 
  # # Freedman-Diacons rule for binwidth:
  # bw_raw <- (2 * IQR(y$d_raw) * (length(y$d_raw)^-(1/3)))
  # bw_corr <- (2 * IQR(y$d_corr) * (length(y$d_corr)^-(1/3)))
  
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
combined_plot3 <- Reduce('+', plots3) + plot_layout(ncol = 5) &
  plot_annotation(theme = theme(plot.background = element_rect(fill = "transparent", colour = "transparent")))
combined_plot3



# res_tab <- ES_rma_df[ES_rma_df$corr == 0,][nn_eff_idx,] %>% 
#   select(MASC, mu, se, pval, tau, QE, k, QEp) %>% 
#   mutate(mu_str = paste0(round(mu, 3), " (", round(se, 3), ")"),
#          QE_str = paste0(round(QE, 3), " (", k-1, ")"))
# 
# MD_tab <- knitr::kable(res_tab %>% 
#                select(MASC, mu_str, pval, tau, QE_str, QEp), 
#              digits = 3, 
#              col.names = c("MASC", "$\\mu_{ES}$ (SE)", "p", "$\\tau$", "$Q_E$ (df)", "p"))
# 
# writeLines(MD_tab, "C:/Users/Lukas/Downloads/markdown_table.txt")
# write.csv(res_tab %>% 
#             mutate_if(is.numeric, round, 3) %>% 
#             mutate(pval = ifelse(pval < .001, yes = "<.001", no = pval),
#                    QEp = ifelse(pval < .001, yes = "<.001", no = QEp)) %>% 
#             select(MASC, mu_str, pval, tau, QE_str, QEp), 
#           "C:/Users/Lukas/Downloads/markdown_table.csv",
#           row.names = F)

res_tab <- data.frame(MASC = ES_rma_df$MASC[which(ES_rma_df$corr == 0)],
                      k = ES_rma_df$k[which(ES_rma_df$corr == 0)],
                      mu_raw = ES_rma_df$mu[which(ES_rma_df$corr == 0)],
                      mu_cor = ES_rma_df$mu[which(ES_rma_df$corr == 1)],
                      mu_SE_raw = ES_rma_df$se[which(ES_rma_df$corr == 0)],
                      mu_SE_cor = ES_rma_df$se[which(ES_rma_df$corr == 1)],
                      tau_raw = ES_rma_df$tau[which(ES_rma_df$corr == 0)],
                      tau_cor = ES_rma_df$tau[which(ES_rma_df$corr == 1)],
                      QEp_raw = ES_rma_df$QEp[which(ES_rma_df$corr == 0)],
                      QEp_cor = ES_rma_df$QEp[which(ES_rma_df$corr == 1)],
                      QE_raw = ES_rma_df$QE[which(ES_rma_df$corr == 0)],
                      QE_cor = ES_rma_df$QE[which(ES_rma_df$corr == 1)]) %>% 
  mutate(tau_diff = tau_raw - tau_cor,
         QE_raw_str = paste0(round(QE_raw, 3), " (", k-1, ")"),
         QE_cor_str = paste0(round(QE_cor, 3), " (", k-1, ")"),
         QEp_raw = ifelse(QEp_raw < .001, yes = "<.001", no = round(QEp_raw, 3)),
         QEp_cor = ifelse(QEp_cor < .001, yes = "<.001", no = round(QEp_cor, 3))
         ) %>% 
  mutate_if(is.numeric, round, digits = 3)

write.csv(res_tab[nn_eff_idx,] %>% 
            select(MASC, tau_raw, QE_raw_str, QEp_raw, tau_cor, QE_cor_str, QEp_cor), 
          here("Tables/Heterogeneity_ES.csv"),
          row.names = F)


mean_es_res_tab <- data.frame(
  mu_raw = ES_rma_df$mu[which(ES_rma_df$corr == 0)],
  mu_cor = ES_rma_df$mu[which(ES_rma_df$corr == 1)],
  mu_ll_raw = ES_rma_df$ci.lb[which(ES_rma_df$corr == 0)],
  mu_ul_raw = ES_rma_df$ci.ub[which(ES_rma_df$corr == 0)],
  mu_ll_cor = ES_rma_df$ci.lb[which(ES_rma_df$corr == 1)],
  mu_ul_cor = ES_rma_df$ci.ub[which(ES_rma_df$corr == 1)],
  mu_SE_raw = ES_rma_df$se[which(ES_rma_df$corr == 0)],
  mu_SE_cor = ES_rma_df$se[which(ES_rma_df$corr == 1)],
  pval_raw = ES_rma_df$pval[which(ES_rma_df$corr == 0)],
  pval_cor = ES_rma_df$pval[which(ES_rma_df$corr == 1)],
  MASC = ES_rma_df$MASC[which(ES_rma_df$corr == 0)]
) %>% 
  mutate(mu_str_raw = paste0(sub("^(-?)0.", "\\1.", sprintf("%.3f", .$mu_raw)), " (", sub("^(-?)0.", "\\1.", sprintf("%.3f", .$mu_SE_raw)), ")"),
         mu_str_cor = paste0(sub("^(-?)0.", "\\1.", sprintf("%.3f", .$mu_cor)), " (", sub("^(-?)0.", "\\1.", sprintf("%.3f", .$mu_SE_cor)), ")"),
         pval_raw = ifelse(pval_raw < .001, yes = "<.001", no = sub("^(-?)0.", "\\1.", sprintf("%.3f", .$pval_raw))),
         pval_cor = ifelse(pval_cor < .001, yes = "<.001", no = sub("^(-?)0.", "\\1.", sprintf("%.3f", .$pval_cor)))) %>% 
  select(MASC, mu_str_raw, pval_raw, mu_str_cor, pval_cor)

write.csv(mean_es_res_tab,
          file = here("Tables/Mean_ES.csv"),
          row.names = FALSE)


alpha_res_tab <- B_alpha_rma_df[nn_eff_idx,] %>% 
  mutate_if(is.numeric, round, digits = 3) %>% 
  mutate(mu_alpha_str = paste0(round(mu_alpha, 3), " [", round(mu_alpha_ll, 3), ":", round(mu_alpha_ul, 3), "]"),
         QE_str = paste0(round(QE, 3), " (", k-1, ")"),
         QEp = ifelse(B_alpha_rma_df[nn_eff_idx,]$QEp < .001, yes = "<.001", no = sub("^(-?)0.", "\\1.", sprintf("%.3f", B_alpha_rma_df[nn_eff_idx,]$QEp)))) %>% 
  select(MASC, mu_alpha_str, tau_alpha, QE_str, QEp) 

write.csv(alpha_res_tab,
          here("Tables/Results_RMA_alpha.csv"),
          row.names = FALSE)

forest_plot_rel(ES_rma_df[which(ES_rma_df$MASC == "Nosek_Explicit_Art"),], agg_L[[which(MASC_names == "Nosek_Explicit_Art")]])

ggsave(filename = here("Graphics/forest_Nosek.png"),
       plot = last_plot(),
       width = 6, height = 4)

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

