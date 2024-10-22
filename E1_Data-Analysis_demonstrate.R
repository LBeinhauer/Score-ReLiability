## Score ReLiability

# E1

# Clean & extract data


# library loading and installing as necessary


# relevant libraries required for this script
packages <- c("tidyverse", "here", "metafor")

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



data_files <- list.files(here("Data/Extracted (Project) Data"), full.names = TRUE)


data.list <- lapply(data_files, read.csv)


data.list[[4]]$source <- substr(data.list[[4]]$source, start = nchar(data.list[[4]]$source) - 2, stop = nchar(data.list[[4]]$source))

MASC_names <- substr(list.files(here("Data/Extracted (Project) Data")), 1, nchar(list.files(here("Data/Extracted (Project) Data")))-4) 


# generate aggregates for each lab/country/source within a single meta-analysis. This is done by
#  running an lapply-loop across all MASCs (meta-analyses, see MetaPipeX, Fünderich et al. 2024)
#  and, nested within each loop, a second lapply-loop across all labs/countries/sources (samples)
agg_L <- lapply(data.list, FUN = function(x){
  
  # loop across labs/countries/sources (samples) within a single MASC
  MASC_MD_L <- lapply(unique(x$source), FUN = function(lab){
    
    # data for respective country
    d <- x[x$source == lab,]
    
    # data for respective treatment (1) or control (0) subgroup
    d1 <- na.omit(d[d$group == 1, -c(grep("group", names(d)), grep("source", names(d)))])
    d0 <- na.omit(d[d$group == 0, -c(grep("group", names(d)), grep("source", names(d)))])
    
    # sample size
    n1 <- nrow(d1)
    n0 <- nrow(d0)
    
    # covariance matrix
    C1 <- cov(d1)
    C0 <- cov(d0)
    
    # number of items in scale (identical across groups)
    j <- dim(C1)[1]
    
    # estimate of Cronbach's Alpha
    alpha1 <- (1 - sum(diag(C1))/sum(C1)) * (j/(j - 1))
    alpha0 <- (1 - sum(diag(C0))/sum(C0)) * (j/(j - 1))
    
    # Standard Error for transformed [log(1 - alpha)] Cronbach's Alpha
    SE_Balpha1 <- sqrt((2 * j)/((j - 1) * (n1 - 2))) 
    SE_Balpha0 <- sqrt((2 * j)/((j - 1) * (n1 - 2))) 
    
    # vector of individua scores (means across items)
    mv1 <- rowMeans(d1)
    mv0 <- rowMeans(d0)
    
    # group means
    m1 <- mean(mv1)
    m0 <- mean(mv0)
    
    # group sample standard deviation
    sd1 <- sqrt(mean((mv1 - m1)^2))
    sd0 <- sqrt(mean((mv0 - m0)^2))
    
    # group population standard devation
    sigma1 <- sd(mv1)
    sigma0 <- sd(mv0)
    
    # mean difference
    MD <- m1-m0
    
    # pooled standard deviation
    pooled_sd <- sqrt(((n1-1)*sd1^2 + (n0-1)*sd0^2)/(n1+n0-2))
    
    # pooled standard deviation, corrected for imperfect reliability
    pooled_sd_corr <- sqrt(((n1-1)*((sd1^2)*alpha1) + (n0-1)*((sd0^2)*alpha0))/(n1+n0-2))
    
    # uncorrected estiamte of Cohen's d
    d_raw <- MD/pooled_sd
    
    # corrected estimate of Cohen's d
    d_corr <- MD/pooled_sd_corr
    
    # Standard Error of uncorrected Cohen's d
    SE_d <- sqrt(((n0 + n1)/(n0 * n1)) + ((d_raw^2) / (2 * (n0 + n1))))
    
    # Standard Error of corrected Cohen's d
    SE_d_corr <- sqrt(SE_d^2 * (d_corr/d_raw)^2)
    
    # return a data-frame of all aggregate estimates made
    return(data.frame(n1 = n1, 
                      n0 = n0,
                      j = j,
                      m1 = m1,
                      m0 = m0,
                      sd1 = sd1,
                      sd0 = sd0,
                      sigma1 = sigma1,
                      sigma0 = sigma0,
                      alpha1 = alpha1,
                      alpha0 = alpha0,
                      SE_Balpha1 = SE_Balpha1,
                      SE_Balpha0 = SE_Balpha0,
                      MD = MD,
                      pooled_sd = pooled_sd,
                      pooled_sd_corr = pooled_sd_corr,
                      d_raw = d_raw,
                      d_corr = d_corr,
                      SE_d = SE_d,
                      SE_d_corr = SE_d_corr,
                      source = lab))
    
  })
  
  # prepare a single data-frame for each MASC
  MASC_MD_df <- do.call(rbind, MASC_MD_L)
  
  return(MASC_MD_df)
  

})


saveRDS(agg_L, file = here("Data/Processed/Aggregates_simple.csv"))



# perform meta-analyses of transformed score reliability
B_alpha_rma <- lapply(agg_L, FUN = function(x){
  
  # if two group-structure was used, we need to pool the estimates of cronbach's alpha
  if(!is.na(x$n0[1])){
    
    # pooled cronbach's alpha
    alpha_pooled <- (x$n1 / (x$n0 + x$n1)) * x$alpha1 + (x$n0 / (x$n0 + x$n1)) * x$alpha0
    
    # standard error for pooled cronbach's alpha
    SE_B_alpha_pooled <- sqrt((2 * x$j)/((x$j - 1) * (x$n1 + x$n0 - 2))) 
    
    # use metafor to perform meta-analysis
    BA_rma <- metafor::rma(yi = log(1 - alpha_pooled),
                           sei = SE_B_alpha_pooled,
                           measure = "GEN",
                           method = "REML")
    
  }else{
    
    # use metafor to perform meta-analysis directly, as no pooling is required
    BA_rma <- metafor::rma(yi = log(1 - x$alpha1),
                           sei = x$SE_Balpha1,
                           measure = "GEN",
                           method = "REML")
    
    
  }
  
  
  # return data-frame with all important aggregates
  return(data.frame(tau2_BA = BA_rma$tau2,
                    mu_BA = BA_rma$b[1],
                    tau_alpha = sqrt(var_Bonnett_backtransformed(BA_rma)),
                    mu_alpha = mean_Bonnett_backtransformed(BA_rma),
                    QE = BA_rma$QE,
                    k = nrow(x),
                    QEp = BA_rma$QEp,
                    I2 = BA_rma$I2,
                    H2 = BA_rma$H2))
  
})

# construct a single data-frame from list element
B_alpha_rma_df <- do.call(rbind, B_alpha_rma) %>% 
  mutate(MASC = MASC_names)

write.csv(B_alpha_rma_df, here("Data/Processed/Reliability_analysis.csv"), row.names = FALSE)



# run meta-analyses on Cohen's d (corrected | uncorrected) on all relevant MASCs
d_rma_full_L <- lapply(agg_L, FUN = function(x){
  
  # use metafor to run meta-analysis of uncorrected Cohen's d
  d_rma_raw <- metafor::rma(yi = x$d_raw,
                            sei = x$SE_d,
                            measure = "GEN",
                            method = "REML")
  
  # use metafor to run meta-analysis of corrected Cohen's d
  d_rma_corr <- metafor::rma(yi = x$d_cor,
                             sei = x$SE_d_corr,
                             measure = "GEN",
                             method = "REML")
  
  return(list(rma_raw = d_rma_raw,
              rma_corr = d_rma_corr))
})


# add names to list-elements
names(d_rma_full_L) <- MASC_names[effect_index]


ES_rma_df <- data.frame(mu = unlist(lapply(d_rma_full_L, FUN = function(x){c(x$rma_raw$b[1], x$rma_cor$b[1])})),
                        se = unlist(lapply(d_rma_full_L, FUN = function(x){c(x$rma_raw$se, x$rma_cor$se)})),
                        ci.lb = unlist(lapply(d_rma_full_L, FUN = function(x){c(x$rma_raw$ci.lb, x$rma_cor$ci.lb)})),
                        ci.ub = unlist(lapply(d_rma_full_L, FUN = function(x){c(x$rma_raw$ci.ub, x$rma_cor$ci.ub)})),
                        tau = sqrt(unlist(lapply(d_rma_full_L, FUN = function(x){c(x$rma_raw$tau2, x$rma_cor$tau2)}))),
                        QE = unlist(lapply(d_rma_full_L, FUN = function(x){c(x$rma_raw$QE, x$rma_cor$QE)})),
                        k = unlist(lapply(d_rma_full_L, FUN = function(x){c(x$rma_raw$k, x$rma_cor$k)})),
                        QEp = unlist(lapply(d_rma_full_L, FUN = function(x){c(x$rma_raw$QEp, x$rma_cor$QEp)})),
                        I2 = unlist(lapply(d_rma_full_L, FUN = function(x){c(x$rma_raw$I2, x$rma_cor$I2)})),
                        corr = rep(c(0,1), 4),
                        MASC = rep(MASC_names, each = 2)) %>% 
  mutate(CV = tau/mu)


write.csv(ES_rma_df, here("Data/Processed/Aggregates_ES_analysis.csv"), row.names = FALSE)


