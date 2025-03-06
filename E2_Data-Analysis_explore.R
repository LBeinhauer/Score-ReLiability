## Score ReLiability

# E1.5

# Clean & extract data


# library loading and installing as necessary


# relevant libraries required for this script
packages <- c("magrittr", "dplyr", "boot", "here")

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

MASC_names <- substr(list.files(here("Data/Extracted (Project) Data")), 1, nchar(list.files(here("Data/Extracted (Project) Data")))-4) 




data_files <- list.files(here("Data/Extracted (Project) Data"), full.names = TRUE)


data.list <- lapply(data_files, read.csv)


data.list[[4]]$source <- substr(data.list[[4]]$source, start = nchar(data.list[[4]]$source) - 2, stop = nchar(data.list[[4]]$source))


effect_index <- MASC_names %in% c("Albarracin_Priming_SAT", 
                                  "Carter_Flag_Priming", "Caruso_Currency_Priming",
                                  "Dijksterhuis_trivia", "Finkel_Exit_Forgiveness", 
                                  "Finkel_Neglect_Forgiveness",
                                  "Giessner_Vertical_Position", "Hart_Criminal_Intentionality",    
                                  "Hart_Detailed_Processing", "Hart_Intention_Attribution",       
                                  "Husnu_Imagined_Contact", "Nosek_Explicit_Art",
                                  "Nosek_Explicit_Math", "PSACR001_anxiety_int", 
                                  "PSACR001_behav_int", 
                                  "PSACR002_neg_photo", "Shnabel_Willingness_Reconcile_Rev",
                                  "Shnabel_Willingness_Reconcile_RPP", "Srull_Behaviour_Hostility",
                                  "Srull_Ronald_Hostility", "Tversky_Directionality_Similarity1", 
                                  "Zhong_Desirability_Cleaning"
)


nn_eff_idx <- which(ES_rma_df$pval[ES_rma_df$corr == 0] <= .05)



# generate estimates of ln-error score variance and its associated standard error using bootstrapping
varE_est.L <- lapply(which(effect_index)[nn_eff_idx], FUN = function(x){
  
  # sometimes calculation of SE might lead to issues, as negative variances can not be log-transformed
  #  therefore, function is run within a tryCatch-environment, so the script does not crash
  # apply_Bootstrap_SE_Project.specific is the project-specific function
  tryCatch(apply_Bootstrap_SE_Project.specific(na.omit(data.list[[x]]), component = "E", R = 100),
           
           # print error to console
           error = function(e)(cat("ERROR: ", conditionMessage(e), " - ",
                                   substr(names(data.list), 
                                          (regexpr("Project) Data/", names(data.list)) + 14), 
                                          (nchar(names(data.list))-4))[x], 
                                   " - ", x, "\n")))
})


# Perform random-effects meta-analysis on estimates of ln-error score variance using metafor
varE_rma.list <- lapply(seq_along(varE_est.L), FUN = function(x){
  
  # at times estimates of ln-varE & SE may be NA, as negative estimates can't be transformed
  #  therefore, function is nested in tryCatch, so it doesn't break down with errors
  tryCatch(metafor::rma(measure = "GEN", method = "REML", 
                        yi = varE_est.L[[x]]$var.est, 
                        sei = varE_est.L[[x]]$SE),
           
           # print error to console
           error = function(e)(cat("ERROR: ", conditionMessage(e), " - ",
                                   substr(names(data.list), 
                                          (regexpr("Project) Data/", names(data.list)) + 14), 
                                          (nchar(names(data.list))-4))[x], 
                                   " - ", x, "\n")))
  
})

# add names to list-object
names(varE_rma.list) <- names(data.list)



# generate estimates of ln-observed score variance for each sample & projects
varX_est.L <- lapply(which(effect_index)[nn_eff_idx], FUN = function(x){
  
  # sometimes calculation of SE might lead to issues, as negative variances can not be log-transformed
  #  therefore, function is run within a tryCatch-environment, so the script does not crash
  # apply_Bootstrap_SE_Project.specific is the project-specific function
  tryCatch(apply_Bootstrap_SE_Project.specific(na.omit(data.list[[x]]), component = "X", R = 100),
           
           # print error to console
           error = function(e)(cat("ERROR: ", conditionMessage(e), " - ",
                                   substr(names(data.list), 
                                          (regexpr("Project) Data/", names(data.list)) + 14), 
                                          (nchar(names(data.list))-4))[x], 
                                   " - ", x, "\n")))
})




# Perform random-effects meta-analysis on estimates of ln-observed score variance using metafor
varX_rma.list <- lapply(seq_along(varX_est.L), FUN = function(x){
  
  # at times estimates of ln-varX & SE may be NA, as negative estimates can't be transformed
  #  therefore, function is nested in tryCatch, so it doesn't break down with errors
  tryCatch(metafor::rma(measure = "GEN", method = "REML", 
                        yi = varX_est.L[[x]]$var.est, 
                        sei = varX_est.L[[x]]$SE),
           
           # print error to console
           error = function(e)(cat("ERROR: ", conditionMessage(e), " - ",
                                   substr(names(data.list), 
                                          (regexpr("Project) Data/", names(data.list)) + 14), 
                                          (nchar(names(data.list))-4))[x], 
                                   " - ", x, "\n")))
  
})



vars_rma_L <- lapply(1:length(varE_rma.list), FUN = function(x){
  
  data.frame(mu_E = bt_var_m(varE_rma.list[[x]]),
             mu_X = bt_var_m(varX_rma.list[[x]]),
             tau2_E = bt_var_v(varE_rma.list[[x]]),
             tau2_X = bt_var_v(varX_rma.list[[x]]))
  
})

vars_rma_df <- do.call(rbind, vars_rma_L) %>% 
  mutate(mu_T = mu_X - mu_E,
         tau2_T = tau2_X - tau2_E) %>% 
  mutate(CV_E = sqrt(tau2_E)/mu_E,
         CV_X = sqrt(tau2_X)/mu_X,
         CV_T = sqrt(tau2_T)/mu_T) %>% 
  mutate(R1 = mu_T / mu_X,
         R2 = tau2_T / tau2_X) %>% 
  mutate(MASC = ES_rma_df$MASC[ES_rma_df$corr == 0][nn_eff_idx]) %>% 
  mutate_if(is.numeric, round, 3)





write.csv(vars_rma_df %>% 
            select(MASC, mu_X, mu_T, CV_X, CV_T, R1, R2), here("Tables/Variances_analysis.csv"), row.names = FALSE)


agg_L <- lapply(agg_L, FUN = function(x){
  y <- x %>% mutate(alpha_pooled = (n1 / (n0 + n1)) * alpha1 + (n0 / (n0 + n1)) * alpha0)
  
  return(y)
}) 
  

