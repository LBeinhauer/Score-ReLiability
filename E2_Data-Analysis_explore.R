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


# generate estimates of ln-error score variance and its associated standard error using bootstrapping
varE_est.L <- lapply(seq_along(data.list), FUN = function(x){
  
  # sometimes calculation of SE might lead to issues, as negative variances can not be log-transformed
  #  therefore, function is run within a tryCatch-environment, so the script does not crash
  # apply_Bootstrap_SE_Project.specific is the project-specific function
  tryCatch(apply_Bootstrap_SE_Project.specific(na.omit(data.list[[x]]), component = "E"),
           
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
varX_est.L <- lapply(seq_along(data.list), FUN = function(x){
  
  # sometimes calculation of SE might lead to issues, as negative variances can not be log-transformed
  #  therefore, function is run within a tryCatch-environment, so the script does not crash
  # apply_Bootstrap_SE_Project.specific is the project-specific function
  tryCatch(apply_Bootstrap_SE_Project.specific(na.omit(data.list[[x]]), component = "X"),
           
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



vars_df <- data.frame(mu_varE = sapply(varE_rma.list, FUN = function(x){bt_var_m(x)}),
                      tau2_varE = sapply(varE_rma.list, FUN = function(x){bt_var_v(x)}),
                      mu_varX = sapply(varX_rma.list, FUN = function(x){bt_var_m(x)}),
                      tau2_varX = sapply(varX_rma.list, FUN = function(x){bt_var_v(x)}))

vars_df$tau2_varE / vars_df$tau2_varX
vars_df$mu_varE / vars_df$mu_varX

write.csv(vars_df, here("Data/Processed/Variances_analysis.csv"), row.names = FALSE)

