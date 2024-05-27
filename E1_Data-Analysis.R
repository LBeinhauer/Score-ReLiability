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



# load initial data-file
# "wide, exclusions" - should be the correct, cleaned data-file
psacr002_dat <- read.csv(here("Data/Downloaded/data_wide_exclusion_12.23.2020.csv"))

# identify variables which contain responses on stimuli (positive and negative responses separate)
pos_items <- names(psacr002_dat)[grep("positiveemotion_photo", names(psacr002_dat))]
neg_items <- names(psacr002_dat)[grep("negativeemotion_photo", names(psacr002_dat))]

pos_state_items <- names(psacr002_dat)[grep("C2_positiveemotion_state.*[^0-9]$", names(psacr002_dat))]
neg_state_items <- names(psacr002_dat)[grep("C2_negativeemotion_state.*[^0-9]$", names(psacr002_dat))]

pos_antc_items <- names(psacr002_dat)[grep("C2_positiveemotion_anticipated.*[^0-9]$", names(psacr002_dat))]
neg_antc_items <- names(psacr002_dat)[grep("C2_negativeemotion_anticipated.*[^0-9]$", names(psacr002_dat))]

# identify which countries contain sufficient (more than 20) respondents
samples <- as.data.frame(table(psacr002_dat$country))$Var1[as.data.frame(table(psacr002_dat$country))$Freq >= 200]


#### PHOTOGRAPHS

# formulate data-file containing responses on positive photo items, information on the grouping
#  variable, and the respondent's country
psacr002_pos <- psacr002_dat %>% 
  # The original authors used more than 2 groups. Here we collapse across several groups
  mutate(group = ifelse(condition %in% c("Active Control", "Passive Control"), yes = 0, no = 1)) %>% 
  select(c(all_of(pos_items), "group", "country")) %>% 
  filter(country %in% samples)


# formulate data-file containing responses on negative photo items, information on the grouping
#  variable, and the respondent's country
psacr002_neg <- psacr002_dat %>% 
  # The original authors used more than 2 groups. Here we collapse across several groups
  mutate(group = ifelse(condition %in% c("Active Control", "Passive Control"), yes = 0, no = 1)) %>% 
  select(c(all_of(neg_items), "group", "country")) %>% 
  filter(country %in% samples)



#### STATE EMOTIONS

# formulate data-file containing responses on positive state items, information on the grouping
#  variable, and the respondent's country
psacr002_pos_state <- psacr002_dat %>% 
  # The original authors used more than 2 groups. Here we collapse across several groups
  mutate(group = ifelse(condition %in% c("Active Control", "Passive Control"), yes = 0, no = 1)) %>% 
  select(c(all_of(pos_state_items), "group", "country")) %>% 
  filter(country %in% samples)


# formulate data-file containing responses on negative state items, information on the grouping
#  variable, and the respondent's country
psacr002_neg_state <- psacr002_dat %>% 
  # The original authors used more than 2 groups. Here we collapse across several groups
  mutate(group = ifelse(condition %in% c("Active Control", "Passive Control"), yes = 0, no = 1)) %>% 
  select(c(all_of(neg_state_items), "group", "country")) %>% 
  filter(country %in% samples)



#### ANTICIPATED EMOTIONS

# formulate data-file containing responses on positive anticipation items, information on the grouping
#  variable, and the respondent's country
psacr002_pos_antc <- psacr002_dat %>% 
  # The original authors used more than 2 groups. Here we collapse across several groups
  mutate(group = ifelse(condition %in% c("Active Control", "Passive Control"), yes = 0, no = 1)) %>% 
  select(c(all_of(pos_antc_items), "group", "country")) %>% 
  filter(country %in% samples)


# formulate data-file containing responses on negative anticipation items, information on the grouping
#  variable, and the respondent's country
psacr002_neg_antc <- psacr002_dat %>% 
  # The original authors used more than 2 groups. Here we collapse across several groups
  mutate(group = ifelse(condition %in% c("Active Control", "Passive Control"), yes = 0, no = 1)) %>% 
  select(c(all_of(neg_antc_items), "group", "country")) %>% 
  filter(country %in% samples)





# function to compute (standardized) mean differences, in form of Cohen's d
# the function takes the data-files created above these lines.
# both uncorrected standardized mean differences, as corrected standardized mean differences
#  using score reliability of each group are computed.
compute_SMDs <- function(data){
  
  # using apply to loop across countries k
  SMDs.L <- lapply(unique(data$country), FUN = function(k){
    
    # differentiate between control and reconstrual/repurposing groups
    dat0 <- data[data$"country" == k & data$"group" == 0, 
                         !(names(data) %in% c("country", "group"))] # control groups
    dat1 <- data[data$"country" == k & data$"group" == 1, 
                         !(names(data) %in% c("country", "group"))] # reconstrual/repurposing groups
    
    # generate dependent variable, in form of mean rating across items per individual
    #  separate for groups
    dv0 <- na.omit(rowMeans(dat0))
    dv1 <- na.omit(rowMeans(dat1))
    
    # generate group means
    m0 <- mean(dv0, na.rm = T)
    m1 <- mean(dv1, na.rm = T)
    
    # generate group variances (sample variances)
    var0 <- mean((dv0 - mean(dv0))^2)
    var1 <- mean((dv1 - mean(dv1))^2)
    
    # generate group standard deivations (sample standard deviations)
    sd0 <- sqrt(var0)
    sd1 <- sqrt(var1)
    
    # get group sizes
    n0 <- length(dv0)
    n1 <- length(dv1)
    
    # compute Cronbach's Alpha using the covariance matrix 
    # covariance matrix
    C0 <- cov(dat0, use = "complete.obs")
    # number of items
    j0 <- dim(C0)[1]
    # generate Cronbach's Alpha
    alpha0 <- (1 - sum(diag(C0))/sum(C0)) * (j0/(j0 - 1))
    
    # compute Cronbach's Alpha using the covariance matrix 
    # covariance matrix
    C1 <- cov(dat1, use = "complete.obs")
    # number of items
    j1 <- dim(C1)[1]
    # generate Cronbach's Alpha
    alpha1 <- (1 - sum(diag(C1))/sum(C1)) * (j1/(j1 - 1))
    
    # compute uncorrected Cohen's d
    d_raw <- (m1 - m0)/sqrt(((n0-1)*var0 + (n1-1)*var1)/(n0+n1-2))
    
    # correct the variances using estimates of Cronbach's Alpha
    var0_c <- var0 * alpha0
    var1_c <- var1 * alpha1
    
    # compute corrected Cohen's d
    d_cor <- (m1 - m0)/sqrt(((n0-1)*var0_c + (n1-1)*var1_c)/(n0+n1-2))
    
    # generate sampling variance/standard error for uncorrected Cohen's d
    V_d_raw <- ((n0 + n1)/(n0 * n1)) + ((d_raw^2) / (2 * (n0 + n1)))
    SE_d_raw <- sqrt(V_d_raw)
    
    # generate sampling variance/standard error for corrected Cohen's d
    V_d_cor <- V_d_raw * (d_cor/d_raw)^2
    SE_d_cor <- sqrt(V_d_cor)
    
    # return estimates in standardized data-frame format
    return(data.frame(d_raw = d_raw,
                      SE_d_raw = SE_d_raw,
                      d_cor = d_cor,
                      SE_d_cor,
                      MD = m0 - m1,
                      SE_MD = sqrt(var0/n0 + var1/n1),
                      alpha0 = alpha0,
                      alpha1 = alpha1,
                      j = j0,
                      n0 = n0,
                      n1 = n1,
                      country = substring(k, nchar(k)-2, nchar(k))))
  })
  
  # collapse across list-elements, so a single data.frame contains all parameters
  return(do.call(rbind, SMDs.L))
  
}


var_Bonnett_backtransformed <- function(mean_x, var_x){
  (((-exp(mean_x))^2) * var_x) + (.5*((-exp(mean_x))^2)*(var_x^2)) + ((-exp(mean_x)) * (-exp(mean_x)) * (var_x^2))
}

mean_Bonnett_backtransformed <- function(mean_x, var_x){
  1 - exp(mean_x) + ((-exp(mean_x)) / 2) * var_x
}




##### PHOTOGRAPHS

# fit function to data-set containing positiveness-ratings of photographs
pos_SMDs.df <- compute_SMDs(psacr002_pos)


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

# print Coefficient of Variation as heterogeneity of corrected and uncorrected Cohen's d, for positiveness-ratings
cat("pos raw CV: ", round(sqrt(pos_SMDs_rma_raw$tau2)/pos_SMDs_rma_raw$b[1], 3), "\n", "pos cor CV: ", round(sqrt(pos_SMDs_rma_cor$tau2)/pos_SMDs_rma_cor$b[1], 3), sep = "")
# print raw heterogeneity in tau of corrected and uncorrected Cohen's d, for positiveness-ratings
cat("pos raw tau: ", round(sqrt(pos_SMDs_rma_raw$tau2), 3), "\n", "pos cor tau: ", round(sqrt(pos_SMDs_rma_cor$tau2), 3), sep = "")


# fit a meta-analytic model to estimates of Cronbach's Alpha, for the control groups

pos_rel_0 <- metafor::rma(yi = log(1 - alpha0),
                          sei = (2*j) / ((j - 1)*(n0 - 1)),
                          measure = "GEN",
                          method = "REML",
                          data = pos_SMDs.df)

sqrt(var_Bonnett_backtransformed(pos_rel_0$b[1], pos_rel_0$tau2))
mean_Bonnett_backtransformed(pos_rel_0$b[1], pos_rel_0$tau2)


# for the treatment groups

pos_rel_1 <- metafor::rma(yi = log(1 - alpha1),
                     sei = (2*j) / ((j - 1)*(n1 - 1)),
                     measure = "GEN",
                     method = "REML",
                     data = pos_SMDs.df)

sqrt(var_Bonnett_backtransformed(pos_rel_1$b[1], pos_rel_1$tau2))
mean_Bonnett_backtransformed(pos_rel_1$b[1], pos_rel_1$tau2)




# fit function to data-set containing negativeness-ratings of photographs
neg_SMDs.df <- compute_SMDs(psacr002_neg)


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

# print Coefficient of Variation as heterogeneity of corrected and uncorrected Cohen's d, for negativeness-ratings
cat("neg raw CV: ", round(sqrt(neg_SMDs_rma_raw$tau2)/neg_SMDs_rma_raw$b[1], 3), "\n", "neg cor CV: ", round(sqrt(neg_SMDs_rma_cor$tau2)/neg_SMDs_rma_cor$b[1], 3), sep = "")
# print raw heterogeneity in tau of corrected and uncorrected Cohen's d, for negativeness-ratings
cat("neg raw tau: ", round(sqrt(neg_SMDs_rma_raw$tau2), 3),"\n", "neg cor tau: ", round(sqrt(neg_SMDs_rma_cor$tau2), 3), sep = "")



# fit a meta-analytic model to estimates of Cronbach's Alpha, for the control groups

neg_rel_0 <- metafor::rma(yi = log(1 - alpha0),
                          sei = (2*j) / ((j - 1)*(n0 - 1)),
                          measure = "GEN",
                          method = "REML",
                          data = neg_SMDs.df)

sqrt(var_Bonnett_backtransformed(neg_rel_0$b[1], neg_rel_0$tau2))
mean_Bonnett_backtransformed(neg_rel_0$b[1], neg_rel_0$tau2)


# for the treatment groups

neg_rel_1 <- metafor::rma(yi = log(1 - alpha1),
                          sei = (2*j) / ((j - 1)*(n1 - 1)),
                          measure = "GEN",
                          method = "REML",
                          data = neg_SMDs.df)

sqrt(var_Bonnett_backtransformed(neg_rel_1$b[1], neg_rel_1$tau2))
mean_Bonnett_backtransformed(neg_rel_1$b[1], neg_rel_1$tau2)




##### STATE EMOTIONS

# fit function to data-set containing positiveness-ratings of photographs
pos_state_SMDs.df <- compute_SMDs(psacr002_pos_state)


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

# print Coefficient of Variation as heterogeneity of corrected and uncorrected Cohen's d, for positiveness-ratings
cat("pos state raw CV: ", round(sqrt(pos_state_SMDs_rma_raw$tau2)/pos_state_SMDs_rma_raw$b[1], 3), "\n", "pos state cor CV: ", round(sqrt(pos_state_SMDs_rma_cor$tau2)/pos_state_SMDs_rma_cor$b[1], 3), sep = "")
# print raw heterogeneity in tau of corrected and uncorrected Cohen's d, for positiveness-ratings
cat("pos state raw tau: ", round(sqrt(pos_state_SMDs_rma_raw$tau2), 3), "\n", "pos state cor tau: ", round(sqrt(pos_state_SMDs_rma_cor$tau2), 3), sep = "")


# fit a meta-analytic model to estimates of Cronbach's Alpha, for the control groups

pos_state_rel_0 <- metafor::rma(yi = log(1 - alpha0),
                                sei = (2*j) / ((j - 1)*(n0 - 1)),
                                measure = "GEN",
                                method = "REML",
                                data = pos_state_SMDs.df)

sqrt(var_Bonnett_backtransformed(pos_state_rel_0$b[1], pos_state_rel_0$tau2))
mean_Bonnett_backtransformed(pos_state_rel_0$b[1], pos_state_rel_0$tau2)


# for the treatment groups

pos_state_rel_1 <- metafor::rma(yi = log(1 - alpha1),
                                sei = (2*j) / ((j - 1)*(n1 - 1)),
                                measure = "GEN",
                                method = "REML",
                                data = pos_state_SMDs.df)

sqrt(var_Bonnett_backtransformed(pos_state_rel_1$b[1], pos_state_rel_1$tau2))
mean_Bonnett_backtransformed(pos_state_rel_1$b[1], pos_state_rel_1$tau2)






# fit function to data-set containing negativeness-ratings of photographs
neg_state_SMDs.df <- compute_SMDs(psacr002_neg_state)


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

# print Coefficient of Variation as heterogeneity of corrected and uncorrected Cohen's d, for negativeness-ratings
cat("neg state raw CV: ", round(sqrt(neg_state_SMDs_rma_raw$tau2)/neg_state_SMDs_rma_raw$b[1], 3), "\n", "neg state cor CV: ", round(sqrt(neg_state_SMDs_rma_cor$tau2)/neg_state_SMDs_rma_cor$b[1], 3), sep = "")
# print raw heterogeneity in tau of corrected and uncorrected Cohen's d, for negativeness-ratings
cat("neg state raw tau: ", round(sqrt(neg_state_SMDs_rma_raw$tau2), 3),"\n", "neg state cor tau: ", round(sqrt(neg_state_SMDs_rma_cor$tau2), 3), sep = "")


# fit a meta-analytic model to estimates of Cronbach's Alpha, for the control groups

neg_state_rel_0 <- metafor::rma(yi = log(1 - alpha0),
                                sei = (2*j) / ((j - 1)*(n0 - 1)),
                                measure = "GEN",
                                method = "REML",
                                data = neg_state_SMDs.df)

sqrt(var_Bonnett_backtransformed(neg_state_rel_0$b[1], neg_state_rel_0$tau2))
mean_Bonnett_backtransformed(neg_state_rel_0$b[1], neg_state_rel_0$tau2)


# for the treatment groups

neg_state_rel_1 <- metafor::rma(yi = log(1 - alpha1),
                                sei = (2*j) / ((j - 1)*(n1 - 1)),
                                measure = "GEN",
                                method = "REML",
                                data = neg_state_SMDs.df)

sqrt(var_Bonnett_backtransformed(neg_state_rel_1$b[1], neg_state_rel_1$tau2))
mean_Bonnett_backtransformed(neg_state_rel_1$b[1], neg_state_rel_1$tau2)








##### ANTICIPATED EMOTIONS


# fit function to data-set containing positiveness-ratings of photographs
pos_antc_SMDs.df <- compute_SMDs(psacr002_pos_antc)


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

# print Coefficient of Variation as heterogeneity of corrected and uncorrected Cohen's d, for positiveness-ratings
cat("pos antc raw CV: ", round(sqrt(pos_antc_SMDs_rma_raw$tau2)/pos_antc_SMDs_rma_raw$b[1], 3), "\n", "pos antc cor CV: ", round(sqrt(pos_antc_SMDs_rma_cor$tau2)/pos_antc_SMDs_rma_cor$b[1], 3), sep = "")
# print raw heterogeneity in tau of corrected and uncorrected Cohen's d, for positiveness-ratings
cat("pos antc raw tau: ", round(sqrt(pos_antc_SMDs_rma_raw$tau2), 3), "\n", "pos antc cor tau: ", round(sqrt(pos_antc_SMDs_rma_cor$tau2), 3), sep = "")



# fit a meta-analytic model to estimates of Cronbach's Alpha, for the control groups

pos_antc_rel_0 <- metafor::rma(yi = log(1 - alpha0),
                                sei = (2*j) / ((j - 1)*(n0 - 1)),
                                measure = "GEN",
                                method = "REML",
                                data = pos_antc_SMDs.df)

sqrt(var_Bonnett_backtransformed(pos_antc_rel_0$b[1], pos_antc_rel_0$tau2))
mean_Bonnett_backtransformed(pos_antc_rel_0$b[1], pos_antc_rel_0$tau2)


# for the treatment groups

pos_antc_rel_1 <- metafor::rma(yi = log(1 - alpha1),
                                sei = (2*j) / ((j - 1)*(n1 - 1)),
                                measure = "GEN",
                                method = "REML",
                                data = pos_antc_SMDs.df)

sqrt(var_Bonnett_backtransformed(pos_antc_rel_1$b[1], pos_antc_rel_1$tau2))
mean_Bonnett_backtransformed(pos_antc_rel_1$b[1], pos_antc_rel_1$tau2)





# fit function to data-set containing negativeness-ratings of photographs
neg_antc_SMDs.df <- compute_SMDs(psacr002_neg_antc)


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

# print Coefficient of Variation as heterogeneity of corrected and uncorrected Cohen's d, for negativeness-ratings
cat("neg antc raw CV: ", round(sqrt(neg_antc_SMDs_rma_raw$tau2)/neg_antc_SMDs_rma_raw$b[1], 3), "\n", "neg cor CV: ", round(sqrt(neg_antc_SMDs_rma_cor$tau2)/neg_antc_SMDs_rma_cor$b[1], 3), sep = "")
# print raw heterogeneity in tau of corrected and uncorrected Cohen's d, for negativeness-ratings
cat("neg antc raw tau: ", round(sqrt(neg_antc_SMDs_rma_raw$tau2), 3),"\n", "neg cor tau: ", round(sqrt(neg_antc_SMDs_rma_cor$tau2), 3), sep = "")


# fit a meta-analytic model to estimates of Cronbach's Alpha, for the control groups

neg_antc_rel_0 <- metafor::rma(yi = log(1 - alpha0),
                               sei = (2*j) / ((j - 1)*(n0 - 1)),
                               measure = "GEN",
                               method = "REML",
                               data = neg_antc_SMDs.df)

sqrt(var_Bonnett_backtransformed(neg_antc_rel_0$b[1], neg_antc_rel_0$tau2))
mean_Bonnett_backtransformed(neg_antc_rel_0$b[1], neg_antc_rel_0$tau2)


# for the treatment groups

neg_antc_rel_1 <- metafor::rma(yi = log(1 - alpha1),
                               sei = (2*j) / ((j - 1)*(n1 - 1)),
                               measure = "GEN",
                               method = "REML",
                               data = neg_antc_SMDs.df)

sqrt(var_Bonnett_backtransformed(neg_antc_rel_1$b[1], neg_antc_rel_1$tau2))
mean_Bonnett_backtransformed(neg_antc_rel_1$b[1], neg_antc_rel_1$tau2)





write.csv(pos_SMDs.df, file = here("Data/Processed/SMDs_positive.csv"), row.names = FALSE)
write.csv(neg_SMDs.df, file = here("Data/Processed/SMDs_negative.csv"), row.names = FALSE)
write.csv(pos_state_SMDs.df, file = here("Data/Processed/SMDs_state_positive.csv"), row.names = FALSE)
write.csv(neg_state_SMDs.df, file = here("Data/Processed/SMDs_state_negative.csv"), row.names = FALSE)
write.csv(pos_antc_SMDs.df, file = here("Data/Processed/SMDs_antc_positive.csv"), row.names = FALSE)
write.csv(neg_antc_SMDs.df, file = here("Data/Processed/SMDs_antc_negative.csv"), row.names = FALSE)
