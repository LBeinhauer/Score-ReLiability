




# generate estimates of ln-observed score variance for each sample & projects
varX_est.L <- lapply(which(effect_index)[nn_eff_idx], FUN = function(x){
  
  # sometimes calculation of SE might lead to issues, as negative variances can not be log-transformed
  #  therefore, function is run within a tryCatch-environment, so the script does not crash
  # apply_Bootstrap_SE_Project.specific is the project-specific function
  tryCatch(apply_Bootstrap_SE_Project.specific(na.omit(data.list[[x]]), component = "X", R = 3000),
           
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


# generate estimates of ln-observed score variance for each sample & projects
varE_est.L <- lapply(which(effect_index)[nn_eff_idx], FUN = function(x){
  
  # sometimes calculation of SE might lead to issues, as negative variances can not be log-transformed
  #  therefore, function is run within a tryCatch-environment, so the script does not crash
  # apply_Bootstrap_SE_Project.specific is the project-specific function
  tryCatch(apply_Bootstrap_SE_Project.specific(na.omit(data.list[[x]]), component = "E", R = 3000),
           
           # print error to console
           error = function(e)(cat("ERROR: ", conditionMessage(e), " - ",
                                   substr(names(data.list), 
                                          (regexpr("Project) Data/", names(data.list)) + 14), 
                                          (nchar(names(data.list))-4))[x], 
                                   " - ", x, "\n")))
})




# Perform random-effects meta-analysis on estimates of ln-observed score variance using metafor
varE_rma.list <- lapply(seq_along(varE_est.L), FUN = function(x){
  
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
  
  data.frame(mu_X = bt_var_m(varX_rma.list[[x]]),
             tau2_X = bt_var_v(varX_rma.list[[x]]))
  
})









apply_Bootstrap_SE_Project.specific2 <- function(data, R = 3000, component = "E"){
  
  if(component == "E"){
    stat.f <- bootstrap_SE_varE
  }
  if(component == "T"){
    stat.f <- bootstrap_SE_varT
  }
  if(component == "X"){
    stat.f <- bootstrap_SE_varX
  }
  
  # suppress messages, as coefficientalpha-package can be quite "noisy"
  suppressMessages(
    
    # use apply loop, to loop over replications of a single effect (indicated by the "source"-variable)
    df <- apply(as.matrix(seq_along(unique(data$source))), MARGIN = 1, FUN = function(x){
      
      # use boot-package to generate bootstrap samples of true or error score variance, within a single replication
      bvar <- boot(data = na.omit(data[data$source == unique(data$source)[x],-grep("source", names(data))]),
                   statistic = stat.f,
                   R = R)
      
      # remove "source"-variable from data-set, to compute an estimate of Cronbach's Alpha
      d <- data[data$source == unique(data$source)[x],-grep("source", names(data))]
      
      d1 <- d[d$group == 1, -(which(names(d) %in% c("source", "group")))]
      d0 <- d[d$group == 0, -(which(names(d) %in% c("source", "group")))]
      
      mv1 <- rowMeans(d1)
      mv0 <- rowMeans(d0)
      
      sd1 <- sqrt(mean((mv1 - mean(mv1))^2))
      sd0 <- sqrt(mean((mv0 - mean(mv0))^2))
      
      n1 <- length(d1)
      n0 <- length(d0)
      
      C1 <- cov(d1)
      j1 <- dim(C1)[1]
      alpha1 <- (1 - sum(diag(C1))/sum(C1)) * (j1/(j1 - 1))
      
      # compute Cronbach's Alpha
      C0 <- cov(d0)
      j0 <- dim(C0)[1]
      alpha0 <- (1 - sum(diag(C0))/sum(C0)) * (j0/(j0 - 1))
      
      alpha_pooled <- (n1 / (n0 + n1)) * alpha1 + (n0 / (n0 + n1)) * alpha0
      
      # pooled standard deviation
      var_X <- ((n1-1)*sd1^2 + (n0-1)*sd0^2)/(n1+n0-2)
      
      # estimate true score variance sigma^2_T
      var_T <- as.numeric(alpha_pooled * var_X)
      
      # estimate errorcsore variance sigma^2_E
      var_E <- var_X - var_T
      
      if(component == "E"){
        var_res <- var_E
      }
      if(component == "X"){
        var_res <- var_X
      }
      if(component == "T"){
        var_res <- var_T
      }
      
      
      # return a data.frame, containing bootstrapped estimates of the standard error of the log-transformed
      #  variance component, as well as the bootstrapped mean of the log-transformed variance component. 
      # additionaly, return the empirical estimate using the full replication sample
      return(data.frame(SE = sd(log(sqrt(bvar$t))), 
                        boot.mean = mean(log(bvar$t)),
                        var.emp = log(sqrt(var_res))))
    })
  ) # output of this apply-loop is a List, each list-element consisting of estimates for a single replication
  
  # transform the apply-loop-list into a single data.frame, containing estimates of:
  # - log-variance-component's bootstrapped standard error
  # - log-variance-component's bootstrapped mean
  # - empirical estimate of variance component, untransformed!
  # - the respective replication author.
  df.formatted <- data.frame(SE = sapply(df, FUN = function(x){x$SE}),
                             boot.mean = sapply(df, FUN = function(x){x$boot.mean}),
                             var.est = sapply(df, FUN = function(x){x$var.emp}),
                             source = unique(data$source))
  
}




# generate estimates of ln-observed score variance for each sample & projects
sdX_est.L <- lapply(which(effect_index)[nn_eff_idx], FUN = function(x){
  
  # sometimes calculation of SE might lead to issues, as negative variances can not be log-transformed
  #  therefore, function is run within a tryCatch-environment, so the script does not crash
  # apply_Bootstrap_SE_Project.specific is the project-specific function
  tryCatch(apply_Bootstrap_SE_Project.specific2(na.omit(data.list[[x]]), component = "X", R = 3000),
           
           # print error to console
           error = function(e)(cat("ERROR: ", conditionMessage(e), " - ",
                                   substr(names(data.list), 
                                          (regexpr("Project) Data/", names(data.list)) + 14), 
                                          (nchar(names(data.list))-4))[x], 
                                   " - ", x, "\n")))
})




# Perform random-effects meta-analysis on estimates of ln-observed score variance using metafor
sdX_rma.list <- lapply(seq_along(sdX_est.L), FUN = function(x){
  
  # at times estimates of ln-varX & SE may be NA, as negative estimates can't be transformed
  #  therefore, function is nested in tryCatch, so it doesn't break down with errors
  tryCatch(metafor::rma(measure = "GEN", method = "REML", 
                        yi = sdX_est.L[[x]]$var.est, 
                        sei = sdX_est.L[[x]]$SE),
           
           # print error to console
           error = function(e)(cat("ERROR: ", conditionMessage(e), " - ",
                                   substr(names(data.list), 
                                          (regexpr("Project) Data/", names(data.list)) + 14), 
                                          (nchar(names(data.list))-4))[x], 
                                   " - ", x, "\n")))
  
})


sds_rma_L <- lapply(1:length(varE_rma.list), FUN = function(x){
  
  data.frame(mu_X = bt_var_m(sdX_rma.list[[x]]),
             tau2_X = bt_var_v(sdX_rma.list[[x]]))
  
})


sds <- do.call(rbind, sds_rma_L)
vars <- do.call(rbind, vars_rma_L)

plot(cbind(sqrt(vars$mu_X), sds$mu_X))
abline(a = 0, b = 1)


plot(cbind(.5*(sqrt(vars$tau2_X) / vars$mu_X), (sqrt(sds$tau2_X) / sds$mu_X)))
abline(a = 0, b = 1)


test_df <- data.frame(unique(ES_rma_df$MASC)[nn_eff_idx],
                      sds,
                      delta_mu_X = sqrt(vars$mu_X),
                      delta_tau2_X = vars$tau2_X / (4 * vars$mu_X)) %>%
  mutate_if(is.numeric, round, 5)

test_df

write.csv(test_df, here("Tables/delta_test.csv"), row.names = FALSE)

round(test_df, 5)




liability_inequality <- function(tau_MD, mu_MD, mu_sigma2x, mu_sigma2t, tau_sigma2x, tau_sigma2t){
  
  R1 <- mu_sigma2t/mu_sigma2x
  R2 <- (tau_sigma2t^2)/(tau_sigma2x^2)
  
  ineq <- ((tau_MD^2 )* ((1/R1) - 1)) + (.25 * (mu_MD^2) * ((tau_sigma2x/mu_sigma2x)^2) * ((R2 / (R1^3)) - 1))
  
  return(data.frame(R1 = R1,
                    R2 = R2,
                    inequality = ineq,
                    tau_MD = tau_MD,
                    mu_MD = mu_MD,
                    mu_sigma2x = mu_sigma2x,
                    mu_sigma2t = mu_sigma2t,
                    tau_sigma2x = tau_sigma2x, 
                    tau_sigma2t = tau_sigma2t,
                    CV_sigma2x = tau_sigma2x/mu_sigma2x))
  
}



input_df <- data.frame(tau_MD = .1, 
                       mu_MD = 1, 
                       mu_sigma2x = 1, 
                       mu_sigma2t = c(.5, .5, .5, .7, .7, .7, .9, .9, .9), 
                       tau_sigma2x = .3, 
                       tau_sigma2t = c(.3, .2, .1, .3, .2, .1, .3, .2, .1))

test <- liability_inequality(input_df$tau_MD, input_df$mu_MD, input_df$mu_sigma2x, input_df$mu_sigma2t,
                     input_df$tau_sigma2x, input_df$tau_sigma2t)

test


R1 <- seq(from = .01, to = 1, length.out = 50)
R2 <- seq(from = .01, to = 1, length.out = 50)

R_grid <- expand.grid(R1, R2)
names(R_grid) <- c("R1", "R2")

tau_MD <- .1
mu_MD <- 1
mu_sigma2x <- 1
tau_sigma2x <- .1

ineq_func <- function(R_grid, tau_MD, mu_MD, mu_sigma2x, tau_sigma2x){
  
  val <- ((tau_MD^2 )* ((1/R_grid$R1) - 1)) + (.25 * (mu_MD^2) * ((tau_sigma2x/mu_sigma2x)^2) * ((R_grid$R2 / (R_grid$R1^3)) - 1))
  
  data.frame(R1 = R_grid$R1,
             R2 = R_grid$R2,
             val = val)
}

input_df <- ineq_func(R_grid, tau_MD, mu_MD, mu_sigma2x, tau_sigma2x)

ggplot(input_df) +
  geom_raster(aes(x = R1, y = R2, fill = val))
