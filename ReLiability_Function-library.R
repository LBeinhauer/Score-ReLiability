

## function library





# function to estimate heterogeneity in Cronbach's Alpha, if transformation was used
var_Bonnett_backtransformed <- function(rma_obj){
  (((-exp(rma_obj$b[1]))^2) * rma_obj$tau2) + (.5*((-exp(rma_obj$b[1]))^2)*(rma_obj$tau2^2)) + ((-exp(rma_obj$b[1])) * (-exp(rma_obj$b[1])) * (rma_obj$tau2^2))
}

# function to estimate mean value of Cronbach's Alpha, if transformation was used
mean_Bonnett_backtransformed <- function(rma_obj){
  1 - exp(rma_obj$b[1]) + ((-exp(rma_obj$b[1])) / 2) * rma_obj$tau2
}



bt_var_m <- function(rma_obj){
  exp(rma_obj$b[1] + (.5*rma_obj$tau2))
}

bt_var_v <- function(rma_obj){
  (exp(rma_obj$tau2) - 1) * exp((2 * rma_obj$b[1]) + rma_obj$tau2)
}

bt_var_m2 <- function(rma_obj){
  exp(.5*rma_obj$b[1] + (.5*.25*rma_obj$tau2))
}

bt_var_v2 <- function(rma_obj){
  (exp(.5*rma_obj$tau2) - 1) * exp((2 * .5 * rma_obj$b[1]) + .25*rma_obj$tau2)
}





## a function required for the boot-package, which allows to estimate the error score
##  variance in a bootstrapped sample
bootstrap_SE_varE <- function(data, indices){
  
  # select subset of data.
  d <- data[indices,]
  
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
  
  return(sqrt(var_E))
  
}

## a function required for the boot-package, which allows to estimate the error score
##  variance in a bootstrapped sample
bootstrap_SE_varX <- function(data, indices){
  
  # select subset of data.
  d <- data[indices,]
  
  d1 <- d[d$group == 1, -(which(names(d) %in% c("source", "group")))]
  d0 <- d[d$group == 0, -(which(names(d) %in% c("source", "group")))]
  
  mv1 <- rowMeans(d1)
  mv0 <- rowMeans(d0)
  
  sd1 <- sqrt(mean((mv1 - mean(mv1))^2))
  sd0 <- sqrt(mean((mv0 - mean(mv0))^2))
  
  n1 <- length(d1)
  n0 <- length(d0)

  # pooled standard deviation
  var_X <- ((n1-1)*sd1^2 + (n0-1)*sd0^2)/(n1+n0-2)
  
  return(var_X)
  
}


## a function required for the boot-package, which allows to estimate the true score
##  variance in a bootstrapped sample
bootstrap_SE_varT <- function(data, indices){
  
  # select subset of data.
  d <- data[indices,]
  
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
  
  return(var_T)
  
}


# Function to estimate bootstrapped Standard Errors for either true or error score variance component,
#  specifically made for empirical data (e.g. containing missing data)
# expects data in a data.frame
apply_Bootstrap_SE_Project.specific <- function(data, R = 100, component = "E"){
  
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
      return(data.frame(SE = sd(log(bvar$t)), 
                        boot.mean = mean(log(bvar$t)),
                        var.emp = log(var_res)))
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
