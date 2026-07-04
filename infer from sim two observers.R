rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19);options(digits=3)

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)
library("patchwork")

my_dir <- r"(C:\Users\ND\OneDrive - University of Edinburgh\Dissertation\UOE-diss)"
setwd(my_dir)
source("simulate 2D ground truth.R")
source("observer models.R")

catt <- function(...){ cat(..., "\n")}# pet peeve

###########  detection functions given

hn <- function(distance, sigma){  exp(-0.5 * (distance/sigma)^2 )  }
log_hn <-  function(distance, sigma){ -0.5 * (distance/sigma)^2 }


######## the PDFs for both detection state and distance

log_g_2observer <- function(distance, detected, sigmaA, sigmaB, eps=1e-6){
  
  log_terms <- numeric(length(detected))
  # A alone
  ii <- detected==1 
  log_terms[ii] = log_hn(distance[ii], sigmaA) + log1p(-hn(distance[ii], sigmaB)*(1-eps) )
  # B alone
  ii <- detected==2
  log_terms[ii] = log1p(-hn(distance[ii], sigmaA)*(1-eps) ) + log_hn(distance[ii], sigmaB)
  # A,B
  ii <- detected==3
  log_terms[ii] = log_hn(distance[ii], sigmaA) + log_hn(distance[ii], sigmaB)
  
  log_terms
}



####### the probability of detection at a given distance

detect_func_2observer_sigma <- function(distance, sigmaA, sigmaB){
  #ln( (1-pA)(1-pB) )
  log_terms <- log1p( -hn(distance, sigmaA) ) + log1p( -hn(distance, sigmaB) )
  # cant think of a better way to do this bit
  #  1-terms  =  1- e^log_terms =  -1*( e^logterms - 1 ) 
  # my assumption is that if the probs are high this will be stable for the small log terms
  # and if the probs are low this will be also fine for the bigger log terms
  -expm1(log_terms)
}

######## functions for assessing model fit 

get_spatial_param_width <- function(mdl, param, correct){
  # checks if the model's 95% credible interval contains the true value of a hyperparam
  # if so it returns 1 and the interval width
  cred_int <- mdl$summary.hyperpar[param, c("0.025quant", "0.975quant")]
  cred_int <- unlist(cred_int)
  
  if ( cred_int[1] <= correct & correct <= cred_int[2] ){
    c(1, diff(cred_int) )
  } else {
    c(0, diff(cred_int) )
  }
}

dawid_sebastiani_score <- function(post_pred, true_value){
  E <- post_pred$mean
  V <- post_pred$mean.mc_std_err^2
  ( (true_value - E)^2 / V)  + log(V)
}


compare_merged_vs_two_observer_models <- function(
    merged, two_obs,
    true_sigmaA, true_sigmaB, dist_res,
    do_plot = F,
    true_spatial_range=true_range, true_spatial_stdv = true_sigma_grf
    ){
  
  output <- matrix(NA, 2, ncol=1+4+1)
  
  
  true_detect <- detect_func_2observer_sigma(dist_res, true_sigmaA, true_sigmaB)
  
  output[1:2, i<-1] <- c("merged", "two obs")
  
  output[1, 2:3] <- get_spatial_param_width(merged, "Range", true_spatial_range)
  output[1, 4:5] <- get_spatial_param_width(merged, "Stdev", true_spatial_stdv)
  
  pred_detect_merged_obs <- predict(merged, formula = ~ hn(dist_res, sigma))
  output[1, 6] <- mean(abs( pred_detect_merged_obs$mean - true_detect ))
  
  
  output[2, 2:3] <- get_spatial_param_width(two_obs, "Range", true_spatial_range)
  output[2, 4:5] <- get_spatial_param_width(two_obs, "Stdev", true_spatial_stdv)
  
  pred_detect_two_obs <- predict(two_obs, formula = ~ detect_func_2observer_sigma(dist_res, sigmaA, sigmaB) )
  output[2, 6] <- mean(abs( pred_detect_two_obs$mean - true_detect ))
  
  if( do_plot){
    res <- rbind(
      cbind(pred_detect_two_obs, model="two obs"),
      cbind(pred_detect_merged_obs, model="merged") 
    )
    res$distance <- rep(dist_res, 2)
    res$Truth <- rep(true_detect, 2)
    
    g <- ggplot(res) +
      geom_ribbon(aes(x=distance, ymin=q0.025, max=q0.975, fill=model),  alpha=.3 ) +
      geom_line(aes(x=distance, y=mean,color=model), linewidth=1.5) +
      geom_line(aes(distance, Truth), linewidth=1.5) +
      labs(y="detection probability", title = "Estimates of the detection function \t Black = Truth") +
      ylim(0,1)
    
    print(g)
  }
  
  output
}

########### simulate the ground truth

true_sigmaA <- 4; true_sigmaB <- 2; true_range <- 500; true_sigma_grf <- 1


############## modelling

######## set up before modelling
matern_prior <- inla.spde2.pcmatern(
  mexdolphin_sf$mesh,
  prior.range = c(600, 0.1), # true rho is 500
  prior.sigma = c(.5, 0.5)   # true sigma is 1
)

dists <- seq(0,8, length.out=1000)


set.seed(800)
# simulate a ground truth
sim_info <- simulate_lcgp_distance_thinning(
  detect_func = hn,
  detect_func_paramA = true_sigmaA,
  detect_func_paramB = true_sigmaB,
  true_beta0 = -5,
  true_rho = true_range, true_sigma_GRF=true_sigma_grf
)

#the sampled points are in dd
dd <- sim_info$samples_df

# Model what A and B saw as a combined observer
catt("fitting merged")
start <- Sys.time()
fit_merged_observers <- model_one_observer_hn(dd, sim_info, matern_prior)
end <- Sys.time()
catt(end-start)

# Model what A and B saw as a two observer likelihood
catt("fitting two")
start <- Sys.time()
fit_two_observers <- model_two_observers_hn(dd, sim_info, matern_prior)
end <- Sys.time()
catt(end-start)


get_DS_score_difference <- function(
    sim_info, fit_merged_observers, fit_two_observers
    ){
  
  #both models predict log lambda
  merged_observers_pred <- predict(
    fit_merged_observers, sim_info$interior_vertices,
    formula = ~ {
      list(
        lambda=exp(Intercept + spde),
        log_lambda = Intercept + spde
      )
    },
    n.samples=2000
  )
  two_observers_pred <- predict(
    fit_two_observers, sim_info$interior_vertices,
    formula = ~ {
      list(
        lambda=exp(Intercept + spde),
        log_lambda = Intercept + spde
      )
    },
    n.samples=2000
  )
  # DS score at vertices
  merged_observers_DS <- dawid_sebastiani_score(merged_observers_pred$log_lambda, sim_info$log_lambda)
  two_observers_DS <- dawid_sebastiani_score(two_observers_pred$log_lambda, sim_info$log_lambda)
  
  length(merged_observers_DS); length(two_observers_DS);length(sim_info$log_lambda)
  # get the difference in DS scores relative to a single merged observer model
  DS_score_difference <- merged_observers_DS - two_observers_DS
  
  # integrated DS score difference over the study region
  integrated_DS_difference <- sum(sim_info$interior_ips$weight * DS_score_difference)
  
  integrated_DS_difference
}


get_DS_score_difference(sim_info, fit_merged_observers, fit_two_observers)

# set.seed(1234)
# nsims <- 10
# results <- matrix(NA, nrow=2*nsims, ncol=6)
# colnames(results) <- c("model", "range_correct", "range_CI_width", "stdv_correct", "stdv_CI_width", "MAE_detection")
# 
# for (i in 1:nsims){
#   catt("Simulation", i, "of", nsims)
#   
#   # simulate a ground truth
#   sim_info <- simulate_lcgp_distance_thinning(
#     detect_func = hn,
#     detect_func_paramA = true_sigmaA,
#     detect_func_paramB = true_sigmaB,
#     true_beta0 = -5,
#     true_rho = true_range, true_sigma_GRF=true_sigma_grf
#   )
#   
#   #the sampled points are in dd
#   dd <- sim_info$samples_df
#   
#   # Model what A and B saw as a combined observer
#   catt("fitting merged")
#   start <- Sys.time()
#   fit_merged_observers <- model_one_observer_hn(dd, sim_info, matern_prior)
#   end <- Sys.time()
#   catt(end-start)
#   
#   # Model what A and B saw as a two observer likelihood
#   catt("fitting two")
#   start <- Sys.time()
#   fit_two_observers <- model_two_observers_hn(dd, sim_info, matern_prior)
#   end <- Sys.time()
#   catt(end-start)
#   
#   catt("comparing the two")
#   results[2*i-1:0,] <- compare_merged_vs_two_observer_models(
#     fit_merged_observers, fit_two_observers,
#     true_sigmaA, true_sigmaB, dists,
#     do_plot = T
#   )
# }
# results
# 
# 
# results <- as.data.frame(results)
# # the model column is a character so so did everything else in the matrix  
# results[-1] <- lapply(results[-1], as.numeric)
# 
# 
# 
# save(results, file="-sims.rda")






