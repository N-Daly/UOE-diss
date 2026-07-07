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


dawid_sebastiani_score <- function(post_pred, true_value){
  E <- post_pred$mean
  # i am not sure if this step is correct
  V <- post_pred$sd^2
  ( (true_value - E)^2 / V)  + log(V)
}


get_scoring_differences <- function(
    sim_info, fit_merged_observers, fit_two_observers,
    true_detect_prob = detect_func_2observer_sigma(dists, true_sigmaA, true_sigmaB),
    dists = seq(0,8, length.out=1000)
){
  #for brevity
  loglambda <- sim_info$log_lambda
  
  #both models predict log lambda, lambda, and the probability of (any) detection
  # squeezing these all into one predict call so it'll hopefully be faster
  merged_observer_pred <- predict(
    fit_merged_observers, 
    newdata = list(
      geometry=sim_info$interior_vertices$geometry,
      dists=dists
    ),
    formula = ~ {
      list(
        lambda = exp(Intercept + spde),
        log_lambda = Intercept + spde,
        detection = hn(dists, sigma)
      )
    },
    n.samples=2000,
  )
  # different detection prob function for two observer model
  two_observers_pred <- predict(
    fit_two_observers, 
    newdata = list(
      geometry=sim_info$interior_vertices$geometry,
      dists=dists
    ),
    formula = ~ {
      list(
        lambda = exp(Intercept + spde),
        log_lambda = Intercept + spde,
        detection = detect_func_2observer_sigma(dists, sigmaA, sigmaB)
      )
    },
    n.samples=2000
  )
  
  # DS scores on log lambda at vertices
  merged_observers_DS <- dawid_sebastiani_score(merged_observer_pred$log_lambda, loglambda)
  two_observers_DS <- dawid_sebastiani_score(two_observers_pred$log_lambda, loglambda)
  # integrated DS score difference over the study region
  integrated_DS_difference <- sum( sim_info$interior_ips$weight * (merged_observers_DS - two_observers_DS) )
  
  # MSE scores on the posterior mean for log lambda
  merged_observers_MSE <- mean(merged_observer_pred$log_lambda$mean - loglambda)
  two_observers_MSE <- mean(two_observers_pred$log_lambda$mean - loglambda)
  
  # MAE scores on the posterior median for lambda
  merged_observers_MAE <- mean(abs( merged_observer_pred$lambda$median - exp(loglambda) )) 
  two_observers_MAE <- mean(abs( two_observers_pred$lambda$median - exp(loglambda) )) 
  
  # MAE errors on the posterior mean detection probability
  merged_observers_detect_MAE <- mean(abs( merged_observer_pred$detection$mean - true_detect_prob ))
  two_observers_detection_MAE <- mean(abs( two_observers_pred$detection$mean - true_detect_prob ))
  
  #report the differences in scores relative to a single merged observer model
  list(
    DS = integrated_DS_difference,
    MSE =  merged_observers_MSE - two_observers_MSE,
    MAE = merged_observers_MAE - two_observers_MAE,
    MAE_detection = merged_observers_detect_MAE - two_observers_detection_MAE
  )
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

#the sampled points are in observed_points
observed_points <- sim_info$samples_df

# Model what A and B saw as a combined observer
catt("fitting merged")
start <- proc.time()
fit_merged_observers <- model_one_observer_hn(observed_points, sim_info, matern_prior)
end <- proc.time()
print((end-start)[3])

# Model what A and B saw as a two observer likelihood
catt("fitting two")
start <- proc.time()
fit_two_observers <- model_two_observers_hn(observed_points, sim_info, matern_prior)
end <- proc.time()
print((end-start)[3])

get_scoring_differences(sim_info, fit_merged_observers, fit_two_observers)



set.seed(1234)
nsims <- 20
results <- NULL

for (i in 1:nsims){
  catt("Simulation", i, "of", nsims)

  # simulate a ground truth
  sim_info <- simulate_lcgp_distance_thinning(
    detect_func = hn,
    detect_func_paramA = true_sigmaA,
    detect_func_paramB = true_sigmaB,
    true_beta0 = -5,
    true_rho = true_range, true_sigma_GRF=true_sigma_grf
  )

  #the sampled points are in observed_points
  observed_points <- sim_info$samples_df

  # Model what A and B saw as a combined observer
  catt("fitting merged")
  start <- proc.time()
  fit_merged_observers <- model_one_observer_hn(observed_points, sim_info, matern_prior)
  end <- proc.time()
  print((end-start)[3])

  # Model what A and B saw as a two observer likelihood
  catt("fitting two")
  start <- proc.time()
  fit_two_observers <- model_two_observers_hn(observed_points, sim_info, matern_prior)
  end <- proc.time()
  print((end-start)[3])

  catt("comparing the two")
  start <- proc.time()
  results <- rbind(
    results, get_scoring_differences(sim_info, fit_merged_observers, fit_two_observers)
  )
  end <- proc.time()
  print((end-start)[3])
}
(results <- as.data.frame( apply(results, 2, as.numeric) ) )

# save the results
saveRDS(results, file="sat-4-07-score-results.rda")


# plot the difference in scores, a difference of zero implies the models perform the same wrt to that score
# all scores are negatively orientated so a positive difference means the simpler merged observer model
# performed worse

g <- ggplot(results) + expand_limits(y=0) +  geom_abline(intercept = 0, slope = 0, colour="red", linewidth=2)

dsp <-g + geom_boxplot(aes(y=DS)) + labs(title="integrated DS on log lambda")
msep <- g + geom_boxplot(aes(y=MSE)) + labs(title="integrated MSE on mean log lambda")
maep <- g + geom_boxplot(aes(y=MAE)) + labs(title = "integrated MAE on median lambda")
maedetectp <- g + geom_boxplot(aes(y=MAE_detection)) + labs(title="MAE on detection prob") 

(dsp + msep) / ( maep + maedetectp)








