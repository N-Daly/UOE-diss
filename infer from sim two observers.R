rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19);options(digits=3);invisible(gc())

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)
library(patchwork)

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

spline_detect_func <- function(spline_effect){ exp(-spline_effect) }
ln_spline_detect_func <- function(spline_effect){ -spline_effect }

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
  models,
  ips,
  true_detect,
  true_loglambda
){
  base_mod <- models[[1]] # should always be the one observer hn
  
  # get the scores of this base model
  base_mod$DS <- dawid_sebastiani_score(base_mod$pred_loglambda, true_loglambda)
  
  base_mod$loglambda_SE <- (base_mod$pred_loglambda$mean - true_loglambda)^2
  
  base_mod$lambda_AE <- abs(base_mod$pred_lambda$median - true_loglambda)
  
  base_mod$detect_AE <- abs(base_mod$pred_detect$median - true_detect)
  
  score_diffs <- NULL
  for (i in 2:length(models)){
    mod <- models[[i]]
    
    #get this model' scores 
    mod$DS <- dawid_sebastiani_score(mod$pred_loglambda, true_loglambda)
    mod$loglambda_SE <- (mod$pred_loglambda$mean - true_loglambda)^2
    mod$lambda_AE <- abs(mod$pred_lambda$median - true_loglambda)
    mod$detect_AE <- abs(mod$pred_detect$median - true_detect)
    
    # now take the difference in scores between this and the base
    # and integrate it over the domain
    mod_diff <- list(
      DS = sum(ips$weight * (base_mod$DS - mod$DS)),
      loglambda_SE = sum(ips$weight * (base_mod$loglambda_SE - mod$loglambda_SE)),
      lambda_AE = sum(ips$weight * (base_mod$lambda_AE - mod$lambda_AE)),
      detect_AE = mean(base_mod$detect_AE - mod$detect_AE)
    )
    
    mod_diff$model <- mod$name
    
    score_diffs <- rbind(score_diffs, mod_diff)
  }
  
  as.data.frame(score_diffs, row.names = F)
}



########### simulate the ground truth

true_sigmaA <- 4; true_sigmaB <- 2; true_range <- 500; true_sigma_grf <- 1

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

############## modelling

######## set up before modelling
matern_prior <- inla.spde2.pcmatern(
  mexdolphin_sf$mesh,
  prior.range = c(600, 0.1), # true rho is 500
  prior.sigma = c(.5, 0.5)   # true sigma is 1
)

detect_mesh <- fm_mesh_1d(
  loc = c(0,2,4,6,8), 
  boundary = c("dirichlet", "free"),
  degree = 1
)

# alpha=2, rho=2, sigma=.25
detect_matern <- inla.spde2.pcmatern(
  detect_mesh,
  alpha = 2,
  prior.range = c(2, 0.99), # P(rho < val) = alpha
  prior.sigma = c(0.25, 0.01) # P(sigma > val) = alpha
)

dists <- seq(0,8, length.out=1000)

ips <- st_as_sf( readRDS("ips_interiorTransects_2subdivisions.rda") )
# so the one observer models have only a geometry dimension
# the two observer model needs a detection state dimension - who spotted the animal
# here i naively expand the existing ips uniformly across this dimension
# fingers crossed
ips_with_detection_states <- ips[rep(1:nrow(ips), each=3),]
ips_with_detection_states$detected <- rep(1:3, times=nrow(ips))
ips_with_detection_states$weight = ips_with_detection_states$weight/3


# for model scoring later
true_detect_prob = detect_func_2observer_sigma(dists, true_sigmaA, true_sigmaB)
true_loglambda_at_ip <- c( fm_evaluate(sim_info$the_mesh, field = sim_info$log_lambda, loc = ips$geometry) )
list_of_models <- list()


######## fitting models

# Model what A and B saw as a combined observer with a hn detection function
catt("fitting merged")
start <- proc.time()
fit_merged_observers <- model_one_observer_hn(observed_points, sim_info, matern_prior, ips=ips)
end <- proc.time()
print((end-start)[3])

start <- proc.time()
list_of_models[[1]] <- get_preds_from_one_observer_hn_fit(fit_merged_observers)
end <- proc.time()
print((end-start)[3])


# Model what A and B saw as a two observer likelihood with hn detection functions
catt("fitting two")
start <- proc.time()
fit_two_observers <- model_two_observers_hn(
  observed_points, sim_info, matern_prior, ips = ips_with_detection_states
)
end <- proc.time()
print((end-start)[3])

start <- proc.time()
list_of_models[[2]] <- get_preds_from_two_observers_hn_fit(fit_two_observers)
end <- proc.time()
print((end-start)[3])


# Model what A and B saw as a combined observer with a spline(like) detection function
catt("fitting spline")
start <- proc.time()
fit_one_spline <- model_one_observer_spline(
  observed_points, sim_info, matern_prior, detect_matern, ips=ips
)
end <- proc.time()
print((end-start)[3])

start <- proc.time()
list_of_models[[3]] <- get_preds_from_one_observer_spline_fit(fit_one_spline)
end <- proc.time()
print((end-start)[3])


initial_results <- get_scoring_differences(
  list_of_models,
  ips,
  true_detect_prob,
  true_loglambda_at_ip
)
initial_results

############## Repeated simulations
set.seed(1234)
nsims <- 2
results <- NULL
# for saving the results
timestamp <- format(Sys.Date())
file_name <- paste(timestamp, "simulation results.rda")

for (i in 1:nsims){
  catt("Simulation", i, "of", nsims)
  
  # need to save space and these model fits are large
  rm(list_of_models, fit_merged_observers, fit_two_observers, fit_one_spline)
  invisible(gc())

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
  list_of_models = list()

  # Model what A and B saw as a combined observer with a hn detection function
  catt("fitting merged")
  start <- proc.time()
  fit_merged_observers <- model_one_observer_hn(observed_points, sim_info, matern_prior)
  end <- proc.time()
  print((end-start)[3])
  
  start <- proc.time()
  list_of_models[[1]] <- get_preds_from_one_observer_hn_fit(fit_merged_observers)
  end <- proc.time()
  print((end-start)[3])
  
  
  # Model what A and B saw as a two observer likelihood with hn detection functions
  catt("fitting two")
  start <- proc.time()
  fit_two_observers <- model_two_observers_hn(observed_points, sim_info, matern_prior)
  end <- proc.time()
  print((end-start)[3])
  
  start <- proc.time()
  list_of_models[[2]] <- get_preds_from_two_observers_hn_fit(fit_two_observers)
  end <- proc.time()
  print((end-start)[3])
  
  
  # Model what A and B saw as a combined observer with a spline(like) detection function
  catt("fitting spline")
  start <- proc.time()
  fit_one_spline <- model_one_observer_spline(
    observed_points, sim_info, matern_prior, detect_matern, ips=ips
  )
  end <- proc.time()
  print((end-start)[3])
  
  start <- proc.time()
  list_of_models[[3]] <- get_preds_from_one_observer_spline_fit(fit_one_spline)
  end <- proc.time()
  print((end-start)[3])
  
  catt("comparing the models' scores")
  start <- proc.time()
  some_results <- get_scoring_differences(
    list_of_models,
    ips,
    true_detect_prob,
    true_loglambda_at_ip
  )
  end <- proc.time()
  print((end-start)[3])
  
  #save results early just in case
  results <- rbind(results, some_results)
  # wrangling with classes
  results$model <- as.character(results$model)
  results[-5] <- apply(results[-5], 2, as.numeric)
  
  saveRDS(results, file=file_name)
}



# plot the difference in scores, a difference of zero implies the models perform the same wrt to that score
# all scores are negatively orientated so a positive difference means the simpler merged observer model
# performed worse

g <- ggplot(results, aes(fill=model)) + expand_limits(y=0) +  geom_abline(intercept = 0, slope = 0, colour="red", linewidth=2)

dsp <-g + geom_boxplot(aes(y=DS)) + labs(title="integrated DS on log lambda")
msep <- g + geom_boxplot(aes(y=loglambda_SE)) + labs(title="integrated SE on mean log lambda")
maep <- g + geom_boxplot(aes(y=lambda_AE)) + labs(title = "integrated AE on median lambda")
maedetectp <- g + geom_boxplot(aes(y=detect_AE)) + labs(title="mean AE on detection prob")

(dsp + msep) / ( maep + maedetectp)



