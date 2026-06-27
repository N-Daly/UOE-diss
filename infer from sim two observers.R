rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19)

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

########### PDFs for detection functions

hn <- function(distance, sigma){  exp(-0.5 * (distance/sigma)^2 )  }
log_hn <-  function(distance, sigma){ -0.5 * (distance/sigma)^2 }

hn_tau <- function(distance, tau){ exp( -.5*tau*(distance^2) )  }
log_hn_tau <- function(distance, tau){  -.5*tau*(distance^2) }

tau_to_sigma <- function(s){ s^-.5 }
sigma_to_tau <- function(t){ t^-2 }

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



####### the PDFs for detection at a given distance

detect_func_2observer_sigma <- function(distance, sigmaA, sigmaB){
  #ln( (1-pA)(1-pB) )
  log_terms <- log1p( -hn(distance, sigmaA) ) + log1p( -hn(distance, sigmaB) )
  # cant think of a better way to do this bit
  #  1-terms  =  1- e^log_terms =  -1*( e^logterms - 1 ) 
  # my assumption is that if the probs are high this will be stable for the small log terms
  # and if the probs are low this will be also fine for the bigger log terms
  -expm1(log_terms)
}


########### simulate the ground truth

true_sigmaA <- 4; true_sigmaB <- 2
set.seed(888)

# simulate a ground truth
ground_truth <- simulate_lcgp_distance_thinning(
  detect_func = hn,
  detect_func_paramA = true_sigmaA,
  detect_func_paramB = true_sigmaB,
  true_beta0 = -5,
  true_rho = 500, true_sigma_GRF=1
)

#the sampled points are in dd
dd <- ground_truth$samples_df

#preprocessing
dd$detected <- 3
dd$detected[ dd$detectA & !dd$detectB ] = 1
dd$detected[ !dd$detectA & dd$detectB ] = 2

overall_lambda <- ground_truth$overall_lambda
nrow(dd); ground_truth$true_abundance;
summary(dd$distance)



############## modelling

######## set up before modelling
matern_prior <- inla.spde2.pcmatern(
  ground_truth$the_mesh,
  prior.range = c(600, 0.1), # true rho is 500
  prior.sigma = c(.5, 0.5)   # true sigma is 1
)


# # Model just what observer A saw
# fit_observer_A <- model_one_observer_hn(dd[dd$detectA, ], ground_truth, matern_prior)
# # Model just what observer B saw
# fit_observer_B <- model_one_observer_hn(dd[dd$detectB, ], ground_truth, matern_prior)

# Model what A and B saw as a combined observer
fit_merged_observers <- model_one_observer_hn(dd, ground_truth, matern_prior)
# Model what A and B saw as a two observer likelihood
fit_two_observers <- model_two_observers_hn(dd, ground_truth, matern_prior)

dists <- seq(1,8, length.out=1000)



pred_detect_two_obs <- predict(
  fit_two_observers,
  formula = ~ detect_func_2observer_sigma(dists, sigmaA, sigmaB)
)

pred_detect_merged_obs <- predict(
  fit_merged_observers,
  formula = ~  hn(dists, sigma)
)

# pred_detect_obs_A <- predict(
#   fit_observer_A,
#   formula = ~ hn(dists, sigma)
# )

# predict_detect_obs_B <- predict(
#   fit_observer_B,
#   formula = ~ hn(dists, sigma)
# )

res <- rbind(
  cbind(pred_detect_two_obs, model="two obs"),
  cbind(pred_detect_merged_obs, model="merged") 
  )
res$distance <- rep(dists, 2)
res$Truth <- rep(detect_func_2observer_sigma(dists, true_sigmaA, true_sigmaB), 2)


ggplot(res) +
  geom_ribbon(aes(x=distance, ymin=q0.025, max=q0.975, fill=model),  alpha=.3 ) +
  geom_line(aes(x=distance, y=mean,color=model), linewidth=1.5) +
  geom_line(aes(distance, Truth), linewidth=1.5) +
  labs(y="detection probability", title = "Estimates of the detection function \t Black = Truth") +
  ylim(0,1)











