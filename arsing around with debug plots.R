rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19);options(digits=3);invisible(gc())

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)
library(patchwork)
library(dplyr)


my_dir <- r"(C:\Users\ND\OneDrive - University of Edinburgh\Dissertation\UOE-diss)"
setwd(my_dir)
source("simulate 2D ground truth.R")
source("observer models.R")
source("function Definitions.R")
source("mesh construction methods.R")

########### simulate the ground truth

true_sigmaA <- 1.75; true_sigmaB <- 2; true_gammaA <- 1; true_gammaB <- 4
# true_sigmaA <- 2; true_sigmaB <- 4
true_range <- 500; true_sigma_grf <- 1

dists <- seq(0,8, length.out=500)
pA <- hr(dists, sigma = true_sigmaA, gamma = true_gammaA)
pB <- hr(dists, sigma = true_sigmaB, gamma = true_gammaB)
# pA <- hn(dists, sigma = true_sigmaA)
# pB <- hn(dists, sigma = true_sigmaB)
pany <- 1-(1-pA)*(1-pB)

# distribution of observable distances - uniform over line transect segments
pdf_dist <- 1#/length(dists)
true_detect_prob_df <- data.frame(
  distance = dists,
  Amarginal = pA,
  Bmarginal = pB,
  "detected1" = (pA*(1-pB)),
  "detected2" = ((1-pA)*pB),
  "detected3" = pA*pB,
  "any" = pany
)



set.seed(8010)
# simulate a ground truth
sim_info <- simulate_lcgp_dual_obs_HR_thinning(
  true_sigmaA, true_sigmaB, true_gammaA, true_gammaB,
  true_beta0 = -4,
)
# sim_info <- simulate_lcgp_distance_thinning(
#   detect_func = hn,
#   detect_func_paramA = true_sigmaA, detect_func_paramB = true_sigmaB,
#   true_beta0 = -4
# )

source("function Definitions.R")
lets_have_a_look_at_you(sim_info, true_detect_prob_df)

# every <- sim_info$unthinned_samples_df
# hist(
#   obs$distance[obs$detected==1],
#   30,
#   freq=F,
#   ylim=0:1
# )
# 
# par(mfrow=c(2,1))
# hist(
#   obs$distance[obs$detectA],
#   30,
#   freq=F,
#   ylim=0:1
# )
# lines(dists, pdf_dist*pA)
# lines(density(obs$distance[obs$detectA], from = 0, to = 8))

# par(old_par)
# 
# d <- density(obs$distance[obs$detectA])
# max(d$y)
# plot(
#   density(obs$distance[obs$detectA], from = 0, to = 8),
#   ylim = 0:1
#   
# )
# lines(dists, pA*pdf_dist)
# 
# hist(obs$distance,30, ylim=0:1, freq=F);abline(h=2/8)
# 
# 


# plot(
#   dists, 0*dists,
#   type = "n",
#   ylim = 0:1
# )
# lines(dists, pdf_dist*pA/pany)



