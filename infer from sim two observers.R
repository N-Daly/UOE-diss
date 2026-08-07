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

########### state the ground truth

true_sigmaA <- 1.75; true_sigmaB <- 2; true_gammaA <- 1; true_gammaB <- 4
true_range <- 500; true_sigma_grf <- 1

# work out the different probabilities of detection
dists <- seq(0,8, length.out=500)

pA <- hr(dists, sigma = true_sigmaA, gamma = true_gammaA)
pB <- hr(dists, sigma = true_sigmaB, gamma = true_gammaB)
pany = detect_func_2_observer_hr(dists, true_sigmaA, true_sigmaB, true_gammaA, true_gammaB)

true_detect_prob_df <- data.frame(
  distance = dists,
  Amarginal = pA,
  Bmarginal = pB,
  "detected1" = pA*(1-pB),
  "detected2" = (1-pA)*pB,
  "detected3" = pA*pB,
  "any" = pany
)


######## set up before modelling
mesh <- make_spatially_varying_mesh3(60, coarse=60)

ips <- readRDS("ips make_spatially_varying_mesh3 hex edge length 2-30.rda")
ips_with_detection_states <- readRDS("ips_with_detect make_spatially_varying_mesh3 hex edge length 2-30.rda")

matern_prior <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(600, 0.1), # true rho is 500
  prior.sigma = c(.5, 0.5)   # true sigma is 1
)

detect_mesh <- fm_mesh_1d(
  loc = c(0,2,4,6,8), 
  boundary = c("dirichlet", "free"),
  degree = 2
)

# alpha=2; rho=2; sigma=1.5
detect_matern <- inla.spde2.pcmatern(
  detect_mesh,
  alpha = 2,
  prior.range = c(2, 0.8), # P(rho < val) = alpha
  prior.sigma = c(2, 0.01), # P(sigma > val) = alpha
  extraconstr = list(
    A=matrix(c(1,0,0,0,0), 1, 5),
    e = matrix(0,1,1)
  )
)

proj_from_DGP_mesh_to_intensity_pred_locs <- fm_evaluator(
  mesh = mexdolphin_sf$mesh, 
  # all models predict on the locations of this integration scheme, with which they are scored.
  # and the scores are integrated using the integration scheme weights
  loc = fm_int(list(geometry = mexdolphin_sf$mesh), samplers=mexdolphin_sf$ppoly)$geometry
)
# for model scoring later
true_detect_prob <- pany


############## Simulation


set.seed(800)
how_verbose = 1
source("one iteration of realisation and fitting.R")
some_results

############## Repeated simulations
set.seed(1234)
nsims <- 20
results <- NULL
how_verbose = 0

# for saving the results
file_name <- paste(format(Sys.time(), "%d-%m-%Y %H-%M"), "simulation results.rda")

for (i in 1:nsims){
  catt("Simulation", i, "of", nsims)
  
  # this script will use fit all the models and store the differenced scores in
  # a `some_results` df
  source("one iteration of realisation and fitting.R")

  results <- bind_rows(results, some_results)

  #save results early just in case
  saveRDS(results, file=file_name)
}


# plot the difference in scores, a difference of zero implies the models perform the same wrt to that score
# all scores are negatively orientated so a positive difference means the simpler merged observer model
# performed worse
# results <- readRDS("28-07-2026 15-22 simulation results.rda")

g <- ggplot(results, aes(fill=model)) +
  expand_limits(y=0) +
  geom_abline(intercept = 0, slope = 0, colour="red", linewidth=2) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
ds_ll <-g + geom_boxplot(aes(y=DS_loglambda)) + labs(title="integrated DS on log lambda")
msep <- g + geom_boxplot(aes(y=loglambda_SE)) + labs(title="integrated SE on mean log lambda")
maep <- g + geom_boxplot(aes(y=lambda_AE)) + labs(title = "integrated AE on median lambda")
maedetectp <- g + geom_boxplot(aes(y=detect_AE)) + labs(title="mean AE on detection prob")
ds_avg_p <- g + geom_boxplot(aes(y=DS_avg_prob)) + labs(title="DS on average detection prob")

(ds_ll + msep) / ( maep + ds_avg_p)

ds_ll
msep + coord_cartesian(ylim=c(NA, 1e07))
maep
ds_avg_p
