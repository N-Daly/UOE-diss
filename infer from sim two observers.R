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
distance_ips <- fm_int(list(distance = fm_mesh_1d(seq(0,8, length.out=250))))
dists <- distance_ips$distance

pA <- hr(dists, sigma = true_sigmaA, gamma = true_gammaA)
pB <- hr(dists, sigma = true_sigmaB, gamma = true_gammaB)
pany = detect_func_2_observer_hr(dists, true_sigmaA, true_sigmaB, true_gammaA, true_gammaB)

# for scoring later
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
  prior.range = c(500, 0.1), # true rho is 250
  prior.sigma = c(.5, 0.5)   # true sigma is 1
)

detect_mesh <- fm_mesh_1d(
  loc = c(0,2,4,6,8), 
  boundary = c("dirichlet", "free"),
  degree = 2
)

# should use
# rho <- 50; sigma <- 2; alpha = 2
detect_matern <- inla.spde2.pcmatern(
  detect_mesh,
  alpha = 2,
  prior.range = c(3, 0.9), # P(rho < val) = alpha
  prior.sigma = c(.75, 0.1), # P(sigma > val) = alpha
  extraconstr = list(
    A=matrix(c(1,0,0,0,0), 1, 5),
    e = matrix(0,1,1)
  )
)

# all models predict on the locations of this integration scheme, with which they are scored.
# and the scores are integrated using the integration scheme weights
pred_loc_ips <- fm_int(
  list(geometry = fm_subdivide(mexdolphin_sf$mesh)),
  samplers=mexdolphin_sf$ppoly
)
proj_from_DGP_mesh_to_intensity_pred_locs <- fm_evaluator(
  mesh = mexdolphin_sf$mesh, 
   loc = pred_loc_ips$geometry
)

# for model scoring later
true_detect_prob <- pany


############## Simulation


# set.seed(800)
# how_verbose = 0
# # source("one iteration of gamma fits.R")
# source("one iteration of realisation and fitting.R")
# some_results

############ Repeated simulations
nsims <- 30 
results <- NULL
how_verbose <- 0
sim_seed <- 1234 - 1

# for saving the results
file_name <- paste(format(Sys.time(), "%d-%m-%Y %H-%M"), "dual HR sens analysis simulation results.rda")

for (i in 1:nsims){
  sim_seed <- sim_seed +1
  set.seed(sim_seed)
  catt("Simulation", i, "of", nsims, "with seed", sim_seed)

  # this script will use fit all the models and store the differenced scores in
  # a `some_results` df
  # source("one iteration of gamma fits.R")
  source("one iteration of realisation and fitting.R")

  results <- bind_rows(results, some_results)

  #save results early just in case
  setwd("sim_results")
  saveRDS(results, file=file_name)
  setwd(my_dir)
  print(some_results)
}