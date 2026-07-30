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
# mesh <- readRDS("mesh20-july qual loc 2-20.rda")
# start <- proc.time()
# ips <- fm_int( mesh, samplers=sim_info$buffered_transects$geometry)
# ips$distance <- dist_to_nearest_line_transect(ips$geometry, mexdolphin_sf$samplers)
# end <- proc.time()
# print((end-start)[3])
# ips

mesh <- mexdolphin_sf$mesh
mesh_sub <- fm_subdivide(mesh, 4)
ips <- fm_int(list(geometry=mesh_sub), samplers=st_buffer(mexdolphin_sf$samplers, 8, endCapStyle = "FLAT") )
ips_with_detection_states <- fm_int(
  list(geometry=mesh_sub, detected = 1:3),
  samplers=st_buffer(mexdolphin_sf$samplers, 8, endCapStyle = "FLAT") )

ips_with_detection_states$distance <- dist_to_nearest_line_transect(ips_with_detection_states$geometry, mexdolphin_sf$samplers)
ips$distance <- dist_to_nearest_line_transect(ips$geometry, mexdolphin_sf$samplers)
rm(mesh_sub)
# 
# mesh <- readRDS("mesh20-july from hex e6 20 only.rda")
# ips <- readRDS("ips20-july from hex e1-5 20 only.rda")
# ips_with_detection_states <- readRDS("ipsdetected20-july from hex e1-5 20 only.rda")

# ips <- st_as_sf( readRDS("ips_interiorTransects_2subdivisions.rda") )

# so the one observer models have only a geometry dimension
# the two observer model needs a detection state dimension - who spotted the animal
# here i naively expand the existing ips uniformly across this dimension
# fingers crossed


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

dists <- seq(0,8, length.out=1000)


# for model scoring later
true_detect_prob = detect_func_2observer_sigma(dists, true_sigmaA, true_sigmaB)
true_loglambda_at_ip <- c( fm_evaluate(mexdolphin_sf$mesh, field = sim_info$log_lambda, loc = ips$geometry) )
list_of_models <- list()
how_verbose = 1


######## fitting models

# # # Model what A and B saw as a combined observer with a hn detection function
# catt("fitting merged HN")
# start <- proc.time()
# fit <- model_one_observer_hn(
#   observed_points, matern_prior, ips=ips, bru_verbose=how_verbose
# )
# end <- proc.time()
# print((end-start)[3])
# 
# start <- proc.time()
# list_of_models$one_obs_hn <- get_preds_from_one_observer_hn_fit(fit, ips)
# list_of_models$one_obs_hn$name <- "one_obs_hn"
# end <- proc.time()
# print((end-start)[3])
# 
# 
# rm(fit); invisible(gc())
# # # Model what A and B saw as a two observer likelihood with hn detection functions
# catt("fitting two HN ")
# start <- proc.time()
# fit <- model_two_observers_hn(
#   observed_points, matern_prior, ips = ips_with_detection_states, bru_verbose=how_verbose
# )
# end <- proc.time()
# print((end-start)[3])
# 
# start <- proc.time()
# list_of_models$two_obs_hn <- get_preds_from_two_observers_hn_fit(fit, ips)
# list_of_models$two_obs_hn$name <- "two_obs_hn"
# end <- proc.time()
# print((end-start)[3])
# 
# 
# # Model what A and B saw as a combined observer with a  hazard-rate detection function
# catt("fitting merged HR")
# start <- proc.time()
# fit <- model_one_observer_hr(
#   observed_points, mtrn_prior =  matern_prior, ips=ips,
#   prior_on_gamma = bm_marginal(qunif, punif, dunif, min=0.0001, max = 10),
#   bru_initial_params = list(
#     sigma = qnorm(pexp(2, rate= 1/8)),
#     gamma = qnorm(punif(2, min=0.0001, max = 10))
#   ),
#   bru_verbose = how_verbose
# )
# end <- proc.time()
# print((end-start)[3])
# 
# start <- proc.time()
# list_of_models$one_obs_HR <- get_preds_from_one_observer_hr_fit(
#   fit, ips
# )
# list_of_models$one_obs_HR$name <- "one_obs_HR"
# end <- proc.time()
# print((end-start)[3])
# 
# 
# rm(fit); invisible(gc())
# 

# # Model what A and B saw as a two observer likelihood with hazard-rate detection functions
# # A unif(.0001, 10) prior on gammaA/B
# catt("fitting two observer HR")
# start <- proc.time()
# fit <- model_two_observers_hr(
#   observed_points,
#   ips = ips_with_detection_states,
#   bru_verbose = how_verbose,
#   mtrn_prior = matern_prior,
#   prior_on_gamma = bm_marginal(qunif, punif, dunif, min=0.0001, max = 10),
#   bru_initial_params = list(
#     gammaA = qnorm(punif(2, min=0.0001, max = 10)),
#     gammaB = qnorm(punif(2, min=0.0001, max = 10))
#   )
# )
# end <- proc.time()
# print((end-start)[3])
# 
# start <- proc.time()
# list_of_models$two_obs_HR <- get_preds_from_two_observers_hn_fit(fit, ips)
# list_of_models$two_obs_HR$name <- "two_obs_HR"
# end <- proc.time()
# print((end-start)[3])
# 
# 
# 
# rm(fit); invisible(gc())
# # Model what A and B saw as a combined observer with a spline(like) detection function
# catt("fitting spline")
# start <- proc.time()
# fit <- model_one_observer_spline(
#   observed_points, matern_prior, detect_matern, ips=ips, bru_verbose=how_verbose
# )
# end <- proc.time()
# print((end-start)[3])
# 
# start <- proc.time()
# list_of_models$one_obs_spline <- get_preds_from_one_observer_spline_fit(fit, ips)
# list_of_models$one_obs_spline$name <- "one_obs_spline"
# end <- proc.time()
# print((end-start)[3])
# 
# rm(fit); invisible(gc())

# initial_results <- get_scoring_differences(
#   list_of_models,
#   ips,
#   true_detect_prob,
#   true_loglambda_at_ip
# )
# initial_results



############## Repeated simulations
set.seed(1234)
nsims <- 20
results <- NULL
how_verbose = 0

# for saving the results
file_name <- paste(format(Sys.time(), "%d-%m-%Y %H-%M"), "simulation results.rda")

for (i in 1:nsims){
  catt("Simulation", i, "of", nsims)

  # need to save space and these model fits are large
  rm(fit); invisible(gc())

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

  
  # # Model what A and B saw as a combined observer with a hn detection function
  catt("fitting merged HN")
  start <- proc.time()
  fit <- model_one_observer_hn(
    observed_points, matern_prior, ips=ips, bru_verbose=how_verbose
  )
  end <- proc.time()
  print((end-start)[3])

  start <- proc.time()
  list_of_models$one_obs_hn <- get_preds_from_one_observer_hn_fit(fit, ips)
  list_of_models$one_obs_hn$name <- "one_obs_hn"
  end <- proc.time()
  print((end-start)[3])


  rm(fit); invisible(gc())
  # # Model what A and B saw as a two observer likelihood with hn detection functions
  catt("fitting two HN ")
  start <- proc.time()
  fit <- model_two_observers_hn(
    observed_points, matern_prior, ips = ips_with_detection_states, bru_verbose=how_verbose
  )
  end <- proc.time()
  print((end-start)[3])

  start <- proc.time()
  list_of_models$two_obs_hn <- get_preds_from_two_observers_hn_fit(fit, ips)
  list_of_models$two_obs_hn$name <- "two_obs_hn"
  end <- proc.time()
  print((end-start)[3])


  # Model what A and B saw as a combined observer with a  hazard-rate detection function
  catt("fitting merged HR")
  start <- proc.time()
  fit <- model_one_observer_hr(
    observed_points, mtrn_prior =  matern_prior, ips=ips,
    prior_on_gamma = bm_marginal(qunif, punif, dunif, min=0.0001, max = 10),
    bru_initial_params = list(
      sigma = qnorm(pexp(2, rate= 1/8)),
      gamma = qnorm(punif(2, min=0.0001, max = 10))
    ),
    bru_verbose = how_verbose
  )
  end <- proc.time()
  print((end-start)[3])

  start <- proc.time()
  list_of_models$one_obs_HR <- get_preds_from_one_observer_hr_fit(
    fit, ips
  )
  list_of_models$one_obs_HR$name <- "one_obs_HR"
  end <- proc.time()
  print((end-start)[3])


  rm(fit); invisible(gc())


  # Model what A and B saw as a two observer likelihood with hazard-rate detection functions
  # A unif(.0001, 10) prior on gammaA/B
  catt("fitting two observer HR")
  start <- proc.time()
  fit <- model_two_observers_hr(
    observed_points,
    ips = ips_with_detection_states,
    bru_verbose = how_verbose,
    mtrn_prior = matern_prior,
    prior_on_gamma = bm_marginal(qunif, punif, dunif, min=0.0001, max = 10),
    bru_initial_params = list(
      gammaA = qnorm(punif(2, min=0.0001, max = 10)),
      gammaB = qnorm(punif(2, min=0.0001, max = 10))
    )
  )
  end <- proc.time()
  print((end-start)[3])

  start <- proc.time()
  list_of_models$two_obs_HR <- get_preds_from_two_observers_hn_fit(fit, ips)
  list_of_models$two_obs_HR$name <- "two_obs_HR"
  end <- proc.time()
  print((end-start)[3])



  rm(fit); invisible(gc())
  # Model what A and B saw as a combined observer with a spline(like) detection function
  catt("fitting spline")
  start <- proc.time()
  fit <- model_one_observer_spline(
    observed_points, matern_prior, detect_matern, ips=ips, bru_verbose=how_verbose
  )
  end <- proc.time()
  print((end-start)[3])

  start <- proc.time()
  list_of_models$one_obs_spline <- get_preds_from_one_observer_spline_fit(fit, ips)
  list_of_models$one_obs_spline$name <- "one_obs_spline"
  end <- proc.time()
  print((end-start)[3])

  rm(fit); invisible(gc())

  some_results <- get_scoring_differences(
    list_of_models, ips, true_detect_prob, true_loglambda_at_ip
  )

  results <- bind_rows(results, some_results)

  #save results early just in case
  saveRDS(results, file=file_name)
}


# plot the difference in scores, a difference of zero implies the models perform the same wrt to that score
# all scores are negatively orientated so a positive difference means the simpler merged observer model
# performed worse
# results <- readRDS("18-07-2026 19-13 simulation results.rda")

g <- ggplot(results, aes(fill=model)) +
  expand_limits(y=0) +
  geom_abline(intercept = 0, slope = 0, colour="red", linewidth=2) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
ds_ll <-g + geom_boxplot(aes(y=DS_loglambda)) + labs(title="integrated DS on log lambda")
msep <- g + geom_boxplot(aes(y=loglambda_SE)) + labs(title="integrated SE on mean log lambda")
maep <- g + geom_boxplot(aes(y=lambda_AE)) + labs(title = "integrated AE on median lambda")
maedetectp <- g + geom_boxplot(aes(y=detect_AE)) + labs(title="mean AE on detection prob")
ds_avg_p <- g + geom_boxplot(aes(y=detect_AE)) + labs(title="DS on average detection prob")

(ds_ll + msep) / ( maep + ds_avg_p)


