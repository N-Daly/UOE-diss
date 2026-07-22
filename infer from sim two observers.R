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
source("function Definitions.R")


########## calculating the scores
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
# mesh <- readRDS("mesh20-july qual loc 2-20.rda")
# start <- proc.time()
# ips <- fm_int( mesh, samplers=sim_info$buffered_transects$geometry)
# ips$distance <- dist_to_nearest_line_transect(ips$geometry, mexdolphin_sf$samplers)
# end <- proc.time()
# print((end-start)[3])
# ips

# mesh <- mexdolphin_sf$mesh
# ips <- fm_int(list(geometry=mexdolphin_sf$mesh), samplers=st_buffer(mexdolphin_sf$samplers, 8, endCapStyle = "FLAT") )
# ips$distance <- dist_to_nearest_line_transect(ips$geometry, mexdolphin_sf$samplers)

mesh <- readRDS("mesh20-july from hex e6 20 only.rda")
ips <- readRDS("ips20-july from hex e1-5 20 only.rda")

# ips <- st_as_sf( readRDS("ips_interiorTransects_2subdivisions.rda") )

# so the one observer models have only a geometry dimension
# the two observer model needs a detection state dimension - who spotted the animal
# here i naively expand the existing ips uniformly across this dimension
# fingers crossed
ips_with_detection_states <- ips[rep(1:nrow(ips), each=3),]
ips_with_detection_states$detected <- rep(1:3, times=nrow(ips))
ips_with_detection_states$weight = ips_with_detection_states$weight/3


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


######## fitting models

# Model what A and B saw as a combined observer with a hn detection function
catt("fitting merged")
start <- proc.time()
fit_merged_observers <- model_one_observer_hn(
  observed_points, sim_info, matern_prior, ips=ips, bru_verbose=1
)
end <- proc.time()
print((end-start)[3])

start <- proc.time()
list_of_models[[1]] <- get_preds_from_one_observer_hn_fit(fit_merged_observers, ips)
end <- proc.time()
print((end-start)[3])


# Model what A and B saw as a two observer likelihood with hn detection functions
catt("fitting two")
start <- proc.time()
fit_two_observers <- model_two_observers_hn(
  observed_points, sim_info, matern_prior, ips = ips_with_detection_states, bru_verbose=1
)
end <- proc.time()
print((end-start)[3])

start <- proc.time()
list_of_models[[2]] <- get_preds_from_two_observers_hn_fit(fit_two_observers, ips)
end <- proc.time()
print((end-start)[3])


# Model what A and B saw as a combined observer with a spline(like) detection function
catt("fitting spline")
start <- proc.time()
fit_one_spline <- model_one_observer_spline(
  observed_points, sim_info, matern_prior, detect_matern, ips=ips, bru_verbose=1
)
end <- proc.time()
print((end-start)[3])

start <- proc.time()
list_of_models[[3]] <- get_preds_from_one_observer_spline_fit(fit_one_spline, ips)
end <- proc.time()
print((end-start)[3])


initial_results <- get_scoring_differences(
  list_of_models,
  ips,
  true_detect_prob,
  true_loglambda_at_ip
)
initial_results


####### sanity checking
dp <- list_of_models[[3]]$pred_detect
plot(
  dists,
  true_detect_prob
)
lines(dists, dp$mean, col="red")


s <- predict(
  fit_one_spline,
  data.frame(distance=dists),
  formula = ~ spline_spde,
  n.samples = 500
)
head(s)
plot(dists, s$mean)
###### end sanity checking

# ############## Repeated simulations
# set.seed(1234)
# nsims <- 20
# results <- NULL
# # for saving the results
# file_name <- paste(format(Sys.time(), "%d-%m-%Y %H-%M"), "simulation results.rda")
# 
# for (i in 1:nsims){
#   catt("Simulation", i, "of", nsims)
#   
#   # need to save space and these model fits are large
#   rm(list_of_models, fit_merged_observers, fit_two_observers, fit_one_spline)
#   invisible(gc())
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
#   #the sampled points are in observed_points
#   observed_points <- sim_info$samples_df
#   list_of_models = list()
# 
#   # Model what A and B saw as a combined observer with a hn detection function
#   catt("fitting merged")
#   start <- proc.time()
#   fit_merged_observers <- model_one_observer_hn(observed_points, sim_info, matern_prior, ips=ips)
#   end <- proc.time()
#   print((end-start)[3])
#   
#   start <- proc.time()
#   list_of_models[[1]] <- get_preds_from_one_observer_hn_fit(fit_merged_observers, ips)
#   end <- proc.time()
#   print((end-start)[3])
#   
#   
#   # Model what A and B saw as a two observer likelihood with hn detection functions
#   catt("fitting two")
#   start <- proc.time()
#   fit_two_observers <- model_two_observers_hn(
#     observed_points, sim_info, matern_prior, ips = ips_with_detection_states
#   )
#   end <- proc.time()
#   print((end-start)[3])
#   
#   start <- proc.time()
#   list_of_models[[2]] <- get_preds_from_two_observers_hn_fit(fit_two_observers, ips)
#   end <- proc.time()
#   print((end-start)[3])
#   
#   
#   # Model what A and B saw as a combined observer with a spline(like) detection function
#   catt("fitting spline")
#   start <- proc.time()
#   fit_one_spline <- model_one_observer_spline(
#     observed_points, sim_info, matern_prior, detect_matern, ips=ips
#   )
#   end <- proc.time()
#   print((end-start)[3])
#   
#   start <- proc.time()
#   list_of_models[[3]] <- get_preds_from_one_observer_spline_fit(fit_one_spline, ips)
#   end <- proc.time()
#   print((end-start)[3])
#   
#   
#   some_results <- get_scoring_differences(
#     list_of_models,
#     ips,
#     true_detect_prob,
#     true_loglambda_at_ip
#   )
#   
#   #save results early just in case
#   results <- rbind(results, some_results)
#   # wrangling with classes
#   results$model <- as.character(results$model)
#   results[-5] <- apply(results[-5], 2, as.numeric)
#   
#   saveRDS(results, file=file_name)
# }


# plot the difference in scores, a difference of zero implies the models perform the same wrt to that score
# all scores are negatively orientated so a positive difference means the simpler merged observer model
# performed worse
# results <- readRDS("18-07-2026 19-13 simulation results.rda")

g <- ggplot(results, aes(fill=model)) + 
  expand_limits(y=0) +  
  geom_abline(intercept = 0, slope = 0, colour="red", linewidth=2) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
dsp <-g + geom_boxplot(aes(y=DS)) + labs(title="integrated DS on log lambda")
msep <- g + geom_boxplot(aes(y=loglambda_SE)) + labs(title="integrated SE on mean log lambda")
maep <- g + geom_boxplot(aes(y=lambda_AE)) + labs(title = "integrated AE on median lambda")
maedetectp <- g + geom_boxplot(aes(y=detect_AE)) + labs(title="mean AE on detection prob")

(dsp + msep) / ( maep + maedetectp)



