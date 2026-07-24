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


plot_it <- function(pred_df, main = "", true=true_detect_prob, d= dists){
  plot(
    d,
    pred_df$mean,
    ylim = range(0:1, pred_df$q0.975),
    main = main,
    ylab = "Prob of detection",
    xlab = "distance"
  )
  lines(d, pred_df$q0.975)
  lines(d, pred_df$q0.025)

  lines(d, true, col = "red", lwd =2)
}

true_sigmaA <- 3; true_sigmaB <- 2; true_range <- 500; true_sigma_grf <- 1

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

dists <- seq(0,8, length.out=1000)
true_detect_prob = detect_func_2observer_sigma(dists, true_sigmaA, true_sigmaB)


mesh <- mexdolphin_sf$mesh
ips <- fm_int(
  domain=list(geometry = fm_subdivide(mesh,2), detected=1:3),
  samplers = sim_info$buffered_transects
)
ips$distance = dist_to_nearest_line_transect(ips$geometry)
matern_prior <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(600, 0.1), # true rho is 500
  prior.sigma = c(.5, 0.5)   # true sigma is 1
)
# 
# # fit a hr fixed with a "lucky" choice of gamma
# dists <- seq(0, 8, length.out=1000)
# gammas <- 1:10
# plot(
#   dists,
#   dists*0,
#   ylim = 0:1,
#   type = "n",
#   main = "hazard rate: sigma = 3, gamma varying"
# )
# for (gam in gammas){
#   lines(
#     dists, 
#     hr(dists, 3, gam),
#     col = rainbow(length(gammas))[gam]
#   )
# }
# legend(
#   x="topright",
#   legend = gammas,
#   col = rainbow(length(gammas)),
#   pch = 19
# )
# lines(dists, true_detect_prob, lwd=3)
# mtext("black is true detect prob")

fixed_gamma <- 2

# single observer with gamma fixed
fit_single <- model_one_observer_hr_fixed_gamma(
  observed_points,
  sim_info,
  ips = ips,
  mtrn_prior = matern_prior,
  bru_verbose = 1,
  fixed_gamma_val = 2
)

dp_single <- predict(
  fit_single,
  data.frame(distance=dists),
  formula = ~ hr(distance, sigma, 2)
)
plot_it(dp_single, "1 observer gamma=2")


# 2 observers but gamma still fixed
fit <- model_two_observers_hr_fixed_gamma(
  observed_points,
  sim_info,
  ips = ips,
  mtrn_prior = matern_prior,
  bru_verbose = 1,
  fixed_gamma_val = fixed_gamma
)

dp <- predict(
  fit,
  data.frame(distance=dists),
  formula = ~ detect_func_2_observer_hr(distance, sigmaA, sigmaB, fixed_gamma, fixed_gamma),
  n.samples = 1000
)
plot_it(dp, main = paste("2 observers gamma =", fixed_gamma))
bru_convergence_plot(fit)

# 
# plot(predict(
#   fit, 
#   formula =~ sigmaA,
#   n.samples=1000
# )) + 
# plot(predict(
#   fit, 
#   formula =~ sigmaB
# ))
# trialling variable gamma and sigma, with different priors for gamma

# gamma(shape=2, rate=1)
fit2 <- model_two_observers_hr(
  observed_points,
  sim_info,
  ips = ips,
  prior_on_gamma = bm_marginal(qgamma, pgamma, dgamma, shape=2, rate=1),
  bru_initial_params = list(
    gammaA = qnorm(pgamma(2, shape=2, rate=1)),
    gammaB = qnorm(pgamma(3, shape=2, rate=1))
  ),
  mtrn_prior = matern_prior,
  bru_verbose = 1
)
dp2 <- predict(
  fit2, 
  data.frame(distance=dists),
  formula = ~ detect_func_2_observer_hr(distance, sigmaA, sigmaB, gammaA, gammaB),
  n.samples = 1000
)
plot_it(dp2, main = "2 observers gamma(2,1) prior")


# exp(1)
fit3 <- model_two_observers_hr(
  observed_points,
  sim_info,
  ips = ips,
  prior_on_gamma = bm_marginal(qexp, pexp, dexp, rate=1),
  bru_initial_params = list(
    gammaA = qnorm(pexp(2, rate=1)),
    gammaB = qnorm(pexp(2, rate=1))
  ),
  mtrn_prior = matern_prior,
  bru_verbose = 1
)
dp3 <- predict(
  fit3, 
  data.frame(distance=dists),
  formula = ~ detect_func_2_observer_hr(distance, sigmaA, sigmaB, gammaA, gammaB),
  n.samples = 1000
)
plot_it(dp3, main = "2 observers exp(1) prior")
# bru_convergence_plot(fit3)

# unif(0, 10)
fit4 <- model_two_observers_hr(
  observed_points,
  sim_info,
  ips = ips,
  prior_on_gamma = bm_marginal(qunif, punif, dunif, min=0.0001, max = 10),
  bru_initial_params = list(
    gammaA = qnorm(punif(2, min=0.0001, max = 10)),
    gammaB = qnorm(punif(2, min=0.0001, max = 10))
  ),
  mtrn_prior = matern_prior,
  bru_verbose = 1
)
dp4 <- predict(
  fit4, 
  data.frame(distance=dists),
  formula = ~ detect_func_2_observer_hr(distance, sigmaA, sigmaB, gammaA, gammaB),
  n.samples = 1000
)
plot_it(dp4, main = "2 observers Unif(.0001, 10) prior")






