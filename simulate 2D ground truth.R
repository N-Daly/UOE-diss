rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19)
set.seed(123)
library("MASS")
library(INLA)
library(inlabru)
library(sf)
library(fmesher)

mm <- fm_mesh_2d(
  loc.domain=mexdolphin_sf$ppoly,
  crs=mexdolphin_sf$mesh$crs,
  max.edge= c(3,5),
  cutoff=1
)
mm

plot(mm)


# Matern correlation
cMatern <- function(dist_ij, nu, kappa) {
  # https://becarioprecario.bitbucket.io/spde-gitbook/ch-intro.html#sec:matern
  ifelse(
    dist_ij > 0,
    besselK(dist_ij * kappa, nu) * (dist_ij * kappa)^nu /(gamma(nu) * 2^(nu - 1)),
    1)
}

simulate_lcgp <- function(
    domain_boundary = inlabru::mexdolphin_sf$ppoly,
    true_nu=1.5, true_kappa = 0.8, true_sigmaU=1,
    true_beta0=0.5, approx_sampling_points = 200
){
  #get the points at which we want to sample the grf
  # grabbing them from a mesh of the gorillas dataset because why not
  
  mesh1 <- fm_mesh_2d_inla(
    loc.domain= domain_boundary,
    # max.edge = c(1), 
    # cutoff = .3,
    # max.n = approx_sampling_points,
    crs = gorillas_sf$mesh$crs
  )

  # dist matrix for these points - necessary for grf's matern covariance
  dmat <- as.matrix(dist(mesh1$loc))
  
  # sample from the grf and get the underlying log intensity
  matern_cov <- true_sigmaU * cMatern(dmat, true_nu, true_kappa)
  # ln lambda = beta0 + grf
  log_lambda_true_underlying <- true_beta0 + mvrnorm(mu=rep(0, nrow(dmat)),Sigma=matern_cov)
  
  # sample the count process given the intensity
  samples_df <- sample.lgcp(
    mesh1,
    log_lambda_true_underlying
  )
  
  #inventing the distance covariate
  # this doesnt have the same interpretation outside of distance sampling
  distance_resolution <- seq(0,8, length.out=1000)
  samples_df$distance <- sample(distance_resolution, size=nrow(samples_df), replace=T)
  
  # calculate the true overall_lambda of the entire space
  ips <- fm_int(mesh1)
  overall_lambda = sum(ips$weight * exp(log_lambda_true_underlying) )
  
  list(
    overall_lambda = overall_lambda,
    samples_df = samples_df,
    true_abundance = nrow(samples_df),
    the_mesh = mesh1
  )
  
}


simulate_lcgp_constant_thinning <- function(
    pA, pB, ...
    ){

  output_list <- simulate_lcgp(
    true_nu, true_kappa, true_sigmaU,
    true_beta0, approx_sampling_points
  )
  samples_df <- output_list$samples_df
  

  # thin the process independent of covariates
  samples_df$detectA <- runif(nrow(samples_df)) <= pA
  samples_df$detectB <- runif(nrow(samples_df)) <= pB

  ii <- ( samples_df$detectA | samples_df$detectB ) 
  
  output_list$samples_df <- samples_df[ii,]
  
  output_list

}

simulate_lcgp_distance_thinning <- function(
    detect_func, detect_func_paramA, detect_func_paramB,
    ...
  ){
  # the thinking behind the generic detect_func is that it can be extended as long as the 
  # function takes distance as the first argument.
  
  output_list <- simulate_lcgp(
    true_nu, true_kappa, true_sigmaU,
    true_beta0, approx_sampling_points
  )
  samples_df <- output_list$samples_df
  
  
  # thin the process with probabilities given by the detection function
  pA <- detect_func(samples_df$distance, detect_func_paramA)
  pB <- detect_func(samples_df$distance, detect_func_paramB)
  
  samples_df$detectA <- runif(nrow(samples_df)) <= pA
  samples_df$detectB <- runif(nrow(samples_df)) <= pB
  
  ii <- ( samples_df$detectA | samples_df$detectB ) 
  
  output_list$samples_df <- samples_df[ii,]
  output_list
}

plot(simulate_lcgp()$the_mesh)

# simulate_lcgp_constant_thinning(.5,.7)
# 
# simulate_lcgp_distance_thinning(
#   detect_func = function(distance, sigma){ exp(-0.5 * (distance / sigma)^2 ) },
#   detect_func_paramA = 3, detect_func_paramB = 7
# )
# 


