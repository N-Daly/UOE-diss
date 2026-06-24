rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19)

library("MASS")
library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(units)


dist_to_nearest_line_transect <- function(pts, transects){
  # this should be in kilometres
  distance_to_each_transect <- st_distance(pts, transects)
  #distance to nearest line transect
  apply(distance_to_each_transect, 1, min)
  
  #  preserving the same units causes issues with exponentiation later
  # as_units(dd, units(distance_to_each_transect))
}


simulate_lcgp <- function(
    true_alpha = 2, true_rho = 500, true_sigma_GRF=1,
    true_beta0= -5,
    halfwidth_km = 8
){
  
  # this uses the mexdolphins dataset as a template
  mesh1 <- inlabru::mexdolphin_sf$mesh
  
  line_transects <- inlabru::mexdolphin_sf$samplers$geometry
  # the observable region around the line transects
  polygon_transects <- st_buffer(
    line_transects,
    set_units(halfwidth_km, "km"),
    endCapStyle = "FLAT")

  # sample from the grf and get the underlying log intensity
  grf_samples <- fmesher::fm_matern_sample(
    mesh1,
    alpha = true_alpha, rho = true_rho,
    sigma = true_sigma_GRF
  )
  cat("sampled GRF \n ")
  # ln lambda = beta0 + grf
  log_lambda_true_underlying <- true_beta0 + grf_samples

  # sample the count process, given the intensity, only in the observable regions 
  s <- Sys.time()
  samples_df <- sample.lgcp(
    mesh1,
    log_lambda_true_underlying,
    samplers = polygon_transects
  )
  e <- Sys.time()
  print(e-s)
  cat("sampled point process \n")
  
  #inventing the distance covariate
  samples_df$distance <- dist_to_nearest_line_transect(samples_df$geometry, line_transects)
  
  # calculate the true overall_lambda of the entire space
  ips <- fm_int(mesh1)
  overall_lambda = sum(ips$weight * exp(log_lambda_true_underlying) )
  
  list(
    overall_lambda = overall_lambda,
    lambda = log_lambda_true_underlying,
    samples_df = samples_df,
    true_abundance = nrow(samples_df),
    the_mesh = mesh1,
    line_transects = line_transects
  )
  
}


simulate_lcgp_constant_thinning <- function(
    pA, pB, ...
    ){

  output_list <- simulate_lcgp(...)
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
  
  output_list <- simulate_lcgp(...)
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

# 
# ss <- simulate_lcgp(true_beta0=-5)
# boxplot(ss$samples_df$distance)
# 
# plot(ss$the_mesh)
# plot(ss$samples_df$geometry, add=T)
# nrow(ss$samples_df)
# plot(mexdolphin_sf$samplers$geometry, add = T, col = "red")
# 
# 
# simulate_lcgp_constant_thinning(.5,.7)
# 
# simulate_lcgp_distance_thinning(
#   detect_func = function(distance, sigma){ exp(-0.5 * (distance / sigma)^2 ) },
#   detect_func_paramA = 3, detect_func_paramB = 7
# )



