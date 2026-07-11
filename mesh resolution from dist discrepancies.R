rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19);options(digits=3)

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)
library("patchwork")
library(ggspatial)

my_dir <- r"(C:\Users\ND\OneDrive - University of Edinburgh\Dissertation\UOE-diss)"
setwd(my_dir)
source("simulate 2D ground truth.R")
source("observer models.R")


### get some data
hn <- function(distance, sigma){  exp(-0.5 * (distance/sigma)^2 )  }
log_hn <-  function(distance, sigma){ -0.5 * (distance/sigma)^2 }

set.seed(123)
sim_info <- simulate_lcgp_distance_thinning(
  detect_func = hn,
  detect_func_paramA = 4,
  detect_func_paramB = 0,
  true_beta0 = -5,
  true_rho = 500, true_sigma_GRF=1
)

pts <- sim_info$samples_df
actual_distance <- pts$distance

# reconstruct the points' perpendicular distance using the mesh 
eval_mesh <- function(a_mesh){

  field <- dist_to_nearest_line_transect(fm_vertices(a_mesh), sim_info$line_transects)

   approximated_distance <- fm_evaluate(
    a_mesh, field = field, loc= pts$geometry
  )
  
  approximated_distance
}

# plot the reconstructed distances against actual ones
# and get the mean absolute deviation, as a percentage of the true value
plot_results <- function(reconstructed, actual = actual_distance){

  mean_percent_error <-  mean( abs(actual-reconstructed)/actual )
  
  plot(
    actual, reconstructed,
    xlim = range(0, actual),
    ylim = range(0, actual, reconstructed),
    main = paste(
      "mean abs % error", round(mean_percent_error, 2)
    )
  )
  abline(0:1, col ="red", lwd=2)
  
  mean_percent_error
}

# fit some model to see how long it takes with a given ips from a mesh
fit_model_with_it <- function(mesh, ips){
  
  hn <- function(distance, sigma){  exp(-0.5 * (distance/sigma)^2 )  }
  log_hn <-  function(distance, sigma){ -0.5 * (distance/sigma)^2 }
  
  matern_prior <- inla.spde2.pcmatern(
    sim_info$the_mesh, # the integration scheme is on a finer resolution than the model
    prior.range = c(600, 0.1), # true rho is 500
    prior.sigma = c(.5, 0.5)   # true sigma is 1
  )
  
  half_width <- 8
  
  cmp <- ~ Intercept(1) +
    sigma(1,
          prec.linear = 1,
          marginal = bm_marginal(qexp, pexp, dexp, rate = 1/half_width) #need to think more about this - effect on detection at end of halfwidth
    ) +
    spde(main=geometry, model = matern_prior)
  
  form <- geometry  ~ Intercept  +
    log_hn(distance, sigma) +
    log(2)+ spde
  
  
  fit <- lgcp(
    components = cmp,
    formula = form,
    data =  sim_info$samples_df,
    ips = ips,
    options = list(bru_verbose= 1,
                   # verbose = 4,
                   bru_initial = list(sigma = half_width/4) # need to review this 
    )
  )
  
  fit
}


################ the actual function ############

make_spatially_varying_mesh <- function(param){
  
  # this determines the target max edge length from each initially seeded vertex
  qual_loc <- function(locs){
    dist <- dist_to_nearest_line_transect(locs, sim_info$line_transects)
    #round the distance so that there are less unique edge lengths initially
    pmin(20, pmax(0.5, floor(dist) ))
  }
  
  # a grid of points across the study region, 40km apart
  seed_points <- fm_hexagon_lattice(sim_info$boundary, edge_len = 40)
  
  # also include points on the line transects themselves
  # so that in the following delauney triangulation the transects will be better captured
  pts_near_line <- st_cast(sim_info$line_transects$geometry, "POINT")
  seed_points <- c(seed_points, pts_near_line)
  
  # make a mesh with:
  # an extension outside the study region, around 10% farther
  # no triangle angles smaller than 23 degrees
  # there is some user set number of triangles, once reached the algorithm starts resizing edge lengths
  mesh <- fm_rcdt_2d_inla(
    loc = seed_points
    , extend = list(offset=-.1),
    refine = list(
      min.angle = 23,
      max.n = param*1000
    ),
    quality.spec = list(
      loc = qual_loc(seed_points)
    ),
    crs=sim_info$the_mesh$crs
  )
  mesh
}


#c(100, 50, 30, 25, 20) 
fit_time <- NULL
param <- NULL
error <- NULL

for (i in seq(10,60, by=5) ){
  cat("\n param =", i, "\n")
  start <- proc.time()
  cat("making mesh and ips \n")
  
  m <- make_spatially_varying_mesh(i)
  approx_dist <- eval_mesh(m)

  ips <- fm_int(
    m, samplers = sim_info$buffered_transects
  )
  ips$distance <- dist_to_nearest_line_transect(ips$geometry, sim_info$line_transects)
  
  end <- proc.time()
  print( (end-start)[3])
  
  #plot the error
  mpe <- plot_results(approx_dist)
  mtext(paste("param=",i, "# vertices =", nrow(fm_vertices(m)) ))
  
  #fit a simple model with spatial random effect and hn detection function
  start <- proc.time()

  some_model_fit <- fit_model_with_it(m, ips)

  end <- proc.time()
  print( (end-start)[3])
  
  #record stats
  fit_time <- c(fit_time, (end-start)[[3]])
  param <- c(param, i)
  error <- c(error, mpe)
  
}

#error <- c(129, 85, 77, 65, 63, 62, 56, 49, 46, 43, 40)/100
plot(
  error, 
  fit_time,
  ylim = range(0, fit_time),
  xlim = range(0, error),
  type = "c"
)
text(
  error, 
  fit_time,
  labels = param
)

ggplot() + gg(m) + gg(sim_info$boundary, color="red", alpha=.2) +
  gg(sim_info$buffered_transects, color="brown", alpha=.2)+
  ggspatial::annotation_scale(location="tr")



