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


hn <- function(distance, sigma){  exp(-0.5 * (distance/sigma)^2 )  }
log_hn <-  function(distance, sigma){ -0.5 * (distance/sigma)^2 }


half_width <- 8
dists <- seq(0, half_width, length.out=1000)
sigmaA <- 2; sigmaB <-4
set.seed(123)
sim_info <- simulate_lcgp_distance_thinning(
  hn, sigmaA, sigmaB
  )

dd <- sim_info$samples_df
mesh1 <- sim_info$the_mesh

detect_mesh <- fm_mesh_1d(
  loc = 0:8, 
  boundary = c("dirichlet", "free"),
  degree = 2
)

dists_on_region_mesh <- dist_to_nearest_line_transect(
  fm_vertices(mesh1), sim_info$line_transects
)


############## model stuff

matern_prior <- inla.spde2.pcmatern(
  mesh1,
  prior.range = c(600, 0.1), # true rho is 500
  prior.sigma = c(.5, 0.5)   # true sigma is 1
)

detect_matern <- inla.spde2.pcmatern(
  detect_mesh,
  prior.range = c(1, 0.1),
  prior.sigma = c(1, 0.1)
)

cmp <- ~ Intercept(1) + 



form <- geometry  ~ Intercept  +
  log_hn(
    dist_to_transect,
    sigma
  ) +log(2)+ spde


fit <- lgcp(
  components = cmp,
  formula = form,
  data =  dd,
  domain = list(geometry=fm_subdivide(mesh1, 3) ),
  # the observed regions are line transects expanded outwards by the transect widths, right?
  samplers = construction_info$buffered_transects,
  options = list(bru_verbose= 0,
                 # verbose = 4,
                 bru_initial = list(sigma = half_width/4) # need to review this 
  )
)

fit

summary(dd)
  
nrow(sim_info$buffered_transects)



nrow(dd)

