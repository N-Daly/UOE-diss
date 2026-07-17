rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19);options(digits=2)

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)
library(patchwork)
library(ggspatial)

my_dir <- r"(C:\Users\ND\OneDrive - University of Edinburgh\Dissertation\UOE-diss)"
setwd(my_dir)
source("simulate 2D ground truth.R")
source("observer models.R")


hn <- function(distance, sigma){  exp(-0.5 * (distance/sigma)^2 )  }
log_hn <-  function(distance, sigma){ -0.5 * (distance/sigma)^2 }

detect_func_2observer_sigma <- function(distance, sigmaA, sigmaB){
  #ln( (1-pA)(1-pB) )
  log_terms <- log1p( -hn(distance, sigmaA) ) + log1p( -hn(distance, sigmaB) )
  # cant think of a better way to do this bit
  #  1-terms  =  1- e^log_terms =  -1*( e^logterms - 1 ) 
  # my assumption is that if the probs are high this will be stable for the small log terms
  # and if the probs are low this will be also fine for the bigger log terms
  -expm1(log_terms)
}

spline_detect_func <- function(spline_effect){ exp(-spline_effect) }
ln_spline_detect_func <- function(spline_effect){ -spline_effect }


half_width <- 8
dists <- seq(0, half_width, length.out=1000)
sigmaA <- 2; sigmaB <- 0
set.seed(123)
sim_info <- simulate_lcgp_distance_thinning(
  hn, sigmaA, sigmaB
)

dd <- sim_info$samples_df

detect_mesh <- fm_mesh_1d(
  loc = seq(0,8, by=2), 
  boundary = c("dirichlet", "free"),
  degree = 2
)

ips <- st_as_sf( readRDS("ips_interiorTransects_3subdivisions.rda") )



# ips <- fm_int(
  # make_spatially_varying_mesh(60, sim_info),
  # samplers = sim_info$buffered_transects
# )
ips$distance <- dist_to_nearest_line_transect(ips$geometry, sim_info$line_transects)
# saveRDS(ips, file = "ipsFromVaryingMesh60K.rda")
head(ips)
############## model stuff

matern_prior <- inla.spde2.pcmatern(
  sim_info$the_mesh,
  prior.range = c(600, 0.1), # true rho is 500
  prior.sigma = c(.5, 0.5)   # true sigma is 1
)

# alpha=2, rho=2, sigma=.25
detect_matern <- inla.spde2.pcmatern(
  detect_mesh,
  alpha = 2,
  prior.range = c(2, 0.99), # P(rho < val) = alpha
  prior.sigma = c(0.25, 0.01) # P(sigma > val) = alpha
)

cmp <- ~ Intercept(1) +
  spline_spde(main=distance, model=detect_matern) #+
  # typical_spde(main=geometry, model=matern_prior)

form <- geometry  ~ Intercept  + #typical_spde +
  ln_spline_detect_func(spline_spde)


fit <- lgcp(
  components = cmp,
  formula = form,
  data =  dd,
  ips=ips,
  options = list(bru_verbose= 1
  )
)
fit

pred <- predict(
  fit,
  data.frame(distance=dists),
  ~ spline_detect_func(spline_spde),
  n.samples=1000
)
plot(dists, abs(pred$mean.mc_std_err/pred$mean) );abline(h=.05)


cmp2 <- ~ Intercept(1) +
  sigma(1,
        prec.linear = 1,
        marginal = bm_marginal(qexp, pexp, dexp, rate = 1/8) #need to think more about this - effect on detection at end of halfwidth
  ) +
  spde(main=geometry, model = matern_prior) 



form2 <- geometry  ~ Intercept  +
  log_hn(distance, sigma) +log(2)+ spde


fit2 <- lgcp(
  components = cmp2,
  formula = form2,
  data =  dd,
  ips=ips,
  options = list(bru_verbose= 1)
)


pred2 <- predict(
  fit2,
  data.frame(distance=dists),
  ~ hn(distance, sigma),
  n.samples = 100
)$mean

true_detect <- detect_func_2observer_sigma(dists, sigmaA, sigmaB)
plot(
  dists,
  true_detect,
  ylim = range(0:1, pred$mean),
  col = "red", type = "l"
)
lines(
  dists, pred$mean,
  lwd=2
)
lines(
  dists, pred2, 
  lwd=2, col = "blue"
)
rug(dist_to_nearest_line_transect(dd$geometry))

ggplot(NULL) + 
  geom_ribbon(data = pred, aes(x=distance, y=mean, ymin = q0.025, ymax=q0.975), fill="turquoise") +
  geom_line(data=data.frame(x=dists, y = true_detect), aes(x,y))


