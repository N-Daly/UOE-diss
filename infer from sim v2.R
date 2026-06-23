rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19)

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)
library("reshape2")
library("patchwork")

my_dir <- r"(C:\Users\ND\OneDrive - University of Edinburgh\Dissertation)"
setwd(my_dir)
source("simulate 2D ground truth.R")

true_sigmaA <- 3; true_sigmaB <- 7
set.seed(123)

# simulate a ground truth and sample from it
# in ground_truth we have samples from the thinned point process, the integrated lambda for the 
# underlying process and the total abundance sample from the unthinned point process
ground_truth <- simulate_lcgp_distance_thinning(
  detect_func = function(distance, sigma){ exp(-0.5 * (distance / sigma)^2 ) },
  detect_func_paramA = true_sigmaA,
  detect_func_paramB = true_sigmaB,
  approx_sampling_points = 100
)


dd <- ground_truth$samples_df
overall_lambda <- ground_truth$overall_lambda
nrow(dd)

#preprocessing

dd$detected <- rep(3, length(dd$detectA))
dd$detected[ dd$detectA & !dd$detectB ] = 1
dd$detected[ !dd$detectA & dd$detectB ] = 2



########### PDFs for detection functions


hn <- function(distance, sigma){ exp(-0.5 * (distance / sigma)^2 ) }
log_hn <-  function(distance, sigma){ -0.5 * (distance / sigma)^2 }

######## the PDFs for both detection state and distance

log_g_2observer_hn <- function(detected, distance, sigmaA, sigmaB, eps=1e-6){
  
  log_terms <- numeric(length(detected))
  # A alone
  ii <- detected==1 
  log_terms[ii] = log_hn(distance[ii], sigmaA) + log1p(-hn(distance[ii], sigmaB)*(1-eps) )
  # B alone
  ii <- detected==2
  log_terms[ii] = log1p(-hn(distance[ii], sigmaA)*(1-eps) ) + log_hn(distance[ii], sigmaB)
  # A,B
  ii <- detected==3
  log_terms[ii] = log_hn(distance[ii], sigmaA) + log_hn(distance[ii], sigmaB)
  
  log_terms
}

log_g_1observer_hn <- function(detected, distance, sigmaSingle, eps=1e-6){
  log_hn(distance, sigmaSingle)
}

####### the PDFs for detection at a given distance

detect_func_2observer <- function(distance, sigmaA, sigmaB){
  #ln( (1-pA)(1-pB) )
  log_terms <- log1p( -hn(distance, sigmaA) ) + log1p( -hn(distance, sigmaB) )
  # cant think of a better way to do this bit
  #  1-terms  =  1- e^log_terms =  -1*( e^logterms - 1 ) 
  # my assumption is that if the probs are high this will be stable for the small log terms
  # and if the probs are low this will be also fine for the bigger log terms
  -expm1(log_terms)
}

detect_func_1observer <- function(distance, sigmaSingle){
  #ln( (-pA )
  log_terms <- log1p( -hn(distance, sigmaSingle) ) 
  # cant think of a better way to do this bit
  #  1-terms  =  1- e^log_terms =  -1*( e^logterms - 1 ) 
  # my assumption is that if the probs are high this will be stable for the small log terms
  # and if the probs are low this will be also fine for the bigger log terms
  -expm1(log_terms)
}


######## set up before modelling

distance_grid <- seq(0,8, by = 0.1)
expanded_domain <- list(geometry=mesh1, detected=1:3, distance=fm_mesh_1d(distance_grid) )

# for sampling the posterior detection function later
mc_samples <- 500
results_df <- data.frame(
  distance = rep(distance_grid, times=mc_samples),
  true = rep( detect_func_2observer(distance_grid, true_sigmaA, true_sigmaB), times=mc_samples )
)

#make a mesh out of the datapoints
mesh1 <- fm_mesh_2d(
  loc=dd, 
  max.edge= c (.5, 1), cutoff=.3,
  crs=st_crs(dd) # apparently inlabru wants a valid CRS
)
plot(mesh1);points(dd, col ="red")


# somewhat guessing a prior P(practical range > 6) = 0.01
# lets put a prior P(sigmaU > 1 ) = 0.05 knowing sigmaU=1 to see what happens

matern_prior <- inla.spde2.pcmatern(
  mesh1,
  prior.range = c(6, 0.01),
  prior.sigma = c(1, 0.05)
)


########### one observer model

cmp <- ~ intercept(1) +
  sigmaSingle(1,
              prec.linear = 1,
              marginal = bm_marginal(qexp, pexp, dexp, rate = 1 / 8)
  ) +
  spde(
    main=geometry, model = matern_prior
  )

form <- geometry + detected + distance ~ intercept + log_g_1observer_hn(detected, distance, sigmaSingle) + spde

fit_1observer <- lgcp(
  components = cmp,
  formula=form,
  domain = expanded_domain,
  data=dd,
  options = list(
    control.compute=list(dic=T, return.marginals.predictor=TRUE) # 
  )
)
fit_1observer
gc()

results_df$one_observer <- c(generate(
  fit_1observer,
  formula = ~ detect_func_1observer(distance_grid, sigmaSingle),
  n.samples=mc_samples 
))

########### two observer model

cmp <- ~ intercept(1) + 
  sigmaA(1,
         prec.linear = 1,
         marginal = bm_marginal(qexp, pexp, dexp, rate = 1 / 8)
  ) +
  sigmaB(1,
         prec.linear = 1,
         marginal = bm_marginal(qexp, pexp, dexp, rate = 1 / 8)
  ) + 
  spde(
    main=geometry, model = matern_prior
  )


form <- geometry + detected + distance ~ intercept + log_g_2observer_hn(detected, distance, sigmaA, sigmaB) + spde

fit_2observer <- lgcp(
  components = cmp,
  formula=form,
  domain = expanded_domain,
  data=dd,
  options = list(
    control.compute=list(dic=T, return.marginals.predictor=TRUE)
  )
)
fit_2observer
gc()

results_df$two_observers = c(generate(
  fit_2observer, 
  formula = ~ detect_func_2observer(distance_grid, sigmaA, sigmaB), 
  n.samples=mc_samples
))



######## basic model comparison and inspection of posterior dists

# the two observer model does in fact have a better DIC
deltaIC(fit_1observer, fit_2observer)

# not good enough at ggplot to know how to do this in a better way

# #this may take a minute. using "loess" as the smoother is slow but gam and lm perform poorly at distance=0,1
# results_df <- melt(results_df ,  id.vars = 'distance', variable.name = 'values')
# ggplot(results_df, aes(distance, value)) + coord_cartesian( ylim = 0:1) +
#   geom_smooth(method="loess", se=F) + # need method="loess" if using <1000 posterior samples
#   aes(colour = values) 



compare_posterior_spatial_params <- function(
  mdl1, mdl2,
  spde_params= c("log.variance", "log.range", "matern.correlation", "matern.covariance"),
  spde_comp_name = "spde"
  ){
  # for visually inspecting each model's estimates of the spatial params  
  
  #https://stackoverflow.com/questions/24309910/how-to-get-name-of-variable-in-r-substitute#24310574
  mdl1_name = deparse(substitute(mdl1))
  mdl2_name = deparse(substitute(mdl2))
  
  #plot each param one by one 
  for (spde_param in spde_params){
    
    top = spde.posterior(fit_1observer, spde_comp_name, what= spde_param)
    bottom = spde.posterior(fit_2observer, spde_comp_name, what= spde_param)
    # you need print to show a ggplot inside a for loop https://statisticsglobe.com/print-ggplot2-plot-within-for-loop-in-r
    # is this code beautiful ? no.
    print(
      ( plot(top) + ggtitle(paste(mdl1_name,"    ", spde_param)) ) / 
      ( plot(bottom) + ggtitle(mdl2_name) )
    )
  }
}

# i cant see any difference 
compare_posterior_spatial_params(fit_1observer, fit_2observer)


estimate_posterior_abundance <- function(mdl, the_mesh=mesh1, num_observed_points, exposure=1){
  # giving a large range of possible values relative to observed points due 
  # to unknown variability in the detections
  # multiplictive factor could be informed using some quantile of posterior detection function
  broad_range <- ceiling( num_observed_points * c(0.5, 2) )
  abundance <- broad_range[1]:broad_range[2]
  predict(
    mdl, fm_int(the_mesh),
    formula = ~ data.frame(
      abundance, 
      dpois(
        abundance,
        sum( weight*exposure*exp(intercept + spde) )
      )
    ), n.samples = 1000
  )
}

estimate_lambda <- function(mdl, the_mesh = mesh1, exposure =domain_expansion_factor){
  predict(
    mdl, fm_int(the_mesh),
    formula = ~ exposure * sum(weight * exp(intercept + spde)),
    n.samples = 1000
  )
}

######## integrated (underlying) lambda
overall_lambda


# I dont have a justification for the different exposures but its too coincidently 
# that one is out by a factor of 3 when there are 3 detection states and the distance 
# domain is [0,8]
fit1_overall_lambda <- estimate_lambda(fit_1observer, exposure=3*8)
fit2_overall_lambda <- estimate_lambda(fit_2observer,exposure = 8)

# 2 observer has tighter HDPIs
fit1_overall_lambda; fit2_overall_lambda

fit1_abundance = estimate_posterior_abundance(fit_1observer, mesh1, nrow(dd), exposure = 3*8)
fit2_abundance = estimate_posterior_abundance(fit_2observer, mesh1, nrow(dd), exposure = 8)

tail(fit1_abundance)$sd
fit1_overall_lambda

# reassuringly, the 2 observer model is well alligned with the true abundance
# I didn't have any expectation of the posterior for the mispecified 1 observer model

#https://stackoverflow.com/questions/19622063/adding-vertical-line-in-plot-ggplot
ggplot(NULL, aes(x=abundance, y=mean)) + 
  geom_line(data=fit2_abundance, color="green") +
  geom_line(data=fit1_abundance, color="blue") +
  geom_vline(xintercept=overall_lambda, color="red")


#this is hard to view and i dont know how to interpret the the tightening of the bounds for larger 
# abundances

# https://stackoverflow.com/questions/9109156/ggplot-combining-two-plots-from-different-data-frames
ggplot(NULL, aes(x=abundance, y=mean, ymin=q0.025, ymax=q0.975)) + 
  geom_ribbon(data=fit2_abundance, alpha = .1, fill = "green") +
  geom_line(data=fit2_abundance, color="green", alpha=1) +
  geom_ribbon(data=fit1_abundance, alpha = .1, fill = "blue") +
  geom_line(data=fit1_abundance, color="blue", alpha=1) +
  geom_vline(xintercept=overall_lambda, color="red")

ground_truth$true_abundance; nrow(dd)




