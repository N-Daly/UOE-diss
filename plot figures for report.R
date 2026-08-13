# rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19);options(digits=3);invisible(gc())

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggspatial)
library(latex2exp)


my_dir <- r"(C:\Users\ND\OneDrive - University of Edinburgh\Dissertation\UOE-diss)"
setwd(my_dir)
source("simulate 2D ground truth.R")
source("observer models.R")
source("function Definitions.R")
source("mesh construction methods.R")

# this is the one used for the model's spatial effects
if (FALSE){
mesh <- make_spatially_varying_mesh3(60, coarse=60)

tran_segm <- st_buffer(mexdolphin_sf$samplers, 8, endCapStyle = "FLAT")

(ggplot() +
  gg(
    mesh,
    edge.color = "black"
  ) +
  gg(
    tran_segm,
    alpha = .1,
    colour = "blue"
  ) +
  gg(
    mexdolphin_sf$ppoly,
    alpha = 0.2,
    color="red"
  ) +
  ggtitle("Mesh of the study region") +
  theme(
    plot.title=element_text(size=rel(2),face="bold"),
    axis.title = element_blank()
    )
)
setwd("figs")
ggsave("modelMesh.pdf", height = 10, width = 20)
setwd("..")
}
######### mesh construction plots

if(FALSE){
base <- ggplot() +
  gg(mexdolphin_sf$ppoly, col = "red", alpha= .2) +
  theme(
    plot.title=element_text(size=rel(2),face="bold"),
    axis.title = element_blank()
  )

# demonstrate mesh construction method 1 - spatially varying edge lengths
mm <- make_spatially_varying_mesh(10)

base + gg(mm) +
  ggtitle("A mesh with spatially varying edge lengths")

setwd("figs")
ggsave("meshConstructionMethod1.pdf", height = 15 ,width = 20)
setwd("..")

# demonstrate mesh construction method 2 - line transect contained within edges
nonoverlapping <- get_nonoverlapping_samplers()
mm <- make_spatially_varying_mesh2(param=0, nonoverlapping_line_transects = nonoverlapping)

base + gg(mm) +
  ggtitle("A mesh composed of line transects and a hexagonal lattice")

setwd("figs")
ggsave("meshConstructionMethod2.pdf", height = 15 ,width = 20)
setwd("..")


# demonstrate mesh construct method 3 -  the two hexagonal lattices
mm <- make_spatially_varying_mesh3(5, coarse = 20)

ggplot() + gg(mm) +
  ggtitle("A mesh composed of fine and coarse hexagonal lattices") +
  gg(mexdolphin_sf$ppoly, col = "red", alpha= .2) +
  theme(
    plot.title=element_text(size=rel(2),face="bold"),
    axis.title = element_blank()
  )

setwd("figs")
ggsave("meshConstructionMethod3.pdf", height = 15 ,width = 20)
setwd("..")
}

############# simulation results plots

setwd("sim_results")
# r <- readRDS("09-08-2026 00-35 30simulation results.rda")
r <- readRDS("13-08-2026 11-18 spatrange250 simulation results.rda")
setwd("..")
r

old_models <- c("two_obs_hn", "one_obs_HR", "one_obs_spline" ,  "two_obs_HR")
new_models <- c("Dual HN", "HR", "Spline", "Dual HR")
rename  <- function(name, old_names, new_names){ 
  i <- match(name, old_names)
  new_names[i] 
}

r$model <- rename(r$model, c(old_models, "spline"), c(new_models, "Spline"))

g <- ggplot(
  r, 
  aes(
    fill=model
    )
  ) +
  expand_limits(y=0) +
  geom_abline(intercept = 0, slope = 0, colour="red", linewidth=2) +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.text = element_text(size=rel(2)),
    legend.title = element_text(size=rel(2)),
    axis.title.y  = element_text(size=rel(2)),
    plot.title = element_text(size=rel(2))
  ) 
    
ds_ll <- g + geom_boxplot(aes(y=DS_loglambda)) +
  labs(
    title="Integrated Dawid-Sebastiani score on log intensity",
    y = "Differenced Dawid-Sebastiani score"
  ) 

se_ll <- g + geom_boxplot(aes(y=loglambda_SE)) + 
  labs(
    title="Integrated Squared Error on mean log intensity",
    y = "Differenced Squared Error score"
  )

ae_l <- g + geom_boxplot(aes(y=lambda_AE)) + 
  labs(
    title = "Integrated Absolute Error on median intensity",
    y = "Differenced Absolute Error score"  
  )

ds_avg_prob <- g + geom_boxplot(aes(y=DS_avg_detect_prob)) +
  labs(
    title="Dawid Sebastiani score on average detection probability",
    y = "Differenced Dawid-Sebastiani score"
  )

ds_detect_prob <- g + geom_boxplot(aes(y=DS_detect_prob)) + 
  labs(
    title="Integrated Dawid Sebastiani score on detection probability",
    y = "Differenced Dawid-Sebastiani score"
  )

ae_detect_prob <- g + geom_boxplot(aes(y=detect_AE)) +
  labs(
    title = "Integrated Absolute error on median detection probability",
    y = "Differenced Absolute error"
  )

h <- 9; w<- 9
setwd("figs")
ds_ll #  + coord_cartesian(ylim = c(NA, 5e05))
# ggsave("results_DS_loglambda.pdf", height = h, width = w)

se_ll
# ggsave("results_SE_loglambda.pdf", height = h, width = w)

ae_l #+ coord_cartesian(ylim = c(-2500, NA))
# ggsave("results_AE_lambda.pdf", height = h, width = w)

# ds_avg_prob
# ggsave("results_DS_avg_prob.pdf", height = h, width = w)

# ds_detect_prob  + coord_cartesian(ylim = c(-3e05, 5e05))
# ggsave("results_DS_detect_prob.pdf", height = h, width = w)

ae_detect_prob
# ggsave("results_AE_detect_prob.pdf", height = h, width = w)

setwd("..")


no <- scale_fill_discrete(guide="none")

((ds_ll + no)+(se_ll+no))/((ae_l+no)+(ds_detect_prob+no))

########## artifacts from plotting dual HR sensitivity analysis
# setwd("sim_results")
# r <- readRDS("12-08-2026 17-30 dual HR sens ans.rda.rda")
# setwd("..")

# old_hr_priors <- c(
#   "Gamma Gamma(2,1) Sigma Exp(1/8)",
#   "Gamma Gamma(2,1) Sigma Exp(1/2)" ,         
#   "Gamma Unif(0.01, 10) Sigma Exp(1/2)" ,
#   "Gamma Unif(0.01, 10) Sigma Exp(1/8)" ,     
#   "Gamma Unif(0.01, 10) Sigma Unif(0.001, 15)",
#   "Gamma Gamma(2,1) Sigma Gamma(4, 1)"
# )
# 
# new_hr_priors <- c(
#   "Gamma(2,1), Exp(1/8)" ,
#   "Gamma(2,1), Exp(1/2)" ,         
#   "Unif(0.01, 10), Exp(1/2)"   ,
#   "Unif(0.01, 10), Exp(1/8)"  ,     
#   "Unif(0.01, 10), Unif(0.001, 15)",
#   "Gamma(2,1), Gamma(4, 1)"
# )
# 
# r$model <- rename(r$model, old_hr_priors, new_hr_priors)
# 
# guides(fill = guide_legend(title = TeX("Choice of priors for $\\gamma$ and $\\sigma$"))) 



