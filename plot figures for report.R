rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19);options(digits=3);invisible(gc())

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

rename  <- function(name, old_names, new_names){
  i <- match(name, old_names)
  new_names[i]
}

######### plot the mesh used for the model's spatial effects
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

# setwd("sim_results")
# r <- readRDS("09-08-2026 00-35 30simulation results.rda")
# # r <- readRDS("13-08-2026 12-33 spatrange250 new spline simulation results.rda ")
# setwd("..")
# r

# old_models <- c("two_obs_hn", "one_obs_HR", "one_obs_spline" ,  "two_obs_HR")
# new_models <- c("Dual HN", "HR", "Spline", "Dual HR")
# r$model <- rename(r$model, old_models, new_models); unique(r$model)

g <- ggplot(
  data = r, aes(fill=model)
  ) +
  geom_vline(xintercept = 0, colour="red", linewidth=2) +
  theme(
    legend.position="none",
    axis.title.x  = element_text(size=rel(2)),
    axis.title.y  = element_text(size=rel(2)),
    axis.text.y = element_text(size=rel(2)),
    plot.title = element_text(size=rel(2))
  )
    
ds_ll <- g + geom_boxplot(aes(DS_loglambda, model)) +
  labs(
    title="Integrated Dawid-Sebastiani score on log intensity",
    x = "Differenced Dawid-Sebastiani score"
  ) 

se_ll <- g + geom_boxplot(aes(loglambda_SE, model)) + 
  labs(
    title="Integrated Squared Error on mean log intensity",
    x = "Differenced Squared Error score"
  )

ae_l <- g + geom_boxplot(aes(lambda_AE, model)) + 
  labs(
    title = "Integrated Absolute Error on median intensity",
    x = "Differenced Absolute Error score"  
  )

ae_detect_prob <- g + geom_boxplot(aes(detect_AE, model)) +
  labs(
    title = "Integrated Absolute error on median detection probability",
    x = "Differenced Absolute error"
  )

# h <- 7; w<- 9
# setwd("figs")
# ds_ll
# ggsave("results_DS_loglambda.pdf", height = h, width = w)
# 
# se_ll
# ggsave("results_SE_loglambda.pdf", height = h, width = w)
# 
# ae_l
# ggsave("results_AE_lambda.pdf", height = h, width = w)
# 
# ae_detect_prob
# ggsave("results_AE_detect_prob.pdf", height = h, width = w)
# 
# setwd("..")



########## artifacts from plotting dual HR sensitivity analysis
setwd("sim_results")
r <- readRDS("12-08-2026 17-30 dual HR sens ans.rda")
setwd("..")

old_hr_priors <- c(
  "Gamma Gamma(2,1) Sigma Exp(1/8)",
  "Gamma Gamma(2,1) Sigma Exp(1/2)" ,
  "Gamma Unif(0.01, 10) Sigma Exp(1/2)" ,
  "Gamma Unif(0.01, 10) Sigma Exp(1/8)" ,
  "Gamma Unif(0.01, 10) Sigma Unif(0.001, 15)",
  "Gamma Gamma(2,1) Sigma Gamma(4, 1)"
)

new_hr_priors <- c(
  "Gamma(2,1), Exp(1/8)" ,
  "Gamma(2,1), Exp(1/2)" ,
  "Unif(0.01, 10), Exp(1/2)"   ,
  "Unif(0.01, 10), Exp(1/8)"  ,
  "Unif(0.01, 10), Unif(0.001, 15)",
  "Gamma(2,1), Gamma(4, 1)"
)

r$model <- rename(r$model, old_hr_priors, new_hr_priors); unique(r$model)

new_ylab<-  ylab(TeX("Choice of priors for $\\gamma$ and $\\sigma$"))

h <- 9; w<- 11
setwd("figs")
ds_ll + new_ylab
ggsave("sens_ans_HR_DS_loglambda.pdf", height = h, width = w)

se_ll + new_ylab
ggsave("sens_ans_HR_SE_loglambda.pdf", height = h, width = w)

ae_l + new_ylab
ggsave("sens_ans_HR_AE_lambda.pdf", height = h, width = w)

ae_detect_prob + new_ylab
ggsave("sens_ans_HR_AE_detect_prob.pdf", height = h, width = w)
setwd("..")

########## artifacts from plotting dual HN sensitivity analysis
# setwd("sim_results")
# r <- readRDS("14-08-2026 19-15 dual HN sens analysis simulation results.rda")
# setwd("..")
# 
# unique(r$model)
# 
# new_ylab<-  ylab(TeX("Choice of prior for $\\sigma$"))
# 
# h <- 7; w<- 11
# setwd("figs")
# ds_ll + new_ylab
# ggsave("sens_ans_HN_DS_loglambda.pdf", height = h, width = w)
# 
# se_ll + new_ylab
# ggsave("sens_ans_HN_SE_loglambda.pdf", height = h, width = w)
# 
# ae_l + new_ylab
# ggsave("sens_ans_HN_AE_lambda.pdf", height = h, width = w)
# 
# ae_detect_prob + new_ylab
# ggsave("sens_ans_HN_AE_detect_prob.pdf", height = h, width = w)
# 
# setwd("..")

########## 
