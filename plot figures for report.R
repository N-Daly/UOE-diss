# rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19);options(digits=3);invisible(gc())

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggspatial)


my_dir <- r"(C:\Users\ND\OneDrive - University of Edinburgh\Dissertation\UOE-diss)"
setwd(my_dir)
source("simulate 2D ground truth.R")
source("observer models.R")
source("function Definitions.R")
source("mesh construction methods.R")

# # this is the one used for the model's spatial effects
# mesh <- make_spatially_varying_mesh3(60, coarse=60)
# 
# tran_segm <- st_buffer(mexdolphin_sf$samplers, 8, endCapStyle = "FLAT")
# 
# ggplot() +
#   gg(mesh, edge.color = "black") + 
#   gg(tran_segm, alpha = .1, colour = "blue") + 
#   gg(mexdolphin_sf$ppoly, alpha = 0.2,  color="red") +
#   ggspatial::annotation_scale(data=mexdolphin_sf$ppoly) +
#   ggtitle("Mesh of the study region") + 
#   theme(
#     plot.title=element_text(size=rel(2),face="bold"),
#     plot.margin=grid::unit(c(0,0,0,0), "mm")
#     )
# 
# ggsave("modelMesh.pdf", height = 10, width = 20)

######### mesh construction plots

# base <- ggplot() +
#   gg(mexdolphin_sf$ppoly, col = "red", alpha= .2) +
#   theme(
#     plot.title=element_text(size=rel(2),face="bold"),
#     axis.title = element_blank()
#   )
# 
# # demonstrate mesh construction method 1 - spatially varying edge lengths
# mm <- make_spatially_varying_mesh(10)
# 
# base + gg(mm) +
#   ggtitle("A mesh with spatially varying edge lengths")
# 
# setwd("figs")
# ggsave("meshConstructionMethod1.pdf", height = 15 ,width = 20)
# setwd("..")
# 
# # demonstrate mesh construction method 2 - line transect contained within edges
# nonoverlapping <- get_nonoverlapping_samplers()
# mm <- make_spatially_varying_mesh2(param=0, nonoverlapping_line_transects = nonoverlapping)
# 
# base + gg(mm) +
#   ggtitle("A mesh composed of line transects and a hexagonal lattice")
# 
# setwd("figs")
# ggsave("meshConstructionMethod2.pdf", height = 15 ,width = 20)
# setwd("..")
# 
# 
# # demonstrate mesh construct method 3 -  the two hexagonal lattices
# mm <- make_spatially_varying_mesh3(5, coarse = 20)
# 
# ggplot() + gg(mm) +
#   ggtitle("A mesh composed of fine and coarse hexagonal lattices") +
#   gg(mexdolphin_sf$ppoly, col = "red", alpha= .2) +
#   theme(
#     plot.title=element_text(size=rel(2),face="bold"),
#     axis.title = element_blank()
#   )
# 
# setwd("figs")
# ggsave("meshConstructionMethod3.pdf", height = 15 ,width = 20)
# setwd("..")


############# simulation results plots

setwd("sim_results")
r <- readRDS("12-08-2026 17-30 simulation results.rda")
# r <- readRDS("08-08-2026 compare dual HR priors simulation results.rda")
setwd("..")
r

desired_col_order <- c("two_obs_hn", "one_obs_HR", "one_obs_spline" ,  "two_obs_HR")
desired_col_names <- c("Dual HN", "HR", "Spline", "Dual HR")
rename  <- function(name){ 
  i <- match(name, desired_col_order)
  desired_col_names[i] 
}

# r$model <- rename(r$model)

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
  ) + 
  guides(fill = guide_legend(title = "Choice of priors"))



    
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

h <- 11; w<- 9
# setwd("figs")
ds_ll + coord_cartesian(ylim = c(NA, 5e05))
ggsave("sens_ans_HR_DS_loglambda.pdf", height = h, width = w)

se_ll
ggsave("sens_ans_HR_SE_loglambda.pdf", height = h, width = w)

ae_l + coord_cartesian(ylim = c(-2500, NA))
ggsave("sens_ans_HR_AE_lambda.pdf", height = h, width = w)

# ds_avg_prob
# ggsave("sens_ans_HR_DS_avg_prob.pdf", height = h, width = w)

# ds_detect_prob  + coord_cartesian(ylim = c(-3e05, 5e05))
# ggsave("sens_ans_HR_DS_detect_prob.pdf", height = h, width = w)

ae_detect_prob
ggsave("sens_ans_HR_AE_detect_prob.pdf", height = h, width = w)

# setwd("..")


no <- scale_fill_discrete(guide="none")

((ds_ll + no)+(se_ll+no))/((ae_l+no)+(ds_detect_prob+no))

theme
