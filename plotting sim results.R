library(ggplot2)
library("patchwork")

# the simulation results were saved as "results"
load("Sunday-sims.rda")
res <- results 
names(res)

# look at the calibration of the 95% credible intervals for spatial range
df <- as.data.frame(
  (2/nrow(res)) * xtabs(range_correct ~ model, res) 
)
ggplot(df, aes(x=model, y=Freq, fill=model)) + 
  geom_col() + 
  ylim(0,1) + labs(title="Calibration for spatial range") + 
  geom_abline(intercept=.95, slope=0, linewidth=2)


# look at the calibration of the 95% credible intervals for spatial std dev
df <- as.data.frame(
  (2/nrow(res)) * xtabs(stdv_correct ~ model, res) 
)
ggplot(df, aes(x=model, y=Freq, fill=model)) + 
  geom_col() + 
  ylim(0,1) + labs(title="Calibration for spatial stdev") + 
  geom_abline(intercept=.95, slope=0, linewidth=2)

# look at the 95% credible interval widths for range and std dev
ggplot(res[res$range_correct==1,], aes(model, range_CI_width, fill=model)) + geom_boxplot()

ggplot(res[res$stdv_correct==1,], aes(model, stdv_CI_width, fill=model)) + geom_boxplot()


# look at the MAE of both models' posterior mean for the probability of detection,
# as compared to the true probability of detection, the mean being taken over the distance domain.
ggplot(res, aes(model, MAE_detection, fill=model)) + 
  geom_boxplot() +
  ylim(0,NA) + labs(title="MAE from true probability of detection")
