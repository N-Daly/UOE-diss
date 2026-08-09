rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19);options(digits=3);invisible(gc())

source("function Definitions.R")
library(reshape2)
library(ggplot2)

dd <- seq(0, 8, length.out=1000)

# pa <- hr(dd, 1.75, 1)
# pb <- hr(dd, 2, 4)
pa <- hn(dd, 2)
pb <- hn(dd, 4)
pany <- 1 - (1-pa)*(1-pb)

par(cex=1)

plot(
  dd,
  pany,
  ylim = 0:1,
  ylab = "Detection probability",
  xlab = "Distance from observers, km",
  lwd = 2,
  type="l",
  main = "True probability of detection within transect segment"
)
lines(dd, pa, col = "red", lwd = 2)
lines(dd, pb, col = "blue", lwd = 2)

legend(
  "topright",
  legend = c("Observer A and/or B", "Observer A", "Observer B"),
  col = c("black", "red", "blue"),
  lwd = "3"
)

# https://stackoverflow.com/questions/23635662/editing-legend-text-labels-in-ggplot
df <- data.frame(distance=dd, pa=pa, pb=pb, pany = pany)
dfm <- melt(df, id= "distance")

ggplot(dfm, aes(x=distance, y = value, color = variable)) +
  geom_line(lwd=2) +
  scale_color_manual(
    breaks = c("pany", "pa", "pb"),
    labels = c("Observer A and/or B", "Observer A", "Observer B"),
    values = c("black", "red", "blue")
  ) +
  ylab("Detection probability") +
  xlab("Distance from observers, km") +
  labs(title = "True probability of detection within transect segment") + 
  theme_bw() +
  guides() +
  theme(
    legend.background = element_rect(fill = "white"),
    legend.position.inside  =c(1,1)
    )
