library(ggplot2)
library(ggrepel)
library(gridExtra)

g1<- ggplot(comparison, aes(MINDIST, KLAVG2, label=Dataset))  + 
  geom_point(aes(colour=std_diff), size = 3) +
  geom_text_repel(colour="black", size=2.5) +
  scale_colour_gradientn(colours = topo.colors(4)) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0.3,1)) + ylim(c(0.3,1)) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g2<- ggplot(comparison, aes(MINDIST, KLAVG3, label=Dataset))  + 
  geom_point(aes(colour=std_diff), size = 3) +
  geom_text_repel(colour="black", size=2.5) +
  scale_colour_gradientn(colours = topo.colors(4)) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0.3,1)) + ylim(c(0.3,1)) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g3<- ggplot(comparison, aes(APXDIST, KLAVG2, label=Dataset))  + 
  geom_point(aes(colour=std_diff), size = 3) +
  geom_text_repel(colour="black", size=2.5) +
  scale_colour_gradientn(colours = topo.colors(4)) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0.3,1)) + ylim(c(0.3,1)) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g4<- ggplot(comparison, aes(APXDIST, KLAVG3, label=Dataset))  + 
  geom_point(aes(colour=std_diff), size = 3) +
  geom_text_repel(colour="black", size=2.5) +
  scale_colour_gradientn(colours = topo.colors(4)) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0.3,1)) + ylim(c(0.3,1)) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

grid.arrange(g1, g2, nrow=2)
grid.arrange(g3, g4, nrow=2)

g4
