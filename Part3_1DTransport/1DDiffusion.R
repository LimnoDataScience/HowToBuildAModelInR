time = 100
space = 100
conc <- matrix(0, nrow = space, ncol = time) 
# our results in a matrix: 100 seconds times 100 m over the depth
K = 0.5 # diffusion coefficient, unit: m2/s
dx = 1 # our spatial step, unit: m
dt = 1 # our time step, unit: s
conc[, 1] = dnorm(seq(1,100,1), mean = 50, sd = 1) * 100
# initial conc. is defined vertically through a normal distribution, unit: -

for (n in 2:ncol(conc)){ # time index
  for (i in 2:(nrow(conc)-1)){ # space index
    conc[i, n] = conc[i, n-1] + K * dt / dx**2 * (conc[i+1, n-1] - 2 * conc[i, n-1] + conc[i-1, n-1]) # our FTCS scheme 
  }
}

# Let's plot this!
library(tidyverse)
library(reshape2)

time =  paste0(seq(1,ncol(conc)))
df <- data.frame(cbind(time, t(conc)) )
colnames(df) <- c("time", as.character(paste0(seq(1,nrow(conc)))))

m.df <- reshape2::melt(df, "time")
m.df$time <- time

ggplot(m.df, aes(as.numeric(time), as.numeric(variable))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(-1,100),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth') +
  labs(fill = 'Conc. [%]')

