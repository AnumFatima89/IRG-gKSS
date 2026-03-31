setwd("...s") #for windows working directory
source("source-codes.R")

n=50
par(mfrow = c(2,2), mar=c(0,0,2,0)+.1)

m.links = 1
power = 1

G = sample_pa(n, power, m=m.links, directed=FALSE)
plot(G, vertex.label = "", main = paste("m = ", m.links ,"and power = ", power))


m.links = 1
power = 2

G = sample_pa(n, power, m=m.links, directed=FALSE)
plot(G, vertex.label = "", main = paste("m = ", m.links ,"and power = ", power))

m.links = 2
power = 1

G = sample_pa(n, power, m=m.links, directed=FALSE)
plot(G, vertex.label = "", main = paste("m = ", m.links ,"and power = ", power))

m.links = 2
power = 2

G = sample_pa(n, power, m=m.links, directed=FALSE)
plot(G, vertex.label = "", main = paste("m = ", m.links ,"and power = ", power))