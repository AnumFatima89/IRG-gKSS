setwd("...") #for windows working directory
source("source-codes.R")
library(patchwork)

g.kernel= CalculateWLKernel
level = 3
alpha = 0.05

size = c(25, 50, 75, 100)
time_test = c()

p_ER = 0.06

M= c(1, 50, 100, 150, 200)

time_test = matrix(0, nrow = length(size), ncol = length(M))

for (k in 1:length(size)) {
  n = size[k]
  C = rep(1, n)
  p = matrix(p_ER, n, n)
  diag(p) = 0
  for (j in 1:length(M)){
    for (i in 1:M[j]) {
      G = sample_IRG(p)
      start.time <- Sys.time()
      generate.one.gKSS(G, C, p, g.kernel,level)$stats.value
      end.time <- Sys.time()
      time_test[k,j] <- time_test[k,j] + as.numeric (end.time - start.time, units = "mins")
    }
    print(paste0("Computation time for n = ",size[k], " and M = ",M[j], " is ", time_test[k,j]))   
  }
}

rownames(time_test) <- c(25, 50, 75, 100)
colnames(time_test) <- c(1, 50, 100, 150, 200)

time <- as.data.frame.table(time_test)
colnames(time) = c("Size", "M", "Time")


T1 = ggplot(time, aes(x = Size, y = Time)) + 
  geom_point(aes(color = M, shape = M, group = M))+
  geom_line(aes(color = M, linetype = M, group = M)) +
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))

T2 = ggplot(time, aes(x = M, y = Time)) + 
  geom_point(aes(color = Size, shape = Size, group = Size))+
  geom_line(aes(color = Size, linetype = Size, group = Size)) +
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))

T1 + T2
