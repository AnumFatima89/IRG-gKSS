setwd("...") #for windows working directory
source("source-codes.R")
source("source-code_graphlet_py.R")

############ Single IRG-gKSS ###########
##### kernel to use

g.kernel= CalculateWLKernel
g.kernel_gr= grakel$GraphletSampling(k=3L)

#### Model parameters ####

p_ER = 0.06

size = c(25, 50, 75, 100)
time_test = matrix(0, nrow = length(size), ncol = 6)

for (k in 1:length(size)) {
  n = size[k]
  C = rep(1, n)
  p = matrix(p_ER, n, n)
  diag(p) = 0
  G = sample_IRG(p)
  
  start.time <- Sys.time()
  generate.one.gKSS(G, C, p, g.kernel,level=2)$stats.value
  end.time <- Sys.time()
  time_test[k,1] <- as.numeric (end.time - start.time, units = "mins")
  
  rm(start.time, end.time)
  
  start.time <- Sys.time()
  generate.one.gKSS(G, C, p, g.kernel,level=3)$stats.value
  end.time <- Sys.time()
  time_test[k,2] <- as.numeric (end.time - start.time, units = "mins")
  
  rm(start.time, end.time)
  
  start.time <- Sys.time()
  generate.one.gKSS(G, C, p, g.kernel,level=5)$stats.value
  end.time <- Sys.time()
  time_test[k,3] <- as.numeric (end.time - start.time, units = "mins")
  
  rm(start.time, end.time)
  
  start.time <- Sys.time()
  generate.one.gKSS(G, C, p, g.kernel,level=7)$stats.value
  end.time <- Sys.time()
  time_test[k,4] <- as.numeric (end.time - start.time, units = "mins")
  
  rm(start.time, end.time)
  
  start.time <- Sys.time()
  generate.one.gKSS(G, C, p, g.kernel = CalculateVertexEdgeHistGaussKernel, level=0.1)$stats.value
  end.time <- Sys.time()
  time_test[k,5] <- as.numeric (end.time - start.time, units = "mins")
  
  rm(start.time, end.time)
  
  start.time <- Sys.time()
  generate.one.gKSS_grakel(G, C, p, g.kernel_gr)$stats.value
  end.time <- Sys.time()
  time_test[k,6] <- as.numeric (end.time - start.time, units = "mins")
  
  rm(start.time, end.time)
  
}

rownames(time_test) <- c(25, 50, 75, 100)
colnames(time_test) <- c("WL2", "WL3", "WL5", "WL7", "VEH", "G3")

time <- as.data.frame.table(time_test)
colnames(time) = c("Size", "Graph_Kernel", "Time")

ggplot(time, aes(x = Size, y = Time)) + 
  geom_point(aes(color = Graph_Kernel, shape = Graph_Kernel, group = Graph_Kernel))+
  geom_line(aes(color = Graph_Kernel, linetype = Graph_Kernel, group = Graph_Kernel)) +
  ylim(0,15)+
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))

