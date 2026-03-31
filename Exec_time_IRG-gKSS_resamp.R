setwd("...") #for windows working directory
source("source-codes.R")

g.kernel= CalculateWLKernel
level = 3
alpha = 0.05

############ Single IRG-gKSS ###########

B = c(5, 25)
p_ER = 0.06

size = c(100, 200, 300, 400)
time_test = matrix(0, nrow = length(size), ncol = length(B))

for (k in 1:length(size)) {
  n = size[k]
  C = rep(1, n)
  p = matrix(p_ER, n, n)
  diag(p) = 0
  G = sample_IRG(p)
  for(j in 1:length(B)){
    s_size = ceiling(choose(n,2)*(B[j]/100))
    sample.index = sample_index(n, s.size = s_size, replace = TRUE)
    start.time <- Sys.time()
    generate.one.gKSS.sampled(G, C, p,sample.index, g.kernel,level)$stats.value
    end.time <- Sys.time()
    time_test[k,j] <- as.numeric (end.time - start.time, units = "mins")
    print(paste0("Computation time for n = ",size[k], " and B = ",B[j], " is ", time_test[k,j], "minutes")) 
  }

}

rownames(time_test) <- c(100, 200, 300, 400)
colnames(time_test) <- c("5%", "25%")

time <- as.data.frame.table(time_test)
colnames(time) = c("Size", "B", "Time")

ggplot(time, aes(x = Size, y = Time)) + 
  geom_point(aes(color = B, shape = B, group = B))+
  geom_line(aes(color = B, linetype = B, group = B)) +
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))




