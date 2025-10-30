######## Methodology ###########
# re-direct to the desired directory
setwd("") #for windows working directory
source("source-codes.R")

#### For spectral test
library(mgcv)
library(RMTstat)
source('utils.R')

case_id = 1

###############################
###### Global Parameters ######
###############################

alpha=0.05
n = 50
C = rep(1,n)

M = 200
m = 50

##### kernel to use

g.kernel= CalculateWLKernel
#CalculateGraphletKernel 
level = 2  # parameter for the kernel setting

{
  if(case_id ==1) p_ER = 0.02
  if(case_id ==2) p_ER = 0.06
  if(case_id ==3) p_ER = 0.10
}

p = matrix(p_ER, n, n)
diag(p) = 0

################Null Parameter set##############

p_0 = p

#######################################################
##############Simulations from null model##############
#######################################################

cores <- detectCores()/2
cl <- makeCluster(cores)
registerDoParallel(cl)
# Parallel sampling of null distribution
statistic <- foreach(i = 1:M, .combine = c, .packages = c("igraph", "graphkernels"),
                     .export = c("symmetricize")) %dopar% {
                       G = sample_IRG(p_0)
                       generate.one.gKSS(G, C, p_0, g.kernel,level)$stats.value
                     }

stopCluster(cl)
critical_value_lower = quantile(statistic,alpha/2)
critical_value_upper = quantile(statistic,(1-alpha/2))
median(statistic)

save.image(paste0("rebuttal/cliques/case_",case_id,".RData"))# path and name of output file for windows

###### Edge resampling #####

s_size = ceiling(choose(n,2)*0.05)

cores <- detectCores()/2
cl <- makeCluster(cores)
registerDoParallel(cl)

# Parallel sampling of null distribution
statistic_5 <- foreach(i = 1:M, .combine = c, .packages = c("igraph", "graphkernels"),
                     .export = c("symmetricize")) %dopar% {
                       G = sample_IRG(p_0)
                       sample.index = sample_index(n, s.size = s_size, replace = TRUE)
                       generate.one.gKSS.sampled(G, C, p_0,sample.index, g.kernel,level)$stats.value
                     }

stopCluster(cl)

critical_value_lower_5 = quantile(statistic_5,alpha/2)
critical_value_upper_5 = quantile(statistic_5,(1-alpha/2))
median(statistic_5)

save.image(paste0("rebuttal/cliques/case_",case_id,".RData"))# path and name of output file for windows

##################

s_size = ceiling(choose(n,2)*0.25)

cores <- detectCores()/2
cl <- makeCluster(cores)
registerDoParallel(cl)

# Parallel sampling of null distribution
statistic_25 <- foreach(i = 1:M, .combine = c, .packages = c("igraph", "graphkernels"),
                       .export = c("symmetricize")) %dopar% {
                         G = sample_IRG(p_0)
                         sample.index = sample_index(n, s.size = s_size, replace = TRUE)
                         generate.one.gKSS.sampled(G, C, p_0,sample.index, g.kernel,level)$stats.value
                       }

stopCluster(cl)

critical_value_lower_25 = quantile(statistic_25,alpha/2)
critical_value_upper_25 = quantile(statistic_25,(1-alpha/2))
median(statistic_25)

save.image(paste0("rebuttal/cliques/case_",case_id,".RData"))# path and name of output file for windows

#######################################################
################### Testing Part ######################
#######################################################

K = seq(0, 25, 1)

l= length(K)

power_gKSS = numeric(l)
power_gKSS_e = numeric(l)
power_gKSS_5 = numeric(l)
power_gKSS_25 = numeric(l)
power_GLR = numeric(l)
power_ST = numeric(l)
total_rej = numeric(l)

for (i in 1:l) {
  
  decision_gKSS = numeric(m)
  decision_gKSS_e = numeric(m)
  decision_gKSS_5 = numeric(m)
  decision_gKSS_25 = numeric(m)
  decision_GLR = numeric(m)
  decision_ST = numeric(m)
  net_rejected = numeric(m)
  
  for (j in 1:m) {
    
    if(K[i] < 2) {
      test_G = sample_IRG(p)} else {
        counter =0
        repeat{
          counter = counter +1
          G = sample_IRG(p)
          if(gsize(G) >= choose(K[i],2)) break
          print(paste0(counter, "tries failed."))
          if(counter > 30) break
        }
        net_rejected[j] = counter-1
        test_G = network_with_clique_nogroup(G, K[i])
        if(gsize(test_G) != gsize(G)) print("Number of edges don't match.")
      }
    
    ########gKSS###########
    
    test_stat_gKSS = generate.one.gKSS(test_G, C, p_0, g.kernel,level)$stats.value
    decision_gKSS[j] = isTRUE(critical_value_lower >= test_stat_gKSS ||test_stat_gKSS >= critical_value_upper)
    #False means do not reject null
    
    #########GLR##########
    
    est = MLE.est(test_G, C)
    Q_h = est$estimates
    E = est$E
    
    #############Test Statistic###################    
    
    test_stat_GLR = -2*(E*log(p_ER) + (choose(n,2) - E)*log(1-p_ER) - E*log(Q_h) - (choose(n,2) - E)*log(1-Q_h))
    decision_GLR[j] = isTRUE(test_stat_GLR > qchisq(1 - alpha, df=1) )
    #False means do not reject null  
    
    ########gKSS with estimated parameters###########
    
    #Q_h = est$estimates
    p_h = matrix(Q_h[C,C], length(C), length(C))
    diag(p_h) = 0
    
    cores <- detectCores()/2
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    # Parallel sampling of null distribution
    statistic_e <- foreach(i = 1:M, .combine = c, .packages = c("igraph", "graphkernels"),
                         .export = c("symmetricize")) %dopar% {
                           g <- sample_IRG(p_h)
                           generate.one.gKSS(g, C, p_h, g.kernel, level)$stats.value
                         }
    
    stopCluster(cl)
    
    critical_value_lower_e = quantile(statistic_e,alpha/2)
    critical_value_upper_e = quantile(statistic_e,(1-alpha/2))
    median(statistic_e)
    
    test_stat_gKSS_e = generate.one.gKSS(test_G, C, p_h, g.kernel,level)$stats.value
    decision_gKSS_e[j] = isTRUE(critical_value_lower_e >= test_stat_gKSS_e ||test_stat_gKSS_e >= critical_value_upper_e)
    #False means do not reject null
    
    ######## Edge resampling ######
    
    s_size = ceiling(choose(n,2)*0.05)
    sample.index_5 = sample_index(n, s.size = s_size, replace = TRUE)
    test_stat_gKSS_5 = generate.one.gKSS.sampled(test_G, C, p_0, sample.index_5, g.kernel,level)$stats.value
    decision_gKSS_5[j] = isTRUE(critical_value_lower_5 >= test_stat_gKSS_5 ||test_stat_gKSS_5 >= critical_value_upper_5)
    #False means do not reject null
    
    s_size = ceiling(choose(n,2)*0.25)
    sample.index_25 = sample_index(n, s.size = s_size, replace = TRUE)
    test_stat_gKSS_25 = generate.one.gKSS.sampled(test_G, C, p_0, sample.index_25, g.kernel,level)$stats.value
    decision_gKSS_25[j] = isTRUE(critical_value_lower_25 >= test_stat_gKSS_25 ||test_stat_gKSS_25 >= critical_value_upper_25)
    #False means do not reject null
    
    #########Spectral test##########
    
    A = as_adjacency_matrix(test_G, type = "both")
    
    #############Test Statistic###################    
    K0 = 1
    test_stat_ST = GoFTestBoot(A, K0, n.dim = K0, n.b = 50, use.boot = T)$test.stat.boot
    decision_ST[j] = isTRUE(max(test_stat_ST) > qtw(alpha / 2, lower.tail = F))
    #False means do not reject null
  }
  total_rej[i] = sum(net_rejected)
  power_gKSS[i] = mean(decision_gKSS)
  power_gKSS_e[i] = mean(decision_gKSS_e)
  power_gKSS_5[i] = mean(decision_gKSS_5)
  power_gKSS_25[i] = mean(decision_gKSS_25)
  power_GLR[i] = mean(decision_GLR)
  power_ST[i] = mean(decision_ST)
  print(i)
  print(power_gKSS_e[i]) 
  
  save(list = c("total_rej","power_gKSS", "power_gKSS_e", "power_gKSS_5", "power_gKSS_25", "power_GLR", "power_ST"), 
       file = paste0(case_id, ".RData"))
}

Q_m = rep(p_ER, length(K))

new_names = c("K", "power", "Q_m", "total_rej", "power_l", "power_u", "Test")

{power.GLR <- data.frame(K, power_GLR, Q_m, total_rej)

pow = power_GLR

power.GLR$power_l = pow - 2*sqrt((pow*(1-pow))/m)
power.GLR$power_l[power.GLR$power_l < 0] <- 0
power.GLR$power_u = pow + 2*sqrt((pow*(1-pow))/m)
power.GLR$power_u[power.GLR$power_u > 1] <- 1

power.GLR$Test <- "GLR"
colnames(power.GLR) =  new_names
}

{power.gKSS <- data.frame(K, power_gKSS, Q_m, total_rej)
  
  pow = power_gKSS
  
  power.gKSS$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.gKSS$power_l[power.gKSS$power_l < 0] <- 0
  power.gKSS$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.gKSS$power_u[power.gKSS$power_u > 1] <- 1
  
  power.gKSS$Test <- "IRG-gKSS"
  colnames(power.gKSS) =  new_names
}

{power.gKSS_e <- data.frame(K, power_gKSS_e, Q_m, total_rej)
  
  pow = power_gKSS_e
  
  power.gKSS_e$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.gKSS_e$power_l[power.gKSS_e$power_l < 0] <- 0
  power.gKSS_e$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.gKSS_e$power_u[power.gKSS_e$power_u > 1] <- 1
  
  power.gKSS_e$Test <- "IRG-gKSS_e"
  colnames(power.gKSS_e) =  new_names
}

{power.gKSS_5 <- data.frame(K, power_gKSS_5, Q_m, total_rej)
  
  pow = power_gKSS_5
  
  power.gKSS_5$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.gKSS_5$power_l[power.gKSS_5$power_l < 0] <- 0
  power.gKSS_5$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.gKSS_5$power_u[power.gKSS_5$power_u > 1] <- 1
  
  power.gKSS_5$Test <- "IRG-gKSS_5"
  colnames(power.gKSS_5) =  new_names
}

{power.gKSS_25 <- data.frame(K, power_gKSS_25, Q_m, total_rej)
  
  pow = power_gKSS_25
  
  power.gKSS_25$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.gKSS_25$power_l[power.gKSS_25$power_l < 0] <- 0
  power.gKSS_25$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.gKSS_25$power_u[power.gKSS_25$power_u > 1] <- 1
  
  power.gKSS_25$Test <- "IRG-gKSS_25"
  colnames(power.gKSS_25) =  new_names
}

{power.ST <- data.frame(K, power_ST, Q_m, total_rej)
  
  pow = power_ST
  
  power.ST$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.ST$power_l[power.ST$power_l < 0] <- 0
  power.ST$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.ST$power_u[power.ST$power_u > 1] <- 1
  
  power.ST$Test <- "ST"
  colnames(power.ST) =  new_names
}


pow = rbind(power.GLR, power.gKSS, power.gKSS_e, power.gKSS_5, power.gKSS_25, power.ST)
#pow = pow %>% filter(K<=8)

legend_order = c("IRG-gKSS", "IRG-gKSS_5", "IRG-gKSS_25", "IRG-gKSS_e", "ST",  "GLR")
my_colours =  c("red", "black" ,"#228833", "#CCBB44", "#009ADE", "#C86DD7")

pow %>% ggplot(aes(x = K, y = power)) + geom_ribbon(aes(ymin = power_l, ymax = power_u, x = K, fill = Test), alpha = 0.2) +
  geom_point(aes(shape = Test, colour = Test)) + geom_line(aes(color = Test)) +
  geom_abline(intercept = alpha, slope = 0, col = "gold") + scale_x_continuous(breaks = seq(0, length(unique(pow$K))-1, by = 2)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="top") +
  ylim(0, 1) + labs(y = "Proportion rejected", x = "K") + scale_color_manual(values=my_colours, breaks=legend_order) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5), breaks=legend_order) + scale_fill_manual(values = my_colours, breaks=legend_order)+
  geom_label(aes(x = K, y = I(1), label = ceiling(total_rej)), angle = 0, size = 3, colour = "black", vjust = 0.2) + coord_cartesian(clip = 'off')
