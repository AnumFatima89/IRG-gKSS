######## Methodology ###########

# re-direct to the desired directory
setwd("...") #for windows working directory
source("source-codes.R")
source("source-code_graphlet_py.R")

#### For spectral test
library(mgcv)
library(RMTstat)
source('utils.R')
###############################
###### Global Parameters ######
###############################

alpha=0.05
n = 30
C = rep(1,n)

M = 200
m = 50

##### kernel to use

g.kernel= CalculateWLKernel
g.kernel_gr= grakel$GraphletSampling(k=3L)
g.kernel_gr5= grakel$GraphletSampling(k=5L)


p_ER = 0.06
p = matrix(p_ER, n, n)
diag(p) = 0

################Null Parameter set##############

p_0 = p

#### Usefull function #####

get.P = function(G, C, p){
  V(G)$name = C
  E(G)$label = edge.labels(G, C)
  
  X = as_adjacency_matrix(G, type = "both")
  S = X[upper.tri(X)]
  n = length(S)
  S.t =  p[upper.tri(p)]
  S.t.vec = abs(S - S.t)
  S.mat = S.t.vec %*% t(S.t.vec)
  
  P = compute.transition.list(X)
  list(P = P, S.mat = S.mat)
}

get.stats = function(P, S, g.kernel=CalculateWLKernel,level=3){
  K = compute.kernel(P, g.kernel, level)
  
  J.kernel = S * K
  
  stats.value = mean(J.kernel)
  return(stats.value)
}

#######################################################
##############Simulations from null model##############
#######################################################

stats_h2 = numeric(M)
stats_h3 = numeric(M)
stats_h5 = numeric(M)
stats_h7 = numeric(M)
stats_g3 = numeric(M)
stats_g5 = numeric(M)
stats_veh = numeric(M)

for (i in 1:M) {
  G = sample_IRG(p_0)
  res = get.P(G, C, p_0)
  P = res$P
  S = res$S.mat
  
  stats_h2[i] = get.stats(P, S, g.kernel,level=2)
  stats_h3[i] = get.stats(P, S, g.kernel,level=3)
  stats_h5[i] = get.stats(P, S, g.kernel,level=5)
  stats_h7[i] = get.stats(P, S, g.kernel,level=7)
  stats_g3[i] = generate.one.gKSS_grakel(G, C, p_0, g.kernel_gr)$stats.value
  stats_g5[i] = generate.one.gKSS_grakel(G, C, p_0, g.kernel_gr5)$stats.value
  stats_veh[i] = get.stats(P, S, g.kernel = CalculateVertexEdgeHistGaussKernel,level=0.1)

}

cvl_h2 = quantile(stats_h2,alpha/2)
cvu_h2 = quantile(stats_h2,(1-alpha/2))
median(stats_h2)

cvl_h3 = quantile(stats_h3,alpha/2)
cvu_h3 = quantile(stats_h3,(1-alpha/2))
median(stats_h3)

cvl_h5 = quantile(stats_h5,alpha/2)
cvu_h5 = quantile(stats_h5,(1-alpha/2))
median(stats_h5)

cvl_h7 = quantile(stats_h7,alpha/2)
cvu_h7 = quantile(stats_h7,(1-alpha/2))
median(stats_h7)

cvl_g3 = quantile(stats_g3,alpha/2)
cvu_g3 = quantile(stats_g3,(1-alpha/2))
median(stats_g3)

cvl_g5 = quantile(stats_g5,alpha/2)
cvu_g5 = quantile(stats_g5,(1-alpha/2))
median(stats_g5)

cvl_veh = quantile(stats_veh,alpha/2)
cvu_veh = quantile(stats_veh,(1-alpha/2))
median(stats_veh)

#######################################################
################### Testing Part ######################
#######################################################

K = seq(0, 15, 1)

l= length(K)

power_gKSS_h2 = numeric(l)
power_gKSS_h3 = numeric(l)
power_gKSS_h5 = numeric(l)
power_gKSS_h7 = numeric(l)
power_gKSS_g3 = numeric(l)
power_gKSS_g5 = numeric(l)
power_gKSS_veh = numeric(l)
total_rej = numeric(l)

for (i in 1:l) {
  
  test_stat.list_gKSS_h2 = numeric(m)
  decision_gKSS_h2 = numeric(m)
  
  test_stat.list_gKSS_h3 = numeric(m)
  decision_gKSS_h3 = numeric(m)
  
  test_stat.list_gKSS_h5 = numeric(m)
  decision_gKSS_h5 = numeric(m)
  
  test_stat.list_gKSS_h7 = numeric(m)
  decision_gKSS_h7 = numeric(m)
  
  test_stat.list_gKSS_g3 = numeric(m)
  decision_gKSS_g3 = numeric(m)
  
  test_stat.list_gKSS_g5 = numeric(m)
  decision_gKSS_g5 = numeric(m)
  
  test_stat.list_gKSS_veh = numeric(m)
  decision_gKSS_veh = numeric(m)
  
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
    
    res = get.P(test_G, C, p_0)
    P = res$P
    S = res$S.mat
    
    test_stat.list_gKSS_h2[j] = get.stats(P, S, g.kernel,level=2)
    decision_gKSS_h2[j] = isTRUE(cvl_h2 >= test_stat.list_gKSS_h2[j] ||test_stat.list_gKSS_h2[j] >= cvu_h2)
    #False means do not reject null
    
    test_stat.list_gKSS_h3[j] = get.stats(P, S, g.kernel,level=3)
    decision_gKSS_h3[j] = isTRUE(cvl_h3 >= test_stat.list_gKSS_h3[j] ||test_stat.list_gKSS_h3[j] >= cvu_h3)
    #False means do not reject null
    
    test_stat.list_gKSS_h5[j] = get.stats(P, S, g.kernel,level=5)
    decision_gKSS_h5[j] = isTRUE(cvl_h5 >= test_stat.list_gKSS_h5[j] ||test_stat.list_gKSS_h5[j] >= cvu_h5)
    #False means do not reject null
    
    test_stat.list_gKSS_h7[j] = get.stats(P, S, g.kernel,level=7)
    decision_gKSS_h7[j] = isTRUE(cvl_h7 >= test_stat.list_gKSS_h7[j] ||test_stat.list_gKSS_h7[j] >= cvu_h7)
    #False means do not reject null
    
    test_stat.list_gKSS_g3[j] = generate.one.gKSS_grakel(test_G, C, p_0, g.kernel_gr)$stats.value
                                #get.stats(P, S, g.kernel = CalculateGraphletKernel,level=3) 
    decision_gKSS_g3[j] = isTRUE(cvl_g3 >= test_stat.list_gKSS_g3[j] ||test_stat.list_gKSS_g3[j] >= cvu_g3)
    #False means do not reject null
    
    test_stat.list_gKSS_g5[j] = generate.one.gKSS_grakel(test_G, C, p_0, g.kernel_gr5)$stats.value
                                #get.stats(P, S, g.kernel = CalculateGraphletKernel,level=5) 
    decision_gKSS_g5[j] = isTRUE(cvl_g5 >= test_stat.list_gKSS_g5[j] ||test_stat.list_gKSS_g5[j] >= cvu_g5)
    #False means do not reject null
    
    test_stat.list_gKSS_veh[j] = get.stats(P, S, g.kernel = CalculateVertexEdgeHistGaussKernel,level=0.1)
    decision_gKSS_veh[j] = isTRUE(cvl_veh >= test_stat.list_gKSS_veh[j] ||test_stat.list_gKSS_veh[j] >= cvu_veh)
    #False means do not reject null
    
    
  }
  total_rej[i] = sum(net_rejected)
  power_gKSS_h2[i] = mean(decision_gKSS_h2)
  power_gKSS_h3[i] = mean(decision_gKSS_h3)
  power_gKSS_h5[i] = mean(decision_gKSS_h5)
  power_gKSS_h7[i] = mean(decision_gKSS_h7)
  power_gKSS_g3[i] = mean(decision_gKSS_g3)
  power_gKSS_g5[i] = mean(decision_gKSS_g5)
  power_gKSS_veh[i] = mean(decision_gKSS_veh)

}

Q_m = rep(p_ER, length(K))

new_names = c("K", "power", "Q_m", "total_rej", "power_l", "power_u", "Test")

{power.h2 <- data.frame(K, power_gKSS_h2, Q_m, total_rej)
  
  pow = power_gKSS_h2
  
  power.h2$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.h2$power_l[power.h2$power_l < 0] <- 0
  power.h2$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.h2$power_u[power.h2$power_u > 1] <- 1
  
  power.h2$Test <- "WL2"
  colnames(power.h2) =  new_names
}

{power.h3 <- data.frame(K, power_gKSS_h3, Q_m, total_rej)
  
  pow = power_gKSS_h3
  
  power.h3$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.h3$power_l[power.h3$power_l < 0] <- 0
  power.h3$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.h3$power_u[power.h3$power_u > 1] <- 1
  
  power.h3$Test <- "WL3"
  colnames(power.h3) =  new_names
}

{power.h5 <- data.frame(K, power_gKSS_h5, Q_m, total_rej)
  
  pow = power_gKSS_h5
  
  power.h5$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.h5$power_l[power.h5$power_l < 0] <- 0
  power.h5$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.h5$power_u[power.h5$power_u > 1] <- 1
  
  power.h5$Test <- "WL5"
  colnames(power.h5) =  new_names
}

{power.h7 <- data.frame(K, power_gKSS_h7, Q_m, total_rej)
  
  pow = power_gKSS_h7
  
  power.h7$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.h7$power_l[power.h7$power_l < 0] <- 0
  power.h7$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.h7$power_u[power.h7$power_u > 1] <- 1
  
  power.h7$Test <- "WL7"
  colnames(power.h7) =  new_names
}

{power.g3 <- data.frame(K, power_gKSS_g3, Q_m, total_rej)
  
  pow = power_gKSS_g3
  
  power.g3$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.g3$power_l[power.g3$power_l < 0] <- 0
  power.g3$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.g3$power_u[power.g3$power_u > 1] <- 1
  
  power.g3$Test <- "G3"
  colnames(power.g3) =  new_names
}

{power.g5 <- data.frame(K, power_gKSS_g5, Q_m, total_rej)
  
  pow = power_gKSS_g5
  
  power.g5$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.g5$power_l[power.g5$power_l < 0] <- 0
  power.g5$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.g5$power_u[power.g5$power_u > 1] <- 1
  
  power.g5$Test <- "G5"
  colnames(power.g5) =  new_names
}

{power.veh <- data.frame(K, power_gKSS_veh, Q_m, total_rej)
  
  pow = power_gKSS_veh
  
  power.veh$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.veh$power_l[power.veh$power_l < 0] <- 0
  power.veh$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.veh$power_u[power.veh$power_u > 1] <- 1
  
  power.veh$Test <- "VEH"
  colnames(power.veh) =  new_names
}

pow = rbind(power.h2, power.h3, power.h5, power.h7, power.g3, power.g5, power.veh)
#pow = pow %>% filter(K<=8)

legend_order = c("WL2", "WL3", "WL5", "WL7", "G3", "G5", "VEH")
my_colours =  c("red", "black" ,"#228833", "#CCBB44", "#009ADE","#E16A86", "#C86DD7")

pow %>% ggplot(aes(x = K, y = power)) + geom_ribbon(aes(ymin = power_l, ymax = power_u, x = K, fill = Test), alpha = 0.2) +
  geom_point(aes(shape = Test, colour = Test)) + geom_line(aes(color = Test)) +
  geom_abline(intercept = alpha, slope = 0, col = "gold") + scale_x_continuous(breaks = seq(0, length(unique(pow$K))-1, by = 2)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="top") +
  ylim(0, 1) + labs(y = "Proportion rejected", x = "K") + scale_color_manual(values=my_colours, breaks=legend_order) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6), breaks=legend_order) + scale_fill_manual(values = my_colours, breaks=legend_order)+
  geom_label(aes(x = K, y = I(1), label = ceiling(total_rej)), angle = 0, size = 3, colour = "black", vjust = 0.2) + coord_cartesian(clip = 'off')

save.image(paste0("rebuttal/cliques/kernel_comp.RData"))# path and name of output file for windows

"#E16A86" "#B88A00" "#50A315" "#00AD9A" "#009ADE" "#C86DD7"