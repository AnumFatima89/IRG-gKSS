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

group_size = 1

alpha=0.05
M = 200
m = 50

##### kernel to use

g.kernel= CalculateWLKernel
g.kernel_gr= grakel$GraphletSampling(k=3L)

### Parameters

if(group_size == 1) bs = c(8, 22)
if(group_size == 2) bs = c(15, 15)

Q = matrix(
    c(
      0.29, 0.01,
      0.01, 0.22
    ),
    ncol = 2,
    byrow = TRUE
  )

C = c(rep(1,bs[1]), rep(2,bs[2]))

n = sum(bs)
L = length(bs)

p = matrix(Q[C,C], length(C), length(C))
diag(p) = 0

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

stat_h2 = numeric(M)
stat_h3 = numeric(M)
stat_h5 = numeric(M)
stat_h7 = numeric(M)
stat_g3 = numeric(M)
stat_veh = numeric(M)

for (i in 1:M) {
  G = sample_IRG(p_0)
  res = get.P(G, C, p_0)
  P = res$P
  S = res$S.mat
  stat_h2[i] = get.stats(P, S, g.kernel,level=2)
  stat_h3[i] = get.stats(P, S, g.kernel,level=3)
  stat_h5[i] = get.stats(P, S, g.kernel,level=5)
  stat_h7[i] = get.stats(P, S, g.kernel,level=7)
  stat_g3[i] = generate.one.gKSS_grakel(G, C, p_0, g.kernel_gr)$stats.value
  stat_veh[i] = get.stats(P, S, g.kernel = CalculateVertexEdgeHistGaussKernel,level=0.1)

}

cvl_h2 = quantile(stat_h2,alpha/2)
cvu_h2 = quantile(stat_h2,(1-alpha/2))
median(stat_h2)

cvl_h3 = quantile(stat_h3,alpha/2)
cvu_h3 = quantile(stat_h3,(1-alpha/2))
median(stat_h3)

cvl_h5 = quantile(stat_h5,alpha/2)
cvu_h5 = quantile(stat_h5,(1-alpha/2))
median(stat_h5)

cvl_h7 = quantile(stat_h7,alpha/2)
cvu_h7 = quantile(stat_h7,(1-alpha/2))
median(stat_h7)

cvl_g3 = quantile(stat_g3,alpha/2)
cvu_g3 = quantile(stat_g3,(1-alpha/2))
median(stat_g3)

cvl_veh = quantile(stat_veh,alpha/2)
cvu_veh = quantile(stat_veh,(1-alpha/2))
median(stat_veh)

##############################################################
################### Testing Part over k ######################
##############################################################

R = 2
k = seq(2, 15, 1)
l = length(k)

power_gKSS_h2 = numeric(l)
power_gKSS_h3 = numeric(l)
power_gKSS_h5 = numeric(l)
power_gKSS_h7 = numeric(l)
power_gKSS_g3 = numeric(l)
power_gKSS_veh = numeric(l)
mean_max_deg = numeric(l)
median_max_deg = numeric(l)
sd_max_deg = numeric(l)
total_mismatches = numeric(l)
avg_attempts = numeric(l)

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
  
  test_stat.list_gKSS_veh = numeric(m)
  decision_gKSS_veh = numeric(m)
  
  max_deg = numeric(m)
  attempts = numeric(m)
  mismatches = numeric(m)
  
  for (j in 1:m) {
    
    g = sample_IRG(p_0)
    modified = planted_hubs(g, R, k = k[i], C)
    test_G = modified$Graph
    attempts[j] = modified$Iter
    max_deg[j] = max(degree(test_G))
    
    if (gsize(test_G) != gsize(g)) {
      mismatches[j] = mismatches[j] + 1
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
    decision_gKSS_g3[j] = isTRUE(cvl_g3 >= test_stat.list_gKSS_g3[j] ||test_stat.list_gKSS_g3[j] >= cvu_g3)
    #False means do not reject null
    
    test_stat.list_gKSS_veh[j] = get.stats(P, S, g.kernel = CalculateVertexEdgeHistGaussKernel,level=0.1)
    decision_gKSS_veh[j] = isTRUE(cvl_veh >= test_stat.list_gKSS_veh[j] ||test_stat.list_gKSS_veh[j] >= cvu_veh)
    #False means do not reject null
    
    
  }
  power_gKSS_h2[i] = mean(decision_gKSS_h2)
  power_gKSS_h3[i] = mean(decision_gKSS_h3)
  power_gKSS_h5[i] = mean(decision_gKSS_h5)
  power_gKSS_h7[i] = mean(decision_gKSS_h7)
  power_gKSS_g3[i] = mean(decision_gKSS_g3)
  power_gKSS_veh[i] = mean(decision_gKSS_veh)
  mean_max_deg[i] = mean(max_deg) 
  median_max_deg[i] = median(max_deg)
  sd_max_deg[i]  = sd(max_deg)
  total_mismatches[i] = sum(mismatches)
  avg_attempts[i] = mean(attempts)

}

if(group_size == 1) n_vec = c(rep("A", length(k)))
if(group_size == 2) n_vec = c(rep("B", length(k)))

new_names = c("k", "power", "n_vec", "mean_max_deg", "power_l", "power_u", "Test")

{power.h2 <- data.frame(k, power_gKSS_h2, n_vec, mean_max_deg)
  
  pow = power_gKSS_h2
  
  power.h2$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.h2$power_l[power.h2$power_l < 0] <- 0
  power.h2$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.h2$power_u[power.h2$power_u > 1] <- 1
  
  power.h2$Test <- "WL2"
  colnames(power.h2) =  new_names
}

{power.h3 <- data.frame(k, power_gKSS_h3, n_vec, mean_max_deg)
  
  pow = power_gKSS_h3
  
  power.h3$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.h3$power_l[power.h3$power_l < 0] <- 0
  power.h3$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.h3$power_u[power.h3$power_u > 1] <- 1
  
  power.h3$Test <- "WL3"
  colnames(power.h3) =  new_names
}

{power.h5 <- data.frame(k, power_gKSS_h5, n_vec, mean_max_deg)
  
  pow = power_gKSS_h5
  
  power.h5$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.h5$power_l[power.h5$power_l < 0] <- 0
  power.h5$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.h5$power_u[power.h5$power_u > 1] <- 1
  
  power.h5$Test <- "WL5"
  colnames(power.h5) =  new_names
}

{power.h7 <- data.frame(k, power_gKSS_h7, n_vec, mean_max_deg)
  
  pow = power_gKSS_h7
  
  power.h7$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.h7$power_l[power.h7$power_l < 0] <- 0
  power.h7$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.h7$power_u[power.h7$power_u > 1] <- 1
  
  power.h7$Test <- "WL7"
  colnames(power.h7) =  new_names
}

{power.g3 <- data.frame(k, power_gKSS_g3, n_vec, mean_max_deg)
  
  pow = power_gKSS_g3
  
  power.g3$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.g3$power_l[power.g3$power_l < 0] <- 0
  power.g3$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.g3$power_u[power.g3$power_u > 1] <- 1
  
  power.g3$Test <- "G3"
  colnames(power.g3) =  new_names
}

{power.veh <- data.frame(k, power_gKSS_veh, n_vec, mean_max_deg)
  
  pow = power_gKSS_veh
  
  power.veh$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.veh$power_l[power.veh$power_l < 0] <- 0
  power.veh$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.veh$power_u[power.veh$power_u > 1] <- 1
  
  power.veh$Test <- "VEH"
  colnames(power.veh) =  new_names
}

power.df_k = rbind(power.h2, power.h3, power.h5, power.h7, power.g3, power.veh)

pow = power.df_k

##############################################################
################### Testing Part over R ######################
##############################################################

R = seq(0, 10, 1)
k = 3
l = length(R)

power_gKSS_h2 = numeric(l)
power_gKSS_h3 = numeric(l)
power_gKSS_h5 = numeric(l)
power_gKSS_h7 = numeric(l)
power_gKSS_g3 = numeric(l)
power_gKSS_veh = numeric(l)
mean_max_deg = numeric(l)
median_max_deg = numeric(l)
sd_max_deg = numeric(l)
total_mismatches = numeric(l)
avg_attempts = numeric(l)

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
  
  test_stat.list_gKSS_veh = numeric(m)
  decision_gKSS_veh = numeric(m)
  
  max_deg = numeric(m)
  attempts = numeric(m)
  mismatches = numeric(m)
  
  for (j in 1:m) {
    
    g = sample_IRG(p_0)
    if(R[i] == 0) {test_G = g
    attempts[j] = 0} else { modified = planted_hubs(g, R[i], k, C)
    test_G =  modified$Graph
    attempts[j] = modified$Iter
    }
    max_deg[j] = max(degree(test_G))
    
    if (gsize(test_G) != gsize(g)) {
      mismatches = mismatches + 1
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
    decision_gKSS_g3[j] = isTRUE(cvl_g3 >= test_stat.list_gKSS_g3[j] ||test_stat.list_gKSS_g3[j] >= cvu_g3)
    #False means do not reject null
    
    test_stat.list_gKSS_veh[j] = get.stats(P, S, g.kernel = CalculateVertexEdgeHistGaussKernel,level=0.1)
    decision_gKSS_veh[j] = isTRUE(cvl_veh >= test_stat.list_gKSS_veh[j] ||test_stat.list_gKSS_veh[j] >= cvu_veh)
    #False means do not reject null
    
    
  }
  power_gKSS_h2[i] = mean(decision_gKSS_h2)
  power_gKSS_h3[i] = mean(decision_gKSS_h3)
  power_gKSS_h5[i] = mean(decision_gKSS_h5)
  power_gKSS_h7[i] = mean(decision_gKSS_h7)
  power_gKSS_g3[i] = mean(decision_gKSS_g3)
  power_gKSS_veh[i] = mean(decision_gKSS_veh)
  mean_max_deg[i] = mean(max_deg) 
  median_max_deg[i] = median(max_deg)
  sd_max_deg[i]  = sd(max_deg)
  total_mismatches[i] = sum(mismatches)
  avg_attempts[i] = mean(attempts)

}

if(group_size == 1) n_vec = c(rep("A", length(R)))
if(group_size == 2) n_vec = c(rep("B", length(R)))

new_names = c("R", "power", "n_vec", "mean_max_deg", "power_l", "power_u", "Test")

{power.h2 <- data.frame(R, power_gKSS_h2, n_vec, mean_max_deg)
  
  pow = power_gKSS_h2
  
  power.h2$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.h2$power_l[power.h2$power_l < 0] <- 0
  power.h2$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.h2$power_u[power.h2$power_u > 1] <- 1
  
  power.h2$Test <- "WL2"
  colnames(power.h2) =  new_names
}

{power.h3 <- data.frame(R, power_gKSS_h3, n_vec, mean_max_deg)
  
  pow = power_gKSS_h3
  
  power.h3$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.h3$power_l[power.h3$power_l < 0] <- 0
  power.h3$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.h3$power_u[power.h3$power_u > 1] <- 1
  
  power.h3$Test <- "WL3"
  colnames(power.h3) =  new_names
}

{power.h5 <- data.frame(R, power_gKSS_h5, n_vec, mean_max_deg)
  
  pow = power_gKSS_h5
  
  power.h5$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.h5$power_l[power.h5$power_l < 0] <- 0
  power.h5$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.h5$power_u[power.h5$power_u > 1] <- 1
  
  power.h5$Test <- "WL5"
  colnames(power.h5) =  new_names
}

{power.h7 <- data.frame(R, power_gKSS_h7, n_vec, mean_max_deg)
  
  pow = power_gKSS_h7
  
  power.h7$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.h7$power_l[power.h7$power_l < 0] <- 0
  power.h7$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.h7$power_u[power.h7$power_u > 1] <- 1
  
  power.h7$Test <- "WL7"
  colnames(power.h7) =  new_names
}

{power.g3 <- data.frame(R, power_gKSS_g3, n_vec, mean_max_deg)
  
  pow = power_gKSS_g3
  
  power.g3$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.g3$power_l[power.g3$power_l < 0] <- 0
  power.g3$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.g3$power_u[power.g3$power_u > 1] <- 1
  
  power.g3$Test <- "G3"
  colnames(power.g3) =  new_names
}

{power.veh <- data.frame(R, power_gKSS_veh, n_vec, mean_max_deg)
  
  pow = power_gKSS_veh
  
  power.veh$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.veh$power_l[power.veh$power_l < 0] <- 0
  power.veh$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.veh$power_u[power.veh$power_u > 1] <- 1
  
  power.veh$Test <- "VEH"
  colnames(power.veh) =  new_names
}

power.df_R = rbind(power.h2, power.h3, power.h5, power.h7, power.g3, power.veh)

save.image(paste0("rebuttal/hubs/kernel_comp.RData"))# path and name of output file for windows

legend_order = c("WL2", "WL3", "WL5", "WL7", "G3", "VEH")
my_colours =  c("red", "black" ,"#228833", "#CCBB44", "#009ADE", "#C86DD7")

p2 = power.df_k %>% ggplot(aes(x = k, y = power)) + geom_ribbon(aes(ymin = power_l, ymax = power_u, x = k, fill = Test), alpha = 0.2) +
  geom_point(aes(shape = Test, colour = Test)) + geom_line(aes(color = Test)) +
  geom_abline(intercept = alpha, slope = 0, col = "gold") + scale_x_continuous(breaks = seq(2, 15, by = 2)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="top", axis.title.y=element_blank()) +
  ylim(0, 1) + labs( x = "k") + scale_color_manual(values=my_colours, breaks=legend_order) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5), breaks=legend_order) + scale_fill_manual(values = my_colours, breaks=legend_order)+
  geom_label(aes(x = k, y = I(1), label = ceiling(mean_max_deg)), angle = 0, size = 3, colour = "black", vjust = 0.2) + coord_cartesian(clip = 'off')


p1 = power.df_R %>% ggplot(aes(x = R, y = power)) + geom_ribbon(aes(ymin = power_l, ymax = power_u, x = R, fill = Test), alpha = 0.2) +
  geom_point(aes(shape = Test, colour = Test)) + geom_line(aes(color = Test)) +
  geom_abline(intercept = alpha, slope = 0, col = "gold") + scale_x_continuous(breaks = seq(0, length(unique(power.df_R$R))-1, by = 2)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="top", axis.title.y=element_blank()) +
  ylim(0, 1) + labs( x = "R") + scale_color_manual(values=my_colours, breaks=legend_order) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5), breaks=legend_order) + scale_fill_manual(values = my_colours, breaks=legend_order)+
  geom_label(aes(x = R, y = I(1), label = ceiling(mean_max_deg)), angle = 0, size = 3, colour = "black", vjust = 0.2) + coord_cartesian(clip = 'off')


multi <- ggarrange(p1,p2, #plots that are going to be included in this multipanel figure
                   #labels = c("A", "B", "C","D"), #labels given each panel 
                   ncol = 2, nrow = 1, #adjust plot space 
                   common.legend = T, widths = c(0.8, 1)) 
annotate_figure(multi, 
                left = "Proportion rejected")
