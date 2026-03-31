######## Methodology ###########

# re-direct to the desired directory
setwd("...") #for windows working directory
source("source-codes.R")

#### For spectral test
library(mgcv)
library(RMTstat)

# and the functions
source('utils.R')

###############################
###### Global Parameters ######
###############################

alpha=0.05
M = 200
m = 50
######### Paramter setting ##############

group_size = 1

##### kernel to use ########

g.kernel= CalculateWLKernel
            #CalculateGraphletKernel 
level = 3  # parameter for the kernel setting

################Model Parameters##############

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


#######################################################
##############Simulations from null model##############
#######################################################

statistic = numeric(M)
for (i in 1:M) {
  G = sample_IRG(p_0)
  statistic[i] = generate.one.gKSS(G, C, p_0, g.kernel,level)$stats.value
}

critical_value_lower = quantile(statistic,alpha/2)
critical_value_upper = quantile(statistic,(1-alpha/2))
median(statistic)

##########################################################################
##Gooddness-of-fit testing for networks simulated from alternative model##
##########################################################################

K = seq(1, 25, 1)

l= length(K)
power_gKSS = numeric(l)
power_gKSS_e = numeric(l)
power_GLR = numeric(l)
power_ST = numeric(l)
total_rej = numeric(l)

for (i in 1:l) {
  
  decision_gKSS = numeric(m)
  
  decision_gKSS_e = numeric(m)
  
  decision_GLR = numeric(m)
  
  decision_ST = numeric(m)
  
  net_rejected = numeric(m)
  
  for (j in 1:m) {
    
    if(K[i] < 2) test_G = sample_IRG(p_0) else {
      counter =0
      repeat{
        counter = counter +1
        G = sample_IRG(p_0)
        test_G = network_with_clique(G, C, K[i])$network
        if(gsize(G) == gsize(test_G)) {
          if(sum(MLE.est(G, C)$E == MLE.est(test_G, C)$E) == L*L) break
        }
        print(paste0(counter, "tries failed."))
        if(counter > 100) break
      }
      net_rejected[j] = counter-1
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
    
    stat_GLR = -2*((E[1,1]*log(Q[1,1]) + (choose(bs[1],2) - E[1,1])*log(1 - Q[1,1]) +
                      E[1,2]*log(Q[1,2]) + (bs[1]*bs[2] - E[1,2])*log(1 - Q[1,2]) +
                      E[2,2]*log(Q[2,2]) + (choose(bs[2],2) - E[2,2])*log(1 - Q[2,2])) - (E[1,1]*log(Q_h[1,1]) + (choose(bs[1],2) - E[1,1])*log(1 - Q_h[1,1]) +
                                                                                            E[1,2]*log(Q_h[1,2]) + (bs[1]*bs[2] - E[1,2])*log(1 - Q_h[1,2]) +
                                                                                            E[2,2]*log(Q_h[2,2]) + (choose(bs[2],2) - E[2,2])*log(1 - Q_h[2,2])))
    
    decision_GLR[j] = isTRUE(stat_GLR > qchisq(1 - alpha, df = 3))
    # df = #free parameters in H1 − #free parameters in H0
    #False means do not reject null  
    
    #### Spectral test
    A = as_adjacency_matrix(test_G, type = "both")
    
    K0 = 2
    stat_ST = GoFTestBoot_truecom(A, K0, n.dim = K0, clusters = C, n.b = 50, use.boot = T)$test.stat.truecom
    decision_ST[j] = isTRUE(stat_ST > qtw(alpha / 2, lower.tail = F))
    #False means do not reject null
    
    ########gKSS with estimated parameters###########
    
    p_h = matrix(Q_h[C,C], length(C), length(C))
    diag(p_h) = 0
    
    statistic_e = numeric(M)
    for (r in 1:M) {
      g <- sample_IRG(p_h)
      statistic_e[r] = generate.one.gKSS(g, C, p_h, g.kernel, level)$stats.value
    }
    
    critical_value_lower_e = quantile(statistic_e,alpha/2)
    critical_value_upper_e = quantile(statistic_e,(1-alpha/2))
    median(statistic_e)
    
    test_stat_gKSS_e = generate.one.gKSS(test_G, C, p_h, g.kernel,level)$stats.value
    decision_gKSS_e[j] = isTRUE(critical_value_lower_e >= test_stat_gKSS_e ||test_stat_gKSS_e >= critical_value_upper_e)
    #False means do not reject null
    
  }
  total_rej[i] = sum(net_rejected)
  power_gKSS[i] = mean(decision_gKSS)
  power_gKSS_e[i] = mean(decision_gKSS_e)
  power_GLR[i] = mean(decision_GLR)
  power_ST[i] = mean(decision_ST)

}

################## Constructing a data frame ######################

if(group_size == 1) n_vec = c(rep("A", length(K)))
if(group_size == 2) n_vec = c(rep("B", length(K)))

new_names = c("K", "power","n_vec", "total_rej", "power_l", "power_u", "Test")

{power.GLR <- data.frame(K, power_GLR, n_vec, total_rej)
  
  pow = power_GLR
  
  power.GLR$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.GLR$power_l[power.GLR$power_l < 0] <- 0
  power.GLR$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.GLR$power_u[power.GLR$power_u > 1] <- 1
  
  power.GLR$Test <- "GLR"
  colnames(power.GLR) =  new_names
}

{power.gKSS <- data.frame(K, power_gKSS, n_vec, total_rej)
  
  pow = power_gKSS
  
  power.gKSS$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.gKSS$power_l[power.gKSS$power_l < 0] <- 0
  power.gKSS$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.gKSS$power_u[power.gKSS$power_u > 1] <- 1
  
  power.gKSS$Test <- "IRG-gKSS"
  colnames(power.gKSS) =  new_names
}

{power.gKSS_e <- data.frame(K, power_gKSS_e, n_vec, total_rej)
  
  pow = power_gKSS_e
  
  power.gKSS_e$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.gKSS_e$power_l[power.gKSS_e$power_l < 0] <- 0
  power.gKSS_e$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.gKSS_e$power_u[power.gKSS_e$power_u > 1] <- 1
  
  power.gKSS_e$Test <- "IRG-gKSS_e"
  colnames(power.gKSS_e) =  new_names
}

{power.ST <- data.frame(K, power_ST, n_vec, total_rej)
  
  pow = power_ST
  
  power.ST$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.ST$power_l[power.ST$power_l < 0] <- 0
  power.ST$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.ST$power_u[power.ST$power_u > 1] <- 1
  
  power.ST$Test <- "ST"
  colnames(power.ST) =  new_names
}

pow = rbind(power.GLR, power.gKSS, power.gKSS_e, power.ST)
#pow = pow %>% filter(K<=11)

legend_order = c("IRG-gKSS", "IRG-gKSS_e", "ST",  "GLR")
my_colours =  c("red", "black" ,"#228833", "#CCBB44", "#009ADE", "#C86DD7")

pow %>% ggplot(aes(x = K, y = power)) + geom_ribbon(aes(ymin = power_l, ymax = power_u, x = K, fill = Test), alpha = 0.2) +
  geom_point(aes(shape = Test, colour = Test)) + geom_line(aes(color = Test)) +
  geom_abline(intercept = alpha, slope = 0, col = "gold") + scale_x_continuous(breaks = seq(0, length(unique(pow$K))-1, by = 2)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="top") +
  ylim(0, 1) + labs(y = "Proportion rejected", x = "K") + scale_color_manual(values=my_colours, breaks=legend_order) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5), breaks=legend_order) + scale_fill_manual(values = my_colours, breaks=legend_order)+
  geom_label(aes(x = K, y = I(1), label = ceiling(total_rej)), angle = 0, size = 3, colour = "black", vjust = 0.2) + coord_cartesian(clip = 'off')

