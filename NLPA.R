# re-direct to the desired directory
setwd("...") #for windows working directory
source("source-codes.R")
source("source-code_graphlet_py.R")

#### For spectral test
library(mgcv)
library(RMTstat)

# and the functions
source('utils.R')

###################################################

n = 20
power = 1
m.links = 1


########### Testing setup ########

alpha = 0.05
M = 200
m=50

########### Input graph ############

G = sample_pa(n, power, m=m.links, directed=FALSE)
p_h = gsize(G)/choose(gorder(G),2)
C = rep(1, n)

################ Null Model Parameters ##############

p = matrix(p_h, n, n)
diag(p) = 0
p_0 = p

##############Simulations from null model ERMM-gKSS ##############

statistic_WL2 = numeric(M)

for (i in 1:M) {
  g <- sample_IRG(p_0)
  statistic_WL2[i] = generate.one.gKSS(g, C, p_0, CalculateWLKernel, 2)$stats.value
}

critical_value_lower_WL2 = quantile(statistic_WL2, alpha/2)
critical_value_upper_WL2 = quantile(statistic_WL2, (1-alpha/2))
median(statistic_WL2)

###

statistic_WL3 = numeric(M)

for (i in 1:M) {
  g <- sample_IRG(p_0)
  statistic_WL3[i] = generate.one.gKSS(g, C, p_0, CalculateWLKernel, 3)$stats.value
}


critical_value_lower_WL3 = quantile(statistic_WL3, alpha/2)
critical_value_upper_WL3 = quantile(statistic_WL3, (1-alpha/2))
median(statistic_WL3)

###

g.kernel= grakel$GraphletSampling(k=3L)

statistic_G = c()
for (i in 1:M){
  g = sample_IRG(p_0)
  statistic_G[i] = generate.one.gKSS_grakel(g, C, p_0, g.kernel)$stats.value
}

critical_value_lower_G = quantile(statistic_G, alpha/2)
critical_value_upper_G = quantile(statistic_G, (1-alpha/2))
median(statistic_G)

############ Testing networks simulated using BA model ###########

decision_gKSS_WL2 = numeric(m)
decision_gKSS_WL3 = numeric(m)
decision_gKSS_G = numeric(m)
decision_ST = numeric(m)

for (j in 1:m) {
  
  test_G = sample_pa(n, power, m=m.links, directed=FALSE)
  
  ########gKSS###########

  test_stat_gKSS_WL2 = generate.one.gKSS(test_G, C, p_0, CalculateWLKernel, 2)$stats.value
  decision_gKSS_WL2[j] = isTRUE(critical_value_lower_WL2 >= test_stat_gKSS_WL2 ||test_stat_gKSS_WL2 >= critical_value_upper_WL2)
  #False means do not reject null
  
  test_stat_gKSS_WL3 = generate.one.gKSS(test_G, C, p_0, CalculateWLKernel, 3)$stats.value
  decision_gKSS_WL3[j] = isTRUE(critical_value_lower_WL3 >= test_stat_gKSS_WL3 ||test_stat_gKSS_WL3 >= critical_value_upper_WL3)
  #False means do not reject null
  
  test_stat_gKSS_G = generate.one.gKSS_grakel(test_G, C, p_0, g.kernel)$stats.value
  decision_gKSS_G[j] = isTRUE(critical_value_lower_G >= test_stat_gKSS_G ||test_stat_gKSS_G >= critical_value_upper_G)
  #False means do not reject null
  
  #### Spectral test
  A = as_adjacency_matrix(test_G, type = "both")
  K0 = 1
  
  stat_ST = GoFTestBoot_truecom(A, K0, n.dim = K0, clusters = C, n.b = 50, use.boot = T)$test.stat.truecom
  decision_ST[j] = isTRUE(stat_ST > qtw(alpha / 2, lower.tail = F))
  #False means do not reject null
  
  print(j)
}

power_gKSS_WL2 = mean(decision_gKSS_WL2)
power_gKSS_WL3 = mean(decision_gKSS_WL3)
power_gKSS_G = mean(decision_gKSS_G)
power_ST = mean(decision_ST)

print(power_gKSS_WL2)
print(power_gKSS_WL3)
print(power_gKSS_G)
print(power_ST)
