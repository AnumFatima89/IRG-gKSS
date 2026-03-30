######## Methodology ###########

# re-direct to the desired directory
setwd("...") #for windows working directory
source("source-codes.R")
source("source-code_graphlet_py.R")

#### For spectral test
library(mgcv)
library(RMTstat)
source('utils.R')

##### kernel to use ########

g.kernel= CalculateWLKernel
level = 3  # parameter for the kernel setting

g.kernel_gr= grakel$GraphletSampling(k=3L)

###############################
###### Global Parameters ######
###############################

alpha=0.05
M = 200
m = 50
n=27

###############################

a = matrix(
    c(
      log(0.6), log(0.1),
      log(0.1), log(0.3)
    ),
    ncol = 2,
    byrow = TRUE
  )

# a = matrix(
#     c(
#       log(0.6), log(0.2),
#       log(0.2), log(0.6)
#     ),
#     ncol = 2,
#     byrow = TRUE
#   )

beta <- log(runif(n, 0, 1))

C = sample(c(1,2),n,replace = TRUE)
g = outer(beta, beta, `+`) + matrix(a[C,C], length(C), length(C)) 
p = exp(g) 
diag(p) = 0

############### Testing the fit of an ERMM ###########

g_0 = matrix(a[C,C], length(C), length(C)) 
p_0 = 0.25*exp(g_0)
diag(p_0) = 0

statistic = numeric(M)

for (i in 1:M) {
  g <- sample_IRG(p_0)
  statistic[i] = generate.one.gKSS(g, C, p_0, g.kernel, level)$stats.value
}

critical_value_lower = quantile(statistic,alpha/2)
critical_value_upper = quantile(statistic,(1-alpha/2))
median(statistic)

# For grakel
statistic_G = c()
for (i in 1:M) {
  g <- sample_IRG(p_0)
  statistic_G[i] = generate.one.gKSS_grakel(g, C, p_0, g.kernel_gr)$stats.value
  print(i)
}

critical_value_lower_G = quantile(statistic_G,alpha/2)
critical_value_upper_G = quantile(statistic_G,(1-alpha/2))
median(statistic_G)

####### Testing Part #########

decision_ST = numeric(m)
decision_gKSS = numeric(m)
decision_gKSS_e = numeric(m)
decision_gKSS_G = numeric(m)
decision_gKSS_Ge = numeric(m)

for (j in 1:m) {
  
  test_G <- sample_IRG(p)
  
  #########Spectral test##########
  
  A = as_adjacency_matrix(test_G, type = "both")
  
  K0 = 2
  test_stat_ST = GoFTestBoot_truecom(A, K0, n.dim = K0, clusters = C, n.b = 50, use.boot = T)$test.stat.truecom #GoFTestBoot(A, K0, n.dim = K0, n.b = 50, use.boot = T)$test.stat.boot[1]
  decision_ST[j] = isTRUE(test_stat_ST > qtw(alpha / 2, lower.tail = F)) #qtw(0.025, beta=1, lower.tail = FALSE, log.p = FALSE))
  #False means do not reject null
  
  ########gKSS_full###########
  
  test_stat_gKSS = generate.one.gKSS(test_G, C, p_0, g.kernel,level)$stats.value
  decision_gKSS[j] = isTRUE(critical_value_lower >= test_stat_gKSS ||test_stat_gKSS >= critical_value_upper)
  #False means do not reject null
  
  test_stat_gKSS_G = generate.one.gKSS_grakel(test_G, C, p_0, g.kernel_gr)$stats.value
  decision_gKSS_G[j] = isTRUE(critical_value_lower_G >= test_stat_gKSS_G ||test_stat_gKSS_G >= critical_value_upper_G)
  #False means do not reject null
  
  ########gKSS_estimated###########
  
  est = MLE.est(test_G, C)
  Q_h = est$estimates
  E = est$E
  
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
  
  # For grakel
  statistic_Ge = c()
  for (r in 1:M) {
    g <- sample_IRG(p_h)
    statistic_Ge[r] = generate.one.gKSS_grakel(g, C, p_h, g.kernel_gr)$stats.value
  }
  
  critical_value_lower_Ge = quantile(statistic_Ge,alpha/2)
  critical_value_upper_Ge = quantile(statistic_Ge,(1-alpha/2))
  median(statistic_Ge)
  
  test_stat_gKSS_Ge = generate.one.gKSS_grakel(test_G, C, p_h, g.kernel_gr)$stats.value
  decision_gKSS_Ge[j] = isTRUE(critical_value_lower_Ge >= test_stat_gKSS_Ge ||test_stat_gKSS_Ge >= critical_value_upper_Ge)
  #False means do not reject null
  
}

power_ST = mean(decision_ST)
power_gKSS = mean(decision_gKSS)
power_gKSS_e = mean(decision_gKSS_e)
power_gKSS_G = mean(decision_gKSS_G)
power_gKSS_Ge = mean(decision_gKSS_Ge)


print(power_ST)
print(power_gKSS)
print(power_gKSS_e)
print(power_gKSS_G)
print(power_gKSS_Ge)
