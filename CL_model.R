######## Methodology ###########

# re-direct to the desired directory
setwd("...") #for windows working directory
source("source-codes.R")
source("source-code_graphlet_py.R")

#### For spectral test
library(mgcv)
library(RMTstat)
source('utils.R')

calculate_cl_probability_matrix <- function(weights, cap_at_one = TRUE) {
  
  # 1. Calculate the sum of all weights (S)
  S <- sum(weights)
  
  # 2. Calculate the numerator: w_i * w_j for all pairs (using outer product)
  # The outer() function returns an N x N matrix where element [i, j] = weights[i] * weights[j]
  numerator_matrix <- outer(weights, weights)
  
  # 3. Divide by the sum of weights (S) to get the probability matrix P_ij
  probability_matrix <- numerator_matrix / S
  
  # 4. Apply the minimum constraint (P_ij <= 1)
  if (cap_at_one) {
    probability_matrix[probability_matrix > 1] <- 1
  }
  
  # Set the diagonal elements to 0 if loops are not allowed.
  # The question implies edge probability between distinct nodes, so we set P_ii = 0
  #diag(probability_matrix) <- 0
  
  return(probability_matrix)
}

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
n=30

###############################

min_degree <- 3

# Maximum expected degree
max_degree <- 10

e <- 3  #### e = 0 means we are testing true null

# --- Generate Uniform Expected Degrees ---

# The runif() function draws 'n' random numbers 
# uniformly between 'min_degree' and 'max_degree'.
expected_degrees <- round(runif(
  n, 
  min = min_degree, 
  max = max_degree
),0)

# Check the sum of weights (Expected total degree: 80*5 + 20*15 = 700)
sum(expected_degrees)

#expected_degrees[!expected_degrees^2 < sum(expected_degrees)]
expected_degrees[expected_degrees^2 > sum(expected_degrees)]

P = calculate_cl_probability_matrix(expected_degrees, cap_at_one = TRUE)
P[P==0]

p_0 = P
diag(p_0) <- 0


G = sample_chung_lu(out.weights = expected_degrees, loops = FALSE)
plot(G)

C = sample(1,n,replace = TRUE)

### Null set calculated using networks simulated from IRG with edge prob matrix p_0

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
}

critical_value_lower_G = quantile(statistic_G,alpha/2)
critical_value_upper_G = quantile(statistic_G,(1-alpha/2))
median(statistic_G)

###### Testing Part #######

decision_SN = numeric(m)
decision_gKSS = numeric(m)
decision_gKSS_G = numeric(m)
decision_gKSS_e = numeric(m)
decision_gKSS_Ge = numeric(m)

for (j in 1:m) {
  
  test_G = sample_chung_lu(out.weights = expected_degrees + e, loops = FALSE)
  
  #########Spectral test##########
  
  A = as_adjacency_matrix(test_G, type = "both")
  
  K0 = 1
  test_stat_SN = GoFTestBoot_truecom(A, K0, n.dim = K0, clusters = C, n.b = 50, use.boot = T)$test.stat.truecom
  decision_SN[j] = isTRUE(test_stat_SN > qtw(alpha / 2, lower.tail = F)) #qtw(0.025, beta=1, lower.tail = FALSE, log.p = FALSE))
  #False means do not reject null
  
  ########gKSS_full###########
  
  test_stat_gKSS = generate.one.gKSS(test_G, C, p_0, g.kernel,level)$stats.value
  decision_gKSS[j] = isTRUE(critical_value_lower >= test_stat_gKSS ||test_stat_gKSS >= critical_value_upper)
  #False means do not reject null
  
  test_stat_gKSS_G = generate.one.gKSS_grakel(test_G, C, p_0, g.kernel_gr)$stats.value
  decision_gKSS_G[j] = isTRUE(critical_value_lower_G >= test_stat_gKSS_G ||test_stat_gKSS_G >= critical_value_upper_G)
  #False means do not reject null
  
  ########gKSS_est###########
  
  obs_deg = degree(test_G)
  obs_deg[obs_deg^2 > sum(obs_deg)]
  
  p_h = calculate_cl_probability_matrix(obs_deg, cap_at_one = TRUE)
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

power_SN = mean(decision_SN)
power_gKSS = mean(decision_gKSS)
power_gKSS_G = mean(decision_gKSS_G)
power_gKSS_e = mean(decision_gKSS_e)
power_gKSS_Ge = mean(decision_gKSS_Ge)


print(power_SN)
print(power_gKSS)
print(power_gKSS_e)
print(power_gKSS_G)
print(power_gKSS_Ge)
