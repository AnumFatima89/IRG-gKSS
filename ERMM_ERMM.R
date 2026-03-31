setwd("...") #for windows working directory
source("source-codes.R")

###############################
###### Global Parameters ######
###############################

alpha=0.05
M = 200
m = 50
######### Paramter setting ##############

Q_mat = 1
group_size = 1

##### kernel to use ########

g.kernel= CalculateWLKernel
          #CalculateGraphletKernel 
level = 3  # parameter for the kernel setting

################Model Parameters##############

if(group_size == 1) bs = c(10, 40)
if(group_size == 2) bs = c(25, 25)

if(Q_mat == 1){
  Q = matrix(
    c(
      0.2, 0.01,
      0.01, 0.2
    ),
    ncol = 2,
    byrow = TRUE
  )
}

if(Q_mat == 2){
  Q = matrix(
    c(
      0.6, 0.1,
      0.1, 0.6
    ),
    ncol = 2,
    byrow = TRUE
  )
}

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

###### Edge resampling #####

s_size = ceiling(choose(n,2)*0.05)
statistic_5 = numeric(M)

for (i in 1:M) {
  G = sample_IRG(p_0)
  sample.index = sample_index(n, s.size = s_size, replace = TRUE)
  statistic_5[i] = generate.one.gKSS.sampled(G, C, p_0,sample.index, g.kernel,level)$stats.value
}

critical_value_lower_5 = quantile(statistic_5,alpha/2)
critical_value_upper_5 = quantile(statistic_5,(1-alpha/2))
median(statistic_5)

##################

s_size = ceiling(choose(n,2)*0.25)
statistic_25 = numeric(M)

for (i in 1:M) {
  G = sample_IRG(p_0)
  sample.index = sample_index(n, s.size = s_size, replace = TRUE)
  statistic_25[i] = generate.one.gKSS.sampled(G, C, p_0,sample.index, g.kernel,level)$stats.value
}

critical_value_lower_25 = quantile(statistic_25,alpha/2)
critical_value_upper_25 = quantile(statistic_25,(1-alpha/2))
median(statistic_25)

##################

s_size = ceiling(choose(n,2)*0.75)
statistic_75 = numeric(M)

for (i in 1:M) {
  G = sample_IRG(p_0)
  sample.index = sample_index(n, s.size = s_size, replace = TRUE)
  statistic_75[i] = generate.one.gKSS.sampled(G, C, p_0,sample.index, g.kernel,level)$stats.value
}

critical_value_lower_75 = quantile(statistic_75,alpha/2)
critical_value_upper_75 = quantile(statistic_75,(1-alpha/2))
median(statistic_75)

###################### Testing part #########################

Q_12 = seq(0, 1, 0.05) 
l= length(Q_12)

power_gKSS = numeric(l)
power_gKSS_5 = numeric(l)
power_gKSS_25 = numeric(l)
power_gKSS_75 = numeric(l)
power_LRT = numeric(l)

for (i in 1:l) {
  
  Q_1 = Q
  Q_1[1,2] = Q_1[2,1] = Q_12[i]
  p_1 = matrix(Q_1[C,C], length(C), length(C))
  diag(p_1) = 0
  
  
  #######################################################


  decision_gKSS = numeric(m)
  decision_gKSS_5 = numeric(m)
  decision_gKSS_25 = numeric(m)
  decision_gKSS_75 = numeric(m)
  decision_LRT = numeric(m)
  
  for (j in 1:m) {
    
    test_G = sample_IRG(p_1)
    
    ########gKSS###########
    test_stat_gKSS = generate.one.gKSS(test_G, C, p_0, g.kernel,level)$stats.value
    decision_gKSS[j] = isTRUE(critical_value_lower >= test_stat_gKSS[j] ||test_stat_gKSS[j] >= critical_value_upper)
    #False means do not reject null
    
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
    
    s_size = ceiling(choose(n,2)*0.75)
    sample.index_75 = sample_index(n, s.size = s_size, replace = TRUE)
    test_stat_gKSS_75 = generate.one.gKSS.sampled(test_G, C, p_0, sample.index_75, g.kernel,level)$stats.value
    decision_gKSS_75[j] = isTRUE(critical_value_lower_75 >= test_stat_gKSS_75 ||test_stat_gKSS_75 >= critical_value_upper_75)
    #False means do not reject null
    
    #########LRT##########
    E = edge_counter_ermm(test_G, C)
    
    LR = ((1-Q[1,2])^(bs[1]*bs[2])*(Q[1,2]/(1-Q[1,2]))^(E[1,2]))/((1-Q_1[1,2])^(bs[1]*bs[2])*(Q_1[1,2]/(1-Q_1[1,2]))^(E[1,2]))
    
    test_stat_LRT = -2*log(LR)
    decision_LRT[j] = isTRUE(test_stat_LRT[j] > qchisq(1 - alpha, df=1) )
    #False means do not reject null
    
  }

  power_gKSS[i] = mean(decision_gKSS)
  power_gKSS_5[i] = mean(decision_gKSS_5)
  power_gKSS_25[i] = mean(decision_gKSS_25)
  power_gKSS_75[i] = mean(decision_gKSS_75)
  power_LRT[i] = mean(decision_LRT)
  print(i)
}

###### Constructing data frame ########

if(Q_mat == 1) Q_m = c(rep("01", length(Q_12)))
if(Q_mat == 2) Q_m = c(rep("02", length(Q_12)))

if(group_size == 1) n_vec = c(rep("A", length(Q_12)))
if(group_size == 2) n_vec = c(rep("B", length(Q_12)))

new_names = c("Q_12", "power","n_vec", "Q_m", "power_l", "power_u", "Test")


{power.LRT <- data.frame(Q_12, power_LRT, n_vec, Q_m)
  
  pow = power_LRT
  
  power.LRT$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.LRT$power_l[power.LRT$power_l < 0] <- 0
  power.LRT$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.LRT$power_u[power.LRT$power_u > 1] <- 1
  
  power.LRT$Test <- "LR"
  colnames(power.LRT) =  new_names
}

{power.gKSS <- data.frame(Q_12, power_gKSS, n_vec, Q_m)
  
  pow = power_gKSS
  
  power.gKSS$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.gKSS$power_l[power.gKSS$power_l < 0] <- 0
  power.gKSS$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.gKSS$power_u[power.gKSS$power_u > 1] <- 1
  
  power.gKSS$Test <- "no edge_resamp"
  colnames(power.gKSS) =  new_names
}

{power.gKSS_5 <- data.frame(Q_12, power_gKSS_5, n_vec, Q_m)
  
  pow = power_gKSS_5
  
  power.gKSS_5$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.gKSS_5$power_l[power.gKSS_5$power_l < 0] <- 0
  power.gKSS_5$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.gKSS_5$power_u[power.gKSS_5$power_u > 1] <- 1
  
  power.gKSS_5$Test <- "5%"
  colnames(power.gKSS_5) =  new_names
}

{power.gKSS_25 <- data.frame(Q_12, power_gKSS_25, n_vec, Q_m)
  
  pow = power_gKSS_25
  
  power.gKSS_25$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.gKSS_25$power_l[power.gKSS_25$power_l < 0] <- 0
  power.gKSS_25$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.gKSS_25$power_u[power.gKSS_25$power_u > 1] <- 1
  
  power.gKSS_25$Test <- "25%"
  colnames(power.gKSS_25) =  new_names
}

{power.gKSS_75 <- data.frame(Q_12, power_gKSS_75, n_vec, Q_m)
  
  pow = power_gKSS_75
  
  power.gKSS_75$power_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.gKSS_75$power_l[power.gKSS_75$power_l < 0] <- 0
  power.gKSS_75$power_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.gKSS_75$power_u[power.gKSS_75$power_u > 1] <- 1
  
  power.gKSS_75$Test <- "75%"
  colnames(power.gKSS_75) =  new_names
}

power.df = rbind(power.LRT, power.gKSS, power.gKSS_5, power.gKSS_25, power.gKSS_75)
assign(paste0("power_", Q_mat,group_size), power.df)

################## Creating data frame for grid plot ####################

power = rbind(power_11, power_12, power_21, power_22)

################### Top and Right facet names ###############################

group_size_names <- c(
  `A` = "Unbalanced_split",
  `B` = "Balanced_split"
)

Q_m_names <- c(
  `01` = "Q['01']",
  `02` = "Q['02']"
)

############## Final plot ##########

legend_order = c("5%", "25%", "75%", "no edge_resamp", "LR")
my_colours = hcl.colors(length(legend_order), "Dark 3")

power %>% ggplot(aes(x = Q_12, y = power)) + geom_point(aes(shape = Test, colour = Test)) + geom_line(aes(color = Test))+
  geom_ribbon(aes(ymin = power_l, ymax = power_u, x = Q_12, fill = Test), alpha = 0.2) +
  geom_abline(intercept = alpha, slope = 0, col = "gold") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="top", strip.background =element_rect(fill="#e4e7eb"), strip.text.x = element_text(size = 12)) +
  ylim(0, 1) + labs(y = "Proportion rejected", x = expression("Q"['12'])) + scale_color_manual(values=my_colours, breaks=legend_order) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4), breaks=legend_order) + scale_fill_manual(values = my_colours, breaks=legend_order) +
  facet_grid(n_vec ~ Q_m, labeller = labeller(n_vec  = as_labeller(group_size_names,  label_parsed), Q_m = as_labeller(Q_m_names, label_parsed)),scales = "free_x")
