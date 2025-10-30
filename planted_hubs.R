######## Methodology ###########

# re-direct to the desired directory
setwd("") #for windows working directory
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

Q_mat = 1
group_size = 2

##### kernel to use ########

g.kernel= CalculateWLKernel
#CalculateGraphletKernel 
level = 2  # parameter for the kernel setting

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

cores <- detectCores()/2
cl <- makeCluster(3)
registerDoParallel(cl)

# Parallel sampling of null distribution
statistic <- foreach(i = 1:M, .combine = c, .packages = c("igraph", "graphkernels"),
                     .export = c("symmetricize")) %dopar% {
                       g <- sample_IRG(p_0)
                       generate.one.gKSS(g, C, p_0, g.kernel, level)$stats.value
                     }

stopCluster(cl)

critical_value_lower = quantile(statistic,alpha/2)
critical_value_upper = quantile(statistic,(1-alpha/2))
median(statistic)

##########################################################################
######### Gooddness-of-fit testing for networks with planted hubs ########
##########################################################################


######### Over k #############

# Define the number of hubs
R = 3
k = seq(2, 15, 1)
l = length(k)

power_gKSS = numeric(l)
power_gKSS_e = numeric(l)
power_GLR = numeric(l)
power_ST = numeric(l)
mean_max_deg = numeric(l)
median_max_deg = numeric(l)
sd_max_deg = numeric(l)
total_mismatches = numeric(l)
avg_attempts = numeric(l)

for (i in 1:l) {
  
  decision_gKSS = numeric(m)
  decision_gKSS_e = numeric(m)
  decision_GLR = numeric(m)
  decision_ST = numeric(m)
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
    
    ### gKSS
    stat_gKSS = generate.one.gKSS(test_G, C, p_0, g.kernel, level)$stats.value
    decision_gKSS[j] = isTRUE(critical_value_lower >= stat_gKSS || stat_gKSS >= critical_value_upper)
    #False means do not reject null
    
    ### GLR
    est = MLE.est(test_G, C)
    Q_h = est$estimates
    E = est$E
    
    stat_GLR = -2*((E[1,1]*log(Q[1,1]) + (choose(bs[1],2) - E[1,1])*log(1 - Q[1,1]) +
                      E[1,2]*log(Q[1,2]) + (bs[1]*bs[2] - E[1,2])*log(1 - Q[1,2]) +
                      E[2,2]*log(Q[2,2]) + (choose(bs[2],2) - E[2,2])*log(1 - Q[2,2])) - (E[1,1]*log(Q_h[1,1]) + (choose(bs[1],2) - E[1,1])*log(1 - Q_h[1,1]) +
                                                                                            E[1,2]*log(Q_h[1,2]) + (bs[1]*bs[2] - E[1,2])*log(1 - Q_h[1,2]) +
                                                                                            E[2,2]*log(Q_h[2,2]) + (choose(bs[2],2) - E[2,2])*log(1 - Q_h[2,2])))
    
    decision_GLR[j] = isTRUE(stat_GLR > qchisq(1 - alpha, df = 3))
    
    #### Spectral test
    A = as_adjacency_matrix(test_G, type = "both")
    
    K0 = 2
    stat_ST = GoFTestBoot(A, K0, n.dim = K0, n.b = 50, use.boot = T)$test.stat.boot 
    decision_ST[j] = isTRUE(max(stat_ST) > qtw(alpha / 2, lower.tail = F))
    
    ########gKSS with estimated parameters###########
    
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
    
    stat_gKSS_e = generate.one.gKSS(test_G, C, p_h, g.kernel,level)$stats.value
    decision_gKSS_e[j] = isTRUE(critical_value_lower_e >= stat_gKSS_e ||stat_gKSS_e >= critical_value_upper_e)
    #False means do not reject null
  }
  
  power_gKSS[i] = mean(decision_gKSS)
  power_gKSS_e[i] = mean(decision_gKSS_e)
  power_GLR[i] = mean(decision_GLR)
  power_ST[i] = mean(decision_ST)
  mean_max_deg[i] = mean(max_deg) 
  median_max_deg[i] = median(max_deg)
  sd_max_deg[i]  = sd(max_deg)
  total_mismatches[i] = sum(mismatches)
  avg_attempts[i] = mean(attempts)
  print(i)
  
  save(list = c("power_gKSS","power_gKSS_e", "power_GLR", "power_ST", "mean_max_deg", "median_max_deg", "sd_max_deg", "total_mismatches", "avg_attempts"), 
       file = paste0(Q_mat,group_size, ".RData"))
}

################## Constructing a data frame ######################

if(Q_mat == 0) Q_m = c(rep("00", length(k)))
if(Q_mat == 1) Q_m = c(rep("01", length(k)))
if(Q_mat == 2) Q_m = c(rep("02", length(k)))

if(group_size == 1) n_vec = c(rep("A", length(k)))
if(group_size == 2) n_vec = c(rep("B", length(k)))


power.df <- data.frame(k, power_gKSS, power_gKSS_e, power_GLR, power_ST, Q_m, n_vec, mean_max_deg, median_max_deg, sd_max_deg, total_mismatches, avg_attempts)

{
  pow = power_gKSS
  
  power.df$power_gKSS_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.df$power_gKSS_l[power.df$power_gKSS_l<0] <- 0
  power.df$power_gKSS_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.df$power_gKSS_u[power.df$power_gKSS_u > 1] <- 1
  
  pow = power_gKSS_e
  
  power.df$power_gKSS_e_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.df$power_gKSS_e_l[power.df$power_gKSS_e_l<0] <- 0
  power.df$power_gKSS_e_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.df$power_gKSS_e_u[power.df$power_gKSS_e_u > 1] <- 1
  
  pow = power_GLR
  
  power.df$power_GLR_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.df$power_GLR_l[power.df$power_GLR_l < 0] <- 0
  power.df$power_GLR_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.df$power_GLR_u[power.df$power_GLR_u > 1] <- 1
  
  pow = power_ST
  
  power.df$power_ST_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.df$power_ST_l[power.df$power_ST_l < 0] <- 0
  power.df$power_ST_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.df$power_ST_u[power.df$power_ST_u > 1] <- 1
  
}

assign("power.df_k", power.df)

######### Over Rep #############

# Define the number of hubs
R = seq(0, 10, 1)
k = 4
l = length(R)

power_gKSS = numeric(l)
power_gKSS_e = numeric(l)
power_GLR = numeric(l)
power_ST = numeric(l)
mean_max_deg = numeric(l)
median_max_deg = numeric(l)
sd_max_deg = numeric(l)
total_mismatches = numeric(l)
avg_attempts = numeric(l)

for (i in 1:l) {
  
  decision_gKSS = numeric(m)
  decision_gKSS_e = numeric(m)
  decision_GLR = numeric(m)
  decision_ST = numeric(m)
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
    
    ### gKSS
    stat_gKSS = generate.one.gKSS(test_G, C, p_0, g.kernel, level)$stats.value
    decision_gKSS[j] = isTRUE(critical_value_lower >= stat_gKSS || stat_gKSS >= critical_value_upper)
    #False means do not reject null
    
    ### GLR
    est = MLE.est(test_G, C)
    Q_h = est$estimates
    E = est$E
    stat_GLR = -2*((E[1,1]*log(Q[1,1]) + (choose(bs[1],2) - E[1,1])*log(1 - Q[1,1]) +
                      E[1,2]*log(Q[1,2]) + (bs[1]*bs[2] - E[1,2])*log(1 - Q[1,2]) +
                      E[2,2]*log(Q[2,2]) + (choose(bs[2],2) - E[2,2])*log(1 - Q[2,2])) - (E[1,1]*log(Q_h[1,1]) + (choose(bs[1],2) - E[1,1])*log(1 - Q_h[1,1]) +
                                                                                            E[1,2]*log(Q_h[1,2]) + (bs[1]*bs[2] - E[1,2])*log(1 - Q_h[1,2]) +
                                                                                            E[2,2]*log(Q_h[2,2]) + (choose(bs[2],2) - E[2,2])*log(1 - Q_h[2,2])))
    
    decision_GLR[j] = isTRUE(stat_GLR > qchisq(1 - alpha, df = 3))
    
    #### Spectral test
    A = as_adjacency_matrix(test_G, type = "both")
    
    K0 = 2
    stat_ST = GoFTestBoot(A, K0, n.dim = K0, n.b = 50, use.boot = T)$test.stat.boot 
    decision_ST[j] = isTRUE(max(stat_ST) > qtw(alpha / 2, lower.tail = F))
    
    ########gKSS with estimated parameters###########

    p_h = matrix(Q_h[C,C], length(C), length(C))
    diag(p_h) = 0
    
    cores <- detectCores()-2
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
    
    stat_gKSS_e = generate.one.gKSS(test_G, C, p_h, g.kernel,level)$stats.value
    decision_gKSS_e[j] = isTRUE(critical_value_lower_e >= stat_gKSS_e ||stat_gKSS_e >= critical_value_upper_e)
    #False means do not reject null
  }
  
  power_gKSS[i] = mean(decision_gKSS)
  power_gKSS_e[i] = mean(decision_gKSS_e)
  power_GLR[i] = mean(decision_GLR)
  power_ST[i] = mean(decision_ST)
  mean_max_deg[i] = mean(max_deg) 
  median_max_deg[i] = median(max_deg)
  sd_max_deg[i]  = sd(max_deg)
  total_mismatches[i] = sum(mismatches)
  avg_attempts[i] = mean(attempts)
  print(i)
  
  save(list = c("power_gKSS","power_gKSS_e", "power_GLR", "power_ST", "mean_max_deg", "median_max_deg", "sd_max_deg", "total_mismatches", "avg_attempts"), 
       file = paste0(Q_mat,group_size, ".RData"))
}

################## Constructing a data frame ######################

if(Q_mat == 0) Q_m = c(rep("00", length(R)))
if(Q_mat == 1) Q_m = c(rep("01", length(R)))
if(Q_mat == 2) Q_m = c(rep("02", length(R)))

if(group_size == 1) n_vec = c(rep("A", length(R)))
if(group_size == 2) n_vec = c(rep("B", length(R)))


power.df <- data.frame(R, power_gKSS, power_gKSS_e, power_GLR, power_ST, Q_m, n_vec, mean_max_deg, median_max_deg, sd_max_deg, total_mismatches, avg_attempts)

{
  pow = power_gKSS
  
  power.df$power_gKSS_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.df$power_gKSS_l[power.df$power_gKSS_l<0] <- 0
  power.df$power_gKSS_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.df$power_gKSS_u[power.df$power_gKSS_u > 1] <- 1
  
  pow = power_gKSS_e
  
  power.df$power_gKSS_e_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.df$power_gKSS_e_l[power.df$power_gKSS_e_l<0] <- 0
  power.df$power_gKSS_e_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.df$power_gKSS_e_u[power.df$power_gKSS_e_u > 1] <- 1
  
  pow = power_GLR
  
  power.df$power_GLR_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.df$power_GLR_l[power.df$power_GLR_l < 0] <- 0
  power.df$power_GLR_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.df$power_GLR_u[power.df$power_GLR_u > 1] <- 1
  
  pow = power_ST
  
  power.df$power_ST_l = pow - 2*sqrt((pow*(1-pow))/m)
  power.df$power_ST_l[power.df$power_ST_l < 0] <- 0
  power.df$power_ST_u = pow + 2*sqrt((pow*(1-pow))/m)
  power.df$power_ST_u[power.df$power_ST_u > 1] <- 1
  
}

assign("power.df_R", power.df)

################## Final Plot ####################

my_colours = c("#EE6677","#4477AA","#228833","#DDAA33", "#AA3377")

p1 <- power.df_R %>% ggplot(aes(x = R, y = power_gKSS)) + geom_point(aes(x = R, y = power_gKSS, shape = "power_gKSS", color = "power_gKSS")) + geom_line(aes(x = R, y = power_gKSS, color = "power_gKSS"))+
  geom_ribbon(aes(ymin = power_gKSS_l, ymax = power_gKSS_u, x = R), alpha = 0.2, fill = my_colours[1]) + #geom_label_repel(aes(x = R, y = power_gKSS, label = mean_max_deg), family = "mono", size = 4, min.segment.length = 0, segment.size = 0.2, force_pull = 0, nudge_y = 0.4, direction = "x", hjust = 0.5, vjust= 0) + #, max.overlaps=15) + 
  geom_point(aes(x = R, y = power_gKSS_e, shape = "power_gKSS_e", color = "power_gKSS_e")) + geom_line(aes(x = R, y = power_gKSS_e, color = "power_gKSS_e"))+
  geom_ribbon(aes(ymin = power_gKSS_e_l, ymax = power_gKSS_e_u, x = R), alpha = 0.2, fill = my_colours[2])+
  geom_point(aes(x = R, y = power_GLR, shape = "power_GLR", color ="power_GLR")) +  geom_line(aes(x = R, y = power_GLR, color = "power_GLR")) +
  geom_ribbon(aes(ymin = power_GLR_l, ymax = power_GLR_u, x = R), alpha = 0.2, fill = my_colours[3]) +
  geom_point(aes(x = R, y = power_ST, shape = "power_ST", color ="power_ST")) +  geom_line(aes(x = R, y = power_ST, color = "power_ST")) +
  geom_ribbon(aes(ymin = power_ST_l, ymax = power_ST_u, x = R), alpha = 0.2, fill = my_colours[4]) +
  geom_abline(intercept = alpha, slope = 0, col = "gold") + 
  scale_color_manual(name = "Test", labels = c("IRG-gKSS","IRG-gKSS_e", "GLR","ST"), values = c(my_colours[1], my_colours[2], my_colours[3], my_colours[4])) + scale_shape(name = "Test", labels = c( "IRG-gKSS","IRG-gKSS_e", "GLR", "ST")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="right", strip.text.y = element_text(size = 12), strip.text.x = element_text(size = 12), strip.background =element_rect(fill="#e4e7eb"), axis.title.y=element_blank()) +
  ylim(0, 1) + scale_x_continuous(breaks = seq(0, 10, by = 2)) + coord_cartesian(clip = 'off') +
  labs(x = "R") + geom_label(aes(x = R, y = I(1), label = ceiling(mean_max_deg)), angle = 0, size = 3, vjust = 0.2) #+ geom_label(aes(x = R, y = I(1), label = round(sd_max_deg,2)), angle = 0, colour = my_colours[4], vjust = 0.5)

p2 <- power.df_k %>% ggplot(aes(x = k, y = power_gKSS)) + geom_point(aes(x = k, y = power_gKSS, shape = "power_gKSS", color = "power_gKSS")) + geom_line(aes(x = k, y = power_gKSS, color = "power_gKSS"))+
  geom_ribbon(aes(ymin = power_gKSS_l, ymax = power_gKSS_u, x = k), alpha = 0.2, fill = my_colours[1]) + #geom_label_repel(aes(x = k, y = power_gKSS, label = mean_max_deg), family = "mono", size = 4, min.segment.length = 0, segment.size = 0.2, force_pull = 0, nudge_y = 0.4, direction = "x", hjust = 0.5, vjust= 0) + #, max.overlaps=15) + 
  geom_point(aes(x = k, y = power_gKSS_e, shape = "power_gKSS_e", color = "power_gKSS_e")) + geom_line(aes(x = k, y = power_gKSS_e, color = "power_gKSS_e"))+
  geom_ribbon(aes(ymin = power_gKSS_e_l, ymax = power_gKSS_e_u, x = k), alpha = 0.2, fill = my_colours[1]) +
  geom_point(aes(x = k, y = power_GLR, shape = "power_GLR", color ="power_GLR")) +  geom_line(aes(x = k, y = power_GLR, color = "power_GLR")) +
  geom_ribbon(aes(ymin = power_GLR_l, ymax = power_GLR_u, x = k), alpha = 0.2, fill = my_colours[3]) +
  geom_point(aes(x = k, y = power_ST, shape = "power_ST", color ="power_ST")) +  geom_line(aes(x = k, y = power_ST, color = "power_ST")) +
  geom_ribbon(aes(ymin = power_ST_l, ymax = power_ST_u, x = k), alpha = 0.2, fill = my_colours[4]) +
  geom_abline(intercept = alpha, slope = 0, col = "gold") + 
  scale_color_manual(name = "Test", labels = c("IRG-gKSS","IRG-gKSS_e", "GLR", "ST"), values = c(my_colours[1],my_colours[2], my_colours[3], my_colours[4])) + scale_shape(name = "Test", labels = c( "IRG-gKSS","IRG-gKSS_e", "GLR", "ST")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="right", strip.text.y = element_text(size = 12), strip.text.x = element_text(size = 12), strip.background =element_rect(fill="#e4e7eb"), axis.title.y=element_blank()) +
  ylim(0, 1) + scale_x_continuous(breaks = seq(0, 15, by = 2)) + coord_cartesian(clip = 'off') +
  labs( x = "k") + geom_label(aes(x = k, y = I(1), label = ceiling(mean_max_deg)), angle = 0, size = 3, vjust = 0.2) #+ geom_label(aes(x = k, y = I(1), label = round(sd_max_deg,2)), angle = 0, colour = my_colours[4], vjust = 0.5)

multi <- ggarrange(p1,p2, #plots that are going to be included in this multipanel figure
                   #labels = c("A", "B", "C","D"), #labels given each panel 
                   ncol = 2, nrow = 1, #adjust plot space 
                   common.legend = T, widths = c(0.8, 1)) 
annotate_figure(multi, 
                left = "Proportion rejected")


