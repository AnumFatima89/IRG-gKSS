setwd("...") #for windows working directory
source("source-codes.R")

########## Paramter setting #########
M = 200
alpha = 0.05

n = 30
bs = c(15, 15)

Q = matrix(
  c(
    0.29, 0.01,
    0.01, 0.22
  ),
  ncol = 2,
  byrow = TRUE
)

C = c(rep(1,bs[1]), rep(2,bs[2]))

p = matrix(Q[C,C], length(C), length(C))
diag(p) = 0

p_0 = p

###################### Kernel setting ##########################

g.kernel = CalculateWLKernel

######## GOF test null sets ##########

statistic_2 = numeric(M)
for (i in 1:M) {
  g <- sample_IRG(p_0)
  statistic_2[i] = generate.one.gKSS(g, C, p_0, g.kernel, level=2)$stats.value
}

critical_value_lower_2 = quantile(statistic_2,alpha/2)
critical_value_upper_2 = quantile(statistic_2,(1-alpha/2))
median(statistic_2)

statistic_3 = numeric(M)
for (i in 1:M) {
  g <- sample_IRG(p_0)
  statistic_3[i] = generate.one.gKSS(g, C, p_0, g.kernel, level=3)$stats.value
}

critical_value_lower_3 = quantile(statistic_3,alpha/2)
critical_value_upper_3 = quantile(statistic_3,(1-alpha/2))
median(statistic_3)

########## Simulated networks ##########

R = 2
k = 3

G = sample_IRG(p)
G_h = planted_hubs(G, R, k, C)$Graph

###### For h=2 #####

TN_test_2 = generate.one.gKSS(G, C, p, g.kernel,level=2)
TN_test_stat_2 = TN_test_2$stats.value
TN_decision_2 = isTRUE(critical_value_lower_2 >= TN_test_stat_2 || TN_test_stat_2 >= critical_value_upper_2)
#False means do not reject null

test_2 = generate.one.gKSS(G_h, C, p, g.kernel,level=2)
test_stat_2 = test_2$stats.value
decision_2 = isTRUE(critical_value_lower_2 >= test_stat_2 || test_stat_2 >= critical_value_upper_2)
#False means do not reject null

###### For h=3 ######

TN_test_3 = generate.one.gKSS(G, C, p, g.kernel,level=3)
TN_test_stat_3 = TN_test_3$stats.value
TN_decision_3 = isTRUE(critical_value_lower_3 >= TN_test_stat_3 || TN_test_stat_3 >= critical_value_upper_3)
#False means do not reject null

test_3 = generate.one.gKSS(G_h, C, p, g.kernel,level=3)
test_stat_3 = test_3$stats.value
decision_3 = isTRUE(critical_value_lower_3 >= test_stat_3 || test_stat_3 >= critical_value_upper_3)
#False means do not reject null

c(TN_decision_2, TN_decision_3, decision_2, decision_3)

######### h(s,s') matrix calculation ###########

Kernel_G_2 = TN_test_2$J.kernel
Kernel_Gh_2 = test_2$J.kernel

Kernel_G_3 = TN_test_3$J.kernel
Kernel_Gh_3 = test_3$J.kernel

################# Heat map plot ########################

Ker_mat = Kernel_G_2

r = nrow(Ker_mat)
rownames(Ker_mat) <- paste("r", 1:r)
rownames(Ker_mat) <- paste(1:r)

library(reshape)
df_G_2 <- melt(Ker_mat)
df_G_2$B <- "G"
df_G_2$h <- 2
colnames(df_G_2) <- c("x", "y", "value","B", "h")

####

Ker_mat = Kernel_Gh_2

r = nrow(Ker_mat)
rownames(Ker_mat) <- paste("r", 1:r)
rownames(Ker_mat) <- paste(1:r)

df_Gh_2 <- melt(Ker_mat)
df_Gh_2$B <- "G_hubs"
df_Gh_2$h <- 2
colnames(df_Gh_2) <- c("x", "y", "value","B", "h")

####

Ker_mat = Kernel_G_3

r = nrow(Ker_mat)
rownames(Ker_mat) <- paste("r", 1:r)
rownames(Ker_mat) <- paste(1:r)

df_G_3 <- melt(Ker_mat)
df_G_3$B <- "G"
df_G_3$h <- 3
colnames(df_G_3) <- c("x", "y", "value","B", "h")

####

Ker_mat = Kernel_Gh_3

r = nrow(Ker_mat)
rownames(Ker_mat) <- paste("r", 1:r)
rownames(Ker_mat) <- paste(1:r)

df_Gh_3 <- melt(Ker_mat)
df_Gh_3$B <- "G_hubs"
df_Gh_3$h <- 3
colnames(df_Gh_3) <- c("x", "y", "value","B", "h")

######### Combined data frames ##############

df2 = rbind(df_G_2, df_Gh_2)
df3 = rbind(df_G_3, df_Gh_3)
df = rbind(df2, df3)

####### Heatmap plots #######

ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile() + scale_fill_gradient(low = "#f1f3f5", high = "black") + 
  coord_fixed() + 
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     strip.text.x = element_text(size = 12), 
                     strip.background =element_rect(fill="#f1f3f5")
  ) +
  facet_grid(h~B) + labs(y = "", x = "")
