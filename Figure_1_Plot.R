library(rdist)
library(ggplot2)
library(KmeansInference)
library(dplyr)
library(latex2exp)
library(ggpubr)

set.seed(54321)
n <- 100
cl <- c(rep(1, 30), rep(2, 30), rep(3, 40))
mu <- rbind(c(0, 0), c(0, 0), c(sqrt(0), 0))
sig <- 1
p <- 2
k <- 3

input_dir  <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/output/"
plot_output_dir <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/plot_output/"


X_sim <- matrix(rnorm(n*p, sd=sig), n, p) + mu[cl, ]
current_kmeans <- kmeans_estimation(X_sim, k,sig,iter.max = 10,seed = 2021)
df_panel_a <- data.frame(X_sim) %>% 
  mutate(clusters = current_kmeans$final_cluster)

centroids <- aggregate(cbind(X1, X2)~clusters, df_panel_a, mean)

p1_a <- ggplot(df_panel_a, aes(x=X1, y=X2, col=as.factor(clusters))) + 
  geom_point(cex=4) + xlab("Feature 1") + ylab("Feature 2") + 
  theme_classic(base_size=18) + theme(legend.position="none") + 
  geom_point(aes(x=X1, y=X2, col=as.factor(clusters)), shape=17, cex=6,data=centroids)+
  #xlim(c(-4.2, 7)) + ylim(c(-4, 4.5)) + 
  scale_colour_manual(values=c("dodgerblue3", "rosybrown", "orange")) + 
  ggtitle("")


png(paste0(plot_output_dir,'figure_1_a','.png'),
    width = 5,height = 5, res=300,units='in')
print(p1_a)
dev.off()


load(paste0(input_dir,"Figure_1_Data.RData"))

p1_b <- data.frame(naive_p_val = naive_vec) %>%
  arrange(naive_p_val) %>%
  mutate(theoretical = c(1:nrow(.))/nrow(.)) %>%
  ggplot() +
  geom_point(aes(y = naive_p_val, x = (theoretical)),size=2.5) +
  ylab('Naive p-value quantiles')+
  xlab('Uniform(0,1) quantiles')+
  #ggtitle(TeX('Wald test p-values'))+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5,size=18),
        #legend.position="bottom",
        #legend.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18,hjust = 0),
        axis.title=element_text(size=18))

png(paste0(plot_output_dir,'figure_1_b','.png'),
    width = 5,height = 5, res=500,units='in')
print(p1_b)
dev.off()

p1_c <- data.frame(selective_p_val = selective_vec) %>%
  arrange(selective_p_val) %>%
  mutate(theoretical = c(1:nrow(.))/nrow(.)) %>%
  ggplot() +
  geom_point(aes(y = selective_p_val, x = (theoretical)),size=2.5) +
  ylab('Selective p-value quantiles')+
  xlab('Uniform(0,1) quantiles')+
  #ggtitle(TeX('Selective test p-values'))+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5,size=18),
        #legend.position="bottom",
        #legend.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18,hjust = 0),
        axis.title=element_text(size=18))


png(paste0(plot_output_dir,'figure_1_c','.png'),
    width = 5,height = 5, res=500,units='in')
print(p1_c)
dev.off()


