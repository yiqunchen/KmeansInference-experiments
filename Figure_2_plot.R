library(rdist)
library(ggplot2)
library(KmeansInference)
library(dplyr)
library(latex2exp)
library(ggpubr)
library(class)

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
current_kmeans <- kmeans_estimation(X_sim, k,sig,iter.max = 20,seed = 2021)
df_panel_a <- data.frame(X_sim) #%>% mutate(clusters = current_kmeans$final_cluster)

#centroids <- aggregate(cbind(X1, X2)~clusters, df_panel_a, mean)


dat.train <- X_sim[1:50, ]
dat.test <- X_sim[51:100, ]

blankdata <- ggplot(df_panel_a) + geom_point(aes(X1, X2), cex=3) +
  xlab("") + ylab("")  +
  theme_bw() + theme(axis.ticks=element_blank(), axis.text=element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.border=element_blank()) 

kmeans_train <- kmeans_estimation(dat.train, k,sig,iter.max = 20,seed = 2021)
dat.train <- data.frame(dat.train) %>% mutate(clusters = kmeans_train$final_cluster)
#hc <- hclust(dist(X[1:50, ])^2, method="average")
#dat.train$clusters <- as.factor(cutree(hc, 3))

plot_train_data <- ggplot(dat.train) +
  geom_point(aes(X1, X2, col=as.factor(clusters)), cex=4, alpha=1, shape="square") +
  xlab("") + ylab("")  + scale_colour_manual(name="Clusters",
                               values=c("dodgerblue3", "rosybrown", "orange")) +
  theme_bw(base_size=18) +theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position="none") +
  xlim(c(-2.8, 2.8)) +  ylim(c(-2.6, 2.8)) +
  #ggtitle("Training set") + 
  xlab("Feature 1") + ylab("Feature 2")+
  theme(plot.title = element_text(hjust = 0.5,size=18),
        #legend.position="bottom",
        #legend.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18,hjust = 0),
        axis.title=element_text(size=18))


png(paste0(plot_output_dir,'figure_2_a','.png'),
    width = 5,height = 5, res=300,units='in')
print(plot_train_data)
dev.off()


dat.test <- data.frame(dat.test) %>% 
  mutate(clusters = knn(X_sim[1:50, ], X_sim[51:100, ], dat.train$clusters, k=3))

plot_test_data <- ggplot(dat.test) +
  geom_point(aes(X1, X2, col=as.factor(clusters)), cex=4, alpha=1, shape="triangle") +
  xlab("") + ylab("")  + scale_colour_manual(name="Clusters",
                                             values=c("dodgerblue3", "rosybrown", "orange")) +
  theme_bw(base_size=18) +theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position="none") +
  xlim(c(-2.8, 2.8)) +  ylim(c(-2.6, 2.8)) +
  #ggtitle("Test set") + 
  xlab("Feature 1") + ylab("Feature 2")+
  theme(plot.title = element_text(hjust = 0.5,size=18),
        #legend.position="bottom",
        #legend.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18,hjust = 0),
        axis.title=element_text(size=18))

png(paste0(plot_output_dir,'figure_2_b','.png'),
    width = 5,height = 5, res=300,units='in')
print(plot_test_data)
dev.off()

load(paste0(input_dir,"Figure_2_Data.RData"))

p2_c <- data.frame(naive_p_val = sample_split_pval_vec) %>%
  arrange(naive_p_val) %>%
  mutate(theoretical = c(1:nrow(.))/nrow(.)) %>%
  ggplot() +
  geom_point(aes(y = naive_p_val, x = (theoretical)),size=2.5) +
  ylab('Sample-split p-value quantiles')+
  xlab('Uniform(0,1) quantiles')+
  #ggtitle(TeX('Sample-split Wald test p-values'))+
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

png(paste0(plot_output_dir,'figure_2_c','.png'),
    width = 5,height = 5, res=500,units='in')
print(p2_c)
dev.off()

