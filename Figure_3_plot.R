library(ggplot2)
library(clusterpval)
library(KmeansInference)
library(latex2exp)
library(ggpubr)

plot_output_dir <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/plot_output/"

norm_vec <- function(x) sqrt(sum(x^2))

set.seed(2022)
n <- 30
cl <- c(rep(1, 10), rep(2, 10), rep(3, 10))
delta <- 5
q <- 2
mu <- rbind(c(delta/2,rep(0,q-1)), c(rep(0,q-1), sqrt(3)*delta/2), c(-delta/2,rep(0,q-1)) )
sig <- 1
X <- matrix(rnorm(n*q, sd=sig), n, q) + mu[cl, ]

k <- 3
cluster_1 <- 1
cluster_2 <- 2

clusters <- kmeans_estimation(X, k,sig,iter.max = 20,seed = 2021)$final_cluster
inference_result <- kmeans_inference(X, k, cluster_1, cluster_2, sig = sig,
                 iter.max = 20, seed = 2021)

p1 <- ggplot(data.frame(X), aes(x=X1, y=X2, col=as.factor(clusters))) +
  geom_point(cex=2) + xlab("Feature 1") + ylab("Feature 2") +
  theme_classic(base_size=18) + theme(legend.position="none") +
  xlim(c(-4.2, 5.5)) + ylim(c(-4, 5.5)) +
  scale_colour_manual(values=c("dodgerblue3", "rosybrown", "orange")) +
  ggtitle(TeX("Original data ($\\phi$=||x^T$\\nu$||=$4.37$)"))+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

print(p1)


png(paste0(plot_output_dir,'figure_3_a','.png'),
    width = 4.5,height = 4.5, res=300,units='in')
print(p1)
dev.off()


k1 <- 1
k2 <- 2
prop_k2 <- sum(clusters == k2)/(sum(clusters == k1) + sum(clusters == k2))
diff_means <- colMeans(X[clusters == k1, ]) - colMeans(X[clusters == k2, ])
phi <- sqrt(0)
stat <- sqrt(sum(diff_means^2))
Xphi <- X
Xphi[clusters == k1, ] <- t(t(X[clusters == k1, ]) + prop_k2*(phi - stat)*diff_means/norm_vec(diff_means))
Xphi[clusters == k2, ] <- t(t(X[clusters == k2, ]) + (prop_k2 - 1)*(phi - stat)*diff_means/norm_vec(diff_means))

p2 <- ggplot(data.frame(Xphi), 
             aes(x=X1, y=X2, col=as.factor(clusters))) +
  geom_point(cex=2) + xlab("Feature 1") + ylab("Feature 2") +
  theme_classic(base_size=18) + theme(legend.position="none") +
  xlim(c(-4.2, 5.5)) + ylim(c(-4, 5.5)) +
  scale_colour_manual(values=c("dodgerblue3", "rosybrown", "orange")) +
  ggtitle(TeX("Perturbed data ($\\phi$=$0$)"))+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))


png(paste0(plot_output_dir,'figure_3_b','.png'),
    width = 4.5,height = 4.5, res=300,units='in')
print(p2)
dev.off()

phi2 <- 6
Xphi2 <- X
Xphi2[clusters == k1, ] <- t(t(X[clusters == k1, ]) + prop_k2*(phi2 - stat)*diff_means/norm_vec(diff_means))
Xphi2[clusters == k2, ] <- t(t(X[clusters == k2, ]) + (prop_k2 - 1)*(phi2 - stat)*diff_means/norm_vec(diff_means))
clusters_Xphi_2 <- kmeans(Xphi2,centers = 3, nstart = 1, algorithm = "Lloyd")$cluster

p3 <- ggplot(data.frame(Xphi2), aes(x=X1, y=X2, col=as.factor(clusters))) +
  geom_point(cex=2) + xlab("Feature 1") + ylab("Feature 2") +
  theme_classic(base_size=18) + theme(legend.position="none") +
  xlim(c(-4.2, 5.5)) + ylim(c(-4, 5.5)) +
  scale_colour_manual(values=c("dodgerblue3", "rosybrown", "orange")) +
  ggtitle(TeX("Perturbed data ($\\phi$=$6$)"))+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))


png(paste0(plot_output_dir,'figure_3_c','.png'),
    width = 4.5,height = 4.5, res=300,units='in')
print(p3)
dev.off()



