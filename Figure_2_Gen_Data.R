library(rdist)
library(ggplot2)
library(KmeansInference)
library(dplyr)
library(latex2exp)
library(ggpubr)
library(class)

norm_vec <- function(x) sqrt(sum(x^2))

multivariate_Z_test <- function(X, cluster_vec, k1, k2, sig) {
  q <- ncol(X)
  diff_means <- colMeans(X[cluster_vec == k1, , drop=F]) -
    colMeans(X[cluster_vec == k2, , drop=F])
  stat <- norm_vec(diff_means)
  n1 <- sum(cluster_vec == k1)
  n2 <- sum(cluster_vec == k2)
  squared_norm_nu <- 1/n1 + 1/n2
  scale_factor <- squared_norm_nu*sig^2
  accurate_pchi <- pchisq(stat^2/scale_factor, df=q, log.p = TRUE,lower.tail=FALSE)
  pval <- exp(accurate_pchi) #1 - pchisq(stat^2/scale_factor, df=q)
  if(is.na(pval)){pval <- 1}
  return(list(stat=stat, pval=pval))
}

set.seed(54321)
n <- 100
cl <- c(rep(1, 30), rep(2, 30), rep(3, 40))
mu <- rbind(c(0, 0), c(0, 0), c(sqrt(0), 0))
sig <- 1
p <- 2
k <- 3

input_dir  <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/output/"
plot_output_dir <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/plot_output/"

### let's do the sample splitting simulation
n_sim <- 2000
sample_split_pval_vec <- rep(NA, times=n_sim)
#selective_vec <- rep(NA, times=n_sim)
#selective_object_vec <- vector("list", length =n_sim)

for (i in c(1:n_sim)){
  set.seed(i)
  if(i%%200==0){
    cat("Iteration", i, "\n")
  }
  X_sim <- matrix(rnorm(n*p, sd=sig), n, p) + mu[cl, ]
  X_train <- X_sim[1:50,]
  X_test <- X_sim[51:100,]
  current_kmeans <- kmeans_estimation(X_train, k, sig, iter.max = 20, seed = i)
  test_clusters <- knn(X_train, X_test, current_kmeans$final_cluster, k=3)
  # naive
  random_pair <- sample(c(1:k),2,replace=FALSE)
  cluster_1 <- random_pair[1]
  cluster_2 <- random_pair[2]
  sample_split_pval_vec[i] <- multivariate_Z_test(X_test, test_clusters,
                                    cluster_1,cluster_2,sig = sig)$pval
}

save(sample_split_pval_vec,file=paste0(input_dir,"Figure_2_Data.RData"))





