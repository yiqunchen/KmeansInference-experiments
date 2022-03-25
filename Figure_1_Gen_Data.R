library(rdist)
library(ggplot2)
library(KmeansInference)

output_dir <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/output/"
n <- 100
cl <- c(rep(1, 30), rep(2, 30), rep(3, 40))
mu <- rbind(c(0, 0), c(0, 0), c(sqrt(0), 0))
sig <- 1
p <- 2
k <- 3

n_sim <- 2000
naive_vec <- rep(NA, times=n_sim)
selective_vec <- rep(NA, times=n_sim)

for (i in c(1:n_sim)){
  set.seed(i) # set seed
  if(i%%200==0){
    cat("Iteration", i, "\n")
  }
  X_sim <- matrix(rnorm(n*p, sd=sig), n, p) + mu[cl, ]
  current_kmeans <- kmeans_estimation(X_sim, k,iter.max = 20,seed = i)
  # naive
  random_pair <- sample(c(1:k),2,replace=FALSE)
  cluster_1 <- random_pair[1]
  cluster_2 <- random_pair[2]
  naive_pval <- multivariate_Z_test(X_sim,
                                    current_kmeans$final_cluster,
                                    cluster_1,cluster_2,sig = sig)$pval
  selective_result <- kmeans_inference(X_sim, k, cluster_1, cluster_2, sig = sig,
                                       iter.max = 20, seed = i)
  selective_pval <- selective_result$pval
  naive_vec[i] <- naive_pval
  selective_vec[i] <- selective_pval
}

save(naive_vec,selective_vec,
     file=paste0(output_dir,"Figure_1_Data.RData"))

