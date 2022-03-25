library(keras)
library(tidyverse)
library(KmeansInference)
library(umap)

plot_output_dir <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/plot_output/"

mnist_df <- dataset_mnist()

train_x <- mnist_df$train$x
train_y <- mnist_df$train$y
table(train_y)
train_x <- array(as.numeric(train_x), dim = c(dim(train_x)[[1]], 784))
# normalize
set.seed(2022)
train_x_noise <- train_x 
train_x_noise <- train_x_noise/255  + matrix(
  rnorm(dim(train_x)[1]*dim(train_x)[2], mean = 0, sd=0.1),
           dim(train_x)[1],dim(train_x)[2])
# train_y
# 0    1    2    3    4    5    6    7    8    9 
# 5923 6742 5958 6131 5842 5421 5918 6265 5851 5949 

### negative example
test_sample_neg <- data.frame(index = c(1:length(train_y)), label = train_y) %>%
  filter(label %in% c(0)) %>%
  group_by(label) %>%
  slice_sample(n=1500)

# we can subsample based on the labels
sampled_index_neg <- test_sample_neg$index
train_x_subsample_neg <- train_x_noise[sampled_index_neg,]
train_y_subsample_neg <- train_y[sampled_index_neg]
cov_mat_X_neg <- coop::covar(train_x_subsample_neg)
# eigen-decmop
eigen_cov_X_neg <- eigen(cov_mat_X_neg)
U_neg <- eigen_cov_X_neg$vectors
D_neg <- diag(1/sqrt(eigen_cov_X_neg$values+0.01))
# whitened data
whiten_train_X_kmeans_neg <- train_x_subsample_neg%*%t(U_neg)%*%(D_neg)%*%(U_neg)

# variance estimation
estimate_MED_zero_inflation <- function(X,tol=0){
  for (j in c(1:ncol(X))){
    X[,j] <- X[,j]-median(X[,j])}
  pre_filter <- X^2
  pre_filter <- pre_filter[abs(pre_filter)>=tol]
  sigma_hat <- sqrt(median(pre_filter)/qchisq(1/2,df=1))
  return(sigma_hat)
}

sig_hat_neg <- estimate_MED_zero_inflation(whiten_train_X_kmeans_neg)


set.seed(1234)
seed_list <- sample(2021,size=30,replace = FALSE)
k_means_list <- lapply(seed_list, function(x)kmeans_estimation(whiten_train_X_kmeans_neg,
                                                               k=6,iter.max = 50,seed = x))
within_ss_list <- lapply(k_means_list, function(x)x$objective[[x$iter]])
best_seed <- seed_list[which.min(unlist(within_ss_list))] 
current_kmeans_k_3_neg <- k_means_list[[which.min(unlist(within_ss_list))]]
final_cluster_k_3_neg <- current_kmeans_k_3_neg$cluster[[current_kmeans_k_3_neg$iter]]
table(final_cluster_k_3_neg,train_y_subsample_neg)
### estimate sigma

set.seed(2021)

final_centroids <- (current_kmeans_k_3_neg$centers[[current_kmeans_k_3_neg$iter]])

png(paste0(plot_output_dir,'Figure_7_a.png'),
    width = 9,height=6,res=600,units='in')

par(mfrow=c(2,3))
image((matrix(t(final_centroids[1,]),nrow = 28,byrow = T)),
      main="Cluster 1", axes = FALSE,cex.main=2)
image((matrix(t(final_centroids[2,]),nrow = 28,byrow = T)),
      main="Cluster 2", axes = FALSE,cex.main=2)
image((matrix(t(final_centroids[3,]),nrow = 28,byrow = T)),
      main="Cluster 3", axes = FALSE,cex.main=2)
image((matrix(t(final_centroids[4,]),nrow = 28,byrow = T)),
      main="Cluster 4", axes = FALSE,cex.main=2)
image((matrix(t(final_centroids[5,]),nrow = 28,byrow = T)),
      main="Cluster 5", axes = FALSE,cex.main=2)
image((matrix(t(final_centroids[6,]),nrow = 28,byrow = T)),
      main="Cluster 6", axes = FALSE,cex.main=2)

dev.off()

pval_1_2 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 1, 2,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)

pval_1_3 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 1, 3,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)

pval_1_4 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 1, 4,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)

pval_1_5 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 1, 5,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)

pval_1_6 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 1, 6,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)
pval_3_4 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 3, 4,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)
pval_3_5 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 3, 5,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)

pval_3_6 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 3, 6,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)
pval_2_5 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 2, 5,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)
pval_2_3 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 2, 3,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)
pval_2_4 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 2, 4,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)

pval_2_6 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 2, 6,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)

pval_4_6 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 4, 6,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)

pval_4_5 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 4, 5,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)

pval_5_6 <- kmeans_inference(whiten_train_X_kmeans_neg, k=6, 5, 6,
                             sig = sig_hat_neg,
                             iter.max = 50,
                             seed = best_seed)

summary(pval_1_2)
summary(pval_1_3)
summary(pval_1_4)
summary(pval_1_5)
summary(pval_1_6)
summary(pval_2_3)
summary(pval_2_4)
summary(pval_2_5)
summary(pval_2_6)
summary(pval_3_4)
summary(pval_3_5)
summary(pval_3_6)
summary(pval_4_5)
summary(pval_4_6)
summary(pval_5_6)
