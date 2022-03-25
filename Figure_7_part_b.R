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
set.seed(2021)
train_x_noise <- train_x 
train_x_noise <- train_x_noise/255  + matrix(
  rnorm(dim(train_x)[1]*dim(train_x)[2], mean = 0, sd=0.1),
  dim(train_x)[1],dim(train_x)[2])
# train_y
# 0    1    2    3    4    5    6    7    8    9
# 5923 6742 5958 6131 5842 5421 5918 6265 5851 5949

# stratify sampling
set.seed(1234)
test_sample <- data.frame(index = c(1:length(train_y)), label = train_y) %>%
  filter(label %in% c(0,1,3,8))%>%
  group_by(label) %>%
  slice_sample(n=500)

# we can subsample based on the labels
sampled_index <- test_sample$index
train_x_subsample <- train_x_noise[sampled_index,]
train_y_subsample <- train_y[sampled_index]

train_X_kmeans <- train_x_subsample
train_y_kmeans <- train_y_subsample

cov_mat_X <- coop::covar(train_X_kmeans)
# eigen-decmop
eigen_cov_X <- eigen(cov_mat_X)
U <- eigen_cov_X$vectors
D <- diag(sqrt(1/(eigen_cov_X$values+0.01))) # inverse sq
whiten_train_X_kmeans <- train_X_kmeans%*%t(U)%*%D%*%(U)

set.seed(2021)

estimate_MED_zero_inflation <- function(X,tol=0){
  for (j in c(1:ncol(X))){
    X[,j] <- X[,j]-median(X[,j])}
  pre_filter <- X^2
  pre_filter <- pre_filter[abs(pre_filter)>=tol]
  sigma_hat <- sqrt(median(pre_filter)/qchisq(1/2,df=1))
  return(sigma_hat)
}

sig_hat <- estimate_MED_zero_inflation(whiten_train_X_kmeans)


####p

set.seed(1234)
# pick best seed?
seed_list <- sample(2021,size=30,replace = FALSE)
k_means_list <- lapply(seed_list, function(x)kmeans_estimation(whiten_train_X_kmeans,
                                                               k=4,
                                                               iter.max = 30,seed = x))
within_ss_list <- lapply(k_means_list, function(x)x$objective[[x$iter]])
best_seed <- seed_list[which.min(unlist(within_ss_list))]
current_kmeans_k_3 <- k_means_list[[which.min(unlist(within_ss_list))]]
final_cluster_k_3 <- current_kmeans_k_3$cluster[[current_kmeans_k_3$iter]]
table(final_cluster_k_3,train_y_kmeans)
table(final_cluster_k_3,train_y_kmeans)/colSums(table(final_cluster_k_3,train_y_kmeans))
mclust::adjustedRandIndex(final_cluster_k_3,train_y_kmeans)

final_centroids_positive <- (current_kmeans_k_3$centers[[current_kmeans_k_3$iter]])

png(paste0(plot_output_dir,'Figure_7_b.png'),
    width = 6,height=5,res=300,units='in')

par(mfrow=c(2,2), mai = c(0.5, 0.5, 0.5, 0.5))
image(t(matrix(t(final_centroids_positive[1,]),nrow = 28,byrow = F)),
      main="Cluster 1", axes = FALSE,cex.main=2)
image(t(matrix(t(final_centroids_positive[2,]),nrow = 28,byrow = F)),
      main="Cluster 2", axes = FALSE,cex.main=2)
image(t(matrix(t(final_centroids_positive[3,]),nrow = 28,byrow = F)),
      main="Cluster 3", axes = FALSE,cex.main=2)
image(t(matrix(t(final_centroids_positive[4,]),nrow = 28,byrow = F)),
      main="Cluster 4", axes = FALSE,cex.main=2)

dev.off()


selective_result_1_3 <- kmeans_inference(whiten_train_X_kmeans, k=4, 1, 3,
                                     sig = sig_hat, iter.max = 20,
                                     seed = best_seed)

selective_result_1_2 <- kmeans_inference(whiten_train_X_kmeans, k=4, 1, 2,
                                         sig = sig_hat, iter.max = 20,
                                         seed = best_seed)

selective_result_1_4 <- kmeans_inference(whiten_train_X_kmeans, k=4, 1, 4,
                                         sig = sig_hat, iter.max = 20,
                                         seed = best_seed)

selective_result_2_3 <- kmeans_inference(whiten_train_X_kmeans, k=4, 2, 3,
                                         sig = sig_hat, iter.max = 20,
                                         seed = best_seed)
selective_result_2_4 <- kmeans_inference(whiten_train_X_kmeans, k=4, 2, 4,
                                         sig = sig_hat, iter.max = 20,
                                         seed = best_seed)

selective_result_3_4 <- kmeans_inference(whiten_train_X_kmeans, k=4, 3, 4,
                                         sig = sig_hat, iter.max = 20,
                                         seed = best_seed)

summary(selective_result_1_2)
summary(selective_result_1_3)
summary(selective_result_1_4)
summary(selective_result_2_3)
summary(selective_result_2_4)
summary(selective_result_3_4)


