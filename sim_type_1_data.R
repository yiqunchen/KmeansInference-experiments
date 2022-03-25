library(rdist)
library(ggplot2)
library(KmeansInference)
library(pracma)

norm_vec <- function(x) {
  ell_2_norm <- sqrt(sum(x^2))
  return(ell_2_norm)
}

delta <- 0
sig_seq <- c(0.25, 0.5, 1)
p_seq <- c(2,10,50,100)
k <- 3#as.numeric(args[4])
output_dir <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/output/"
#cat("delta",delta,"sig",sig,"p",p,"\n")
set.seed(1234)

n <- 150
cl <- c(rep(1, 50), rep(2, 50), rep(3, 50))


estimate_MED <- function(X){
 for (j in c(1:ncol(X))){
   X[,j] <- X[,j]-median(X[,j])
}
  return(median(X^2)/qchisq(1/2,df=1))
}


for (p in p_seq){
  for (sig in sig_seq){
    if(delta==0){
      n_sim <- 3000
    }else{
      n_sim <- 2000*25
    }
    mu <- rbind(c(-delta/2, rep(0, p-1)), 
    c(rep(0, p-1), sqrt(3)*delta/2), 
    c(delta/2, rep(0, p-1))) 

    naive_vec <- rep(NA, times=n_sim)
    selective_vec <- rep(NA, times=n_sim)
    selective_vec_sample <- rep(NA, times=n_sim)
    selective_vec_MED <- rep(NA, times=n_sim)
    selective_object_vec <- vector("list", length =n_sim)
    true_mu_diff <- rep(NA, times=n_sim)
    adjust_rand_index <- rep(NA, times=n_sim)
    nu_norm <- rep(NA, times=n_sim)
    # run sims 
    for (i in c(1:n_sim)){
      set.seed(i)
      cat("Iteration", i, "\n")
      # draw a sample from X
      X_sim <- matrix(rnorm(n*p, sd=sig), n, p) + mu[cl, ]
      current_kmeans <- kmeans_estimation(X_sim, k,iter.max = 30,seed = i)
      # sigma hat estimation
      # if k-means failed -- try again
      if(length(unique(current_kmeans$final_cluster))<k){ next }

      sigma_hat_sample <- sqrt(sum(scale(X_sim, scale = FALSE)^2)/
                                        (nrow(X_sim) * ncol(X_sim)-ncol(X_sim)))
      sigma_hat_MED <- sqrt(estimate_MED(X_sim))

      random_pair <- sample(c(1:k),2,replace=FALSE)
      cluster_1 <- random_pair[1]
      cluster_2 <- random_pair[2]
      # naive
      naive_pval <- multivariate_Z_test(X_sim, current_kmeans$final_cluster,cluster_1,cluster_2,sig = sig)$pval
      selective_result <- kmeans_inference(X_sim, k, cluster_1, cluster_2, sig = sig, iter.max = 30, seed = i)
      selective_pval <- selective_result$pval
      selective_p_hat_MED <- kmeans_inference(X_sim, k, cluster_1, cluster_2, sig = sigma_hat_MED, iter.max = 30, seed = i)
      selective_p_hat_SAMPLE <- kmeans_inference(X_sim, k, cluster_1, cluster_2, sig = sigma_hat_sample, iter.max = 30, seed = i)
      naive_vec[i] <- naive_pval
      selective_vec[i] <- selective_pval
      selective_vec_sample[i] <- selective_p_hat_SAMPLE$pval
      selective_vec_MED[i] <- selective_p_hat_MED$pval 
      selective_object_vec[[i]] <- selective_result
      cat(selective_pval,
        selective_p_hat_SAMPLE$pval,
        selective_p_hat_MED$pval,"\n")
      true_mean_X <- mu[cl, ]
      nu_norm[i] <- sqrt(1/sum(current_kmeans$final_cluster == cluster_1)+
        1/sum(current_kmeans$final_cluster == cluster_2))
      true_mu_diff[i] <- norm_vec(colMeans(true_mean_X[current_kmeans$final_cluster == cluster_1, , drop=F]) -
        colMeans(true_mean_X[current_kmeans$final_cluster == cluster_2, , drop=F]))/sig
      adjust_rand_index[i] <- mclust::adjustedRandIndex(current_kmeans$final_cluster,cl)

    }

    result_df <- data.frame(naive_vec=naive_vec,nu_norm=nu_norm,selective_vec=selective_vec,
    selective_vec_sample=selective_vec_sample,selective_vec_MED=selective_vec_MED,
    true_mu_diff=true_mu_diff,adjust_rand_index=adjust_rand_index,k=k,delta=delta,sigma=sig,p=p)
    save(result_df, naive_vec,nu_norm,selective_vec,selective_vec_sample,selective_vec_MED,true_mu_diff,adjust_rand_index,
         file=paste0(output_dir,"validation_sim","_k_",k,"_delta_",delta,"_sigma_", sig, "_p_",p,".RData"))


  }
}

