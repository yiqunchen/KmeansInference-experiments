library(palmerpenguins)
library(ggplot2)
library(patchwork)
library(KmeansInference)

options(ggplot2.discrete.colour=list(RColorBrewer::brewer.pal(6, "Accent")))


input_dir  <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/output/"
plot_output_dir <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/plot_output/"

# Subset to female penguins in years & bill depth; flipper length variables
penguins <- penguins[complete.cases(penguins), ]
dat <- penguins[penguins$sex == "female", c(1, 4, 5)]

X <- dat[, -c(1)]
X <- as.matrix(X)

estimate_MED <- function(X){
  for (j in c(1:ncol(X))){
    X[,j] <- X[,j]-median(X[,j])}
    sigma_hat <- sqrt(median(X^2)/qchisq(1/2,df=1))
    return(sigma_hat)
  }

sig_hat_MED <- estimate_MED(X)




# Cluster and visualize data; K=4
set.seed(1234)
# pick best rnadom seed out of 30 initial seeds 
seed_list <- sample(2021,size=30,replace = FALSE)
k_means_list <- lapply(seed_list, function(x)kmeans_estimation(X,
                                               k=4,iter.max = 10,seed = x))
within_ss_list <- lapply(k_means_list, function(x)x$objective[[x$iter]])
best_seed <- seed_list[which.min(unlist(within_ss_list))] 
current_kmeans_k_5_peng <- k_means_list[[which.min(unlist(within_ss_list))]]
final_cluster_k_5_peng <- current_kmeans_k_5_peng$cluster[[current_kmeans_k_5_peng$iter]]
table(final_cluster_k_5_peng,dat$species)


p2 <- ggplot(data = penguins[penguins$sex == "female", ]) +
  geom_point(aes(x=flipper_length_mm , y = bill_depth_mm,
                 fill = as.factor(final_cluster_k_5_peng),
                 shape=as.factor(species)), size = 4, colour="black") +
  scale_fill_brewer(name="Clusters",palette="Dark2",
                      guide=guide_legend(ncol=2,
                                        override.aes=list(shape=21))) +
  scale_shape_manual(name="Species", values=c(21, 24, 22),
                     guide=guide_legend(ncol=2,override.aes=list(fill="black"))) +
  ylab("Bill depth (mm)") + xlab("Flipper length (mm)") +
  coord_flip() +
  theme_bw(base_size=12) +
  theme(legend.position="none",  plot.margin=unit(c(0,0,0,0),"mm"))+
  coord_fixed(ratio=4)

png(paste0(plot_output_dir,'Figure_6.png'),
    width = 3.15,height=2.2,res=400,units='in')
p2
dev.off()


# naive and selective p-vals
test_1_2 <- kmeans_inference(X, k=4, 1, 2,
                             sig = sig_hat_MED,
                             iter.max = 10,
                             seed = best_seed)
test_1_3 <- kmeans_inference(X, k=4, 1, 3,
                             sig = sig_hat_MED,
                             iter.max = 10,
                             seed = best_seed)
test_1_4 <- kmeans_inference(X, k=4, 1, 4,
                             sig = sig_hat_MED,
                             iter.max = 10,
                             seed = best_seed)
test_2_3 <- kmeans_inference(X, k=4, 2, 3,
                             sig = sig_hat_MED,
                             iter.max = 10,
                             seed = best_seed)
test_2_4 <- kmeans_inference(X, k=4, 2, 4,
                             sig = sig_hat_MED,
                             iter.max = 10,
                             seed = best_seed)
test_3_4 <- kmeans_inference(X, k=4, 4, 3,
                             sig = sig_hat_MED,
                             iter.max = 10,
                             seed = best_seed)

summary(test_1_2)
summary(test_1_3)
summary(test_1_4)
summary(test_2_3)
summary(test_2_4)
summary(test_3_4)














