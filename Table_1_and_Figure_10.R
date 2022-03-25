library(DropletUtils)
library(scater)
library(ggfortify)
library(patchwork)
library(KmeansInference)
library(umap)

input_dir <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/output/"
plot_output_dir <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/plot_output/"
setwd(input_dir)

##### Pre-processing and splitting data ####
cd4.t <- read10xCounts("./raw/filtered_matrices_mex_memory/hg19")
cd19.b <- read10xCounts( "./raw/filtered_matrices_mex_bcell/hg19")
mono <- read10xCounts( "./raw/filtered_matrices_mex_mono/hg19")
nk <- read10xCounts( "./raw/filtered_matrices_mex_nk/hg19/")
naive.cytotoxic <-  read10xCounts( "./raw/filtered_matrices_mex_naive_cytotoxic/hg19/")

pure.t <- cd4.t
colData(pure.t)$CellType <- c(rep("Memory T cell", ncol(cd4.t)))

mix5 <- cbind(cd4.t, cd19.b, mono, nk, naive.cytotoxic)
colData(mix5)$CellType <- c(rep("Memory T cell", ncol(cd4.t)),
                            rep("B cell", ncol(cd19.b)),
                            rep("Monocyte", ncol(mono)),
                            rep("Nature Killer", ncol(nk)),
                            rep("Naive Cytotoxic", ncol(naive.cytotoxic)) )

process_data <- function(sce) {
  # Define mitochondrial, ribosomal protein genes
  rowData(sce)$Mito <- grepl("^MT", rowData(sce)$Symbol)
  rowData(sce)$ProteinRibo <- grepl("^RP", rowData(sce)$Symbol)

  # Calculate quality statistics
  sce <- addPerCellQC(sce, subsets=list(mt=rowData(sce)$Mito, rbp=rowData(sce)$ProteinRibo))
  sce <- addPerFeatureQC(sce)

  # Delete bad cells (low library size, low gene count, high mitochondrial read proportion)
  libsize.drop <- isOutlier(sce$sum, nmads = 3, type = "lower", log = TRUE)
  feature.drop <- isOutlier(sce$detected, nmads = 3, type = "both", log = TRUE)
  mito.drop <- isOutlier(sce$subsets_mt_percent, nmads = 3, type = "higher", log = TRUE)
  sce <- sce[,!(libsize.drop | feature.drop |mito.drop)]

  # Normalize library sizes, then log2 transformation with pseudo-count of 1
  sce <- logNormCounts(sce)

  # Subset to top 500 average count genes
  X <- t(logcounts(sce))
  X <- as.matrix(X)
  X <- X[,  order(rowData(sce)$mean, decreasing=T)[1:500]]
  return(list(data=X, labels=colData(sce)$CellType))
}

set.seed(1)
processed.t <- process_data(pure.t)
ss.id.null <- sample(1:nrow(processed.t$data), 1000)
X1 <- processed.t$data[ss.id.null, ]
X1.ss <- processed.t$data[setdiff(1:nrow(processed.t$data), ss.id.null), ]

processed.mix5 <- process_data(mix5)
tcell.id <- which(processed.mix5$labels == "Memory T cell")
bcell.id <- which(processed.mix5$labels == "B cell")
mono.id <- which(processed.mix5$labels == "Monocyte")
nk.id <- which(processed.mix5$labels == "Nature Killer")
naive.id <- which(processed.mix5$labels == "Naive Cytotoxic")

ss.id <- c(sample(tcell.id, 400), sample(bcell.id, 400),
           sample(mono.id, 400),
           sample(nk.id, 400), sample(naive.id, 400))

X2 <- processed.mix5$data[ss.id, ]
X2.ss <- processed.mix5$data[setdiff(1:nrow(processed.mix5$data), ss.id), ]

# we first analyze the ``negative control'' data
eigen_cov_X <- eigen(coop::covar(X1))
U <- eigen_cov_X$vectors
D <- diag(sqrt(1/(eigen_cov_X$values+0.01)))
X1_iso <- X1%*%t(U)%*%D%*%(U)

set.seed(1234)
# pick best seed out of 30 random initialziations
seed_list <- sample(2021,size=30,replace = FALSE)
k_means_list <- lapply(seed_list, function(x)kmeans_estimation(X1_iso,
                                                               k=5,1,iter.max = 30,seed = x))
within_ss_list <- lapply(k_means_list, function(x)x$objective[[x$iter]])
best_seed <- seed_list[which.min(unlist(within_ss_list))]
current_kmeans_k_3_neg <- k_means_list[[which.min(unlist(within_ss_list))]]
final_cluster_k_3_neg <- current_kmeans_k_3_neg$cluster[[current_kmeans_k_3_neg$iter]]
table(final_cluster_k_3_neg)
# we have the final cluster

# Visualize data and estimated clusters -- Figure 10 in the paper
# this is just pcr
set.seed(2021)
X1_umap = umap(X1_iso)
umap_data_x1 <- X1_umap$layout
colnames(umap_data_x1) <- c("UMAP1","UMAP2")

png(paste0(plot_output_dir,'Figure_10_a.png'),
    width = 6,height=6,res=600,units='in')

ggplot(umap_data_x1, aes(x=UMAP1, y=UMAP2,
                         colour=as.factor(final_cluster_k_3_neg))) +
  geom_point() +
  #coord_fixed(ratio=1.38/3.33) +
  #ggtitle("Memory T cells only") +
  ylab("UMAP2") +
  xlab("UMAP1") +
  scale_colour_brewer(palette="Dark2", name="Estimated clusters",
                      guide=guide_legend(nrow=2))  +
  theme_bw(base_size=18)  +
 # theme(legend.position="bottom")+
  theme(plot.title = element_text(hjust = 0.5,size=18),
        legend.position="bottom",
        #legend.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18,hjust = 0),
        axis.title=element_text(size=18))


dev.off()

estimate_MED <- function(X){
  for (j in c(1:ncol(X))){
    X[,j] <- X[,j]-median(X[,j])}
  sigma_hat <- sqrt(median(X^2)/qchisq(1/2,df=1))
  return(sigma_hat)
}

sig_hat_neg <- estimate_MED(X1_iso)

neg_pval_1_2 <- kmeans_inference(X=X1_iso, k=5,
                                 cluster_1=1,
                                 cluster_2=2,
                                 sig = sig_hat_neg,
                                 iter.max = 30,
                                 seed = best_seed)

neg_pval_1_3 <- kmeans_inference(X=X1_iso, k=5,
                                 cluster_1=1,
                                 cluster_2=3,
                                 sig = sig_hat_neg,
                                 iter.max = 30,
                                 seed = best_seed)

neg_pval_1_4 <- kmeans_inference(X=X1_iso, k=5,
                                 cluster_1=1,
                                 cluster_2=4,
                                 sig = sig_hat_neg,
                                 iter.max = 30,
                                 seed = best_seed)

neg_pval_1_5 <- kmeans_inference(X=X1_iso, k=5,
                                 cluster_1=1,
                                 cluster_2=5,
                                 sig = sig_hat_neg,
                                 iter.max = 30,
                                 seed = best_seed)

neg_pval_2_3 <- kmeans_inference(X=X1_iso, k=5,
                                 cluster_1=2,
                                 cluster_2=3,
                                 sig = sig_hat_neg,
                                 iter.max = 30,
                                 seed = best_seed)

neg_pval_2_4 <- kmeans_inference(X=X1_iso, k=5,
                                 cluster_1=2,
                                 cluster_2=4,
                                 sig = sig_hat_neg,
                                 iter.max = 30,
                                 seed = best_seed)

neg_pval_2_5 <- kmeans_inference(X=X1_iso, k=5,
                                 cluster_1=2,
                                 cluster_2=5,
                                 sig = sig_hat_neg,
                                 iter.max = 30,
                                 seed = best_seed)


neg_pval_3_4 <- kmeans_inference(X=X1_iso, k=5,
                                 cluster_1=3,
                                 cluster_2=4,
                                 sig = sig_hat_neg,
                                 iter.max = 30,
                                 seed = best_seed)


neg_pval_3_5 <- kmeans_inference(X=X1_iso, k=5,
                                 cluster_1=3,
                                 cluster_2=5,
                                 sig = sig_hat_neg,
                                 iter.max = 30,
                                 seed = best_seed)


neg_pval_4_5 <- kmeans_inference(X=X1_iso, k=5,
                                 cluster_1=4,
                                 cluster_2=5,
                                 sig = sig_hat_neg,
                                 iter.max = 30,
                                 seed = best_seed)

summary(neg_pval_1_2)
summary(neg_pval_1_3)
summary(neg_pval_1_4)
summary(neg_pval_1_5)
summary(neg_pval_2_3)
summary(neg_pval_2_4)
summary(neg_pval_2_5)
summary(neg_pval_3_4)
summary(neg_pval_3_5)
summary(neg_pval_4_5)


## now move on to the ``cluster'' data!
eigen_cov_X <- eigen(coop::covar(X2))
U <- eigen_cov_X$vectors
D <- diag(sqrt(1/(eigen_cov_X$values+0.01)))
X2_iso <- X2%*%t(U)%*%D%*%(U)

set.seed(2021)
X2_umap = umap(X2_iso)
umap_data_x2 <- data.frame(X2_umap$layout)
colnames(umap_data_x2) <- c("UMAP1","UMAP2")
umap_data_x2$CellType <- colData(mix5)$CellType[ss.id]

# median
sig_hat_pos <- estimate_MED(X2_iso)

# look at k means
set.seed(1234)
# pick best seed out of 30 random initializations
seed_list <- sample(2021,size=30,replace = FALSE)
k_means_list <- lapply(seed_list, function(x)
  kmeans_estimation(X2_iso,k=5,iter.max = 20,seed = x))

within_ss_list <- lapply(k_means_list, function(x)x$objective[[x$iter]])
best_seed_pos <- seed_list[which.min(unlist(within_ss_list))]
current_kmeans_k_3_pos <- k_means_list[[which.min(unlist(within_ss_list))]]
final_cluster_k_3_pos <- current_kmeans_k_3_pos$cluster[[current_kmeans_k_3_pos$iter]]

mclust::adjustedRandIndex(final_cluster_k_3_pos, umap_data_x2$CellType)
result_tab <- table(final_cluster_k_3_pos, umap_data_x2$CellType)
result_tab/rowSums(result_tab)

png(paste0(plot_output_dir,'Figure_10_b.png'),
    width = 6,height=6,res=600,units='in')    

ggplot(umap_data_x2, aes(x=UMAP1, y=UMAP2, colour=as.factor(final_cluster_k_3_pos))) +
  geom_point() +
 # coord_fixed(ratio=6.84/18.23) +
 # ggtitle("Five cell type admixture") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  scale_colour_brewer(palette="Dark2", name="Estimated clusters",
                      guide=guide_legend(nrow=2))  +
  theme_bw(base_size=18)  +
  #coord_fixed(ratio=1/2)+
  # theme(legend.position="bottom")+
  theme(plot.title = element_text(hjust = 0.5,size=18),
        legend.position="bottom",
        #legend.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18,hjust = 0),
        axis.title=element_text(size=18))

dev.off()



pos_pval_1_2 <- kmeans_inference(X=X2_iso, k=5,
                                 cluster_1=1,
                                 cluster_2=2,
                                 sig = sig_hat_pos,
                                 iter.max = 20,
                                 seed = best_seed_pos)
pos_pval_1_3 <- kmeans_inference(X=X2_iso, k=5,
                                 cluster_1=1,
                                 cluster_2=3,
                                 sig = sig_hat_pos,
                                 iter.max = 20,
                                 seed = best_seed_pos)
pos_pval_1_4 <- kmeans_inference(X=X2_iso, k=5,
                                 cluster_1=1,
                                 cluster_2=4,
                                 sig = sig_hat_pos,
                                 iter.max = 20,
                                 seed = best_seed_pos)
pos_pval_1_5 <- kmeans_inference(X=X2_iso, k=5,
                                 cluster_1=1,
                                 cluster_2=5,
                                 sig = sig_hat_pos,
                                 iter.max = 20,
                                 seed = best_seed_pos)
pos_pval_2_3 <- kmeans_inference(X=X2_iso, k=5,
                                 cluster_1=2,
                                 cluster_2=3,
                                 sig = sig_hat_pos,
                                 iter.max = 10,
                                 seed = best_seed_pos)
pos_pval_2_4 <- kmeans_inference(X=X2_iso, k=5,
                                 cluster_1=2,
                                 cluster_2=4,
                                 sig = sig_hat_pos,
                                 iter.max = 10,
                                 seed = best_seed_pos)

pos_pval_2_5 <- kmeans_inference(X=X2_iso, k=5,
                                 cluster_1=2,
                                 cluster_2=4,
                                 sig = sig_hat_pos,
                                 iter.max = 10,
                                 seed = best_seed_pos)

pos_pval_3_4 <- kmeans_inference(X=X2_iso, k=5,
                                 cluster_1=3,
                                 cluster_2=4,
                                 sig = sig_hat_pos,
                                 iter.max = 20,
                                 seed = best_seed_pos)
pos_pval_3_5 <- kmeans_inference(X=X2_iso, k=5,
                                 cluster_1=3,
                                 cluster_2=5,
                                 sig = sig_hat_pos,
                                 iter.max = 20,
                                 seed = best_seed_pos)

pos_pval_4_5 <- kmeans_inference(X=X2_iso, k=5,
                                 cluster_1=4,
                                 cluster_2=5,
                                 sig = sig_hat_pos,
                                 iter.max = 20,
                                 seed = best_seed_pos)

summary(pos_pval_1_2)
summary(pos_pval_1_3)
summary(pos_pval_1_4)
summary(pos_pval_1_5)
summary(pos_pval_2_3)
summary(pos_pval_2_4)
summary(pos_pval_2_5)
summary(pos_pval_3_4)
summary(pos_pval_3_5)
summary(pos_pval_4_5)

