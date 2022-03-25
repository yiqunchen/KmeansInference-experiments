library(ggplot2)
library(tidyverse)
library(latex2exp)
library(ggpubr)

input_dir <- "~/Desktop/dissertation/k_means_project/kmeans_inference_experiment/cluster_output/"
setwd(input_dir)
plot_output_dir <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/plot_output/"


sim_files <- list.files(path = input_dir,
                        pattern = glob2rx("validation_sim_k_*delta_0*_p_*.*"),full.names = T)

concat_p_value_list <- list()

for (i in seq_along(sim_files)){
  current_file <- sim_files[i]
  load(file=current_file)
  concat_p_value_list[[i]] <- result_df
  rm(result_df)
}

df_p_value_list <- dplyr::bind_rows(concat_p_value_list)

alpha_thresh <- 0.05
# load power data


q_label_names <- c(
  `2` = "q=2",
  `10` = "q=10",
  `50` = "q=50",
  `100` = "q=100")


p_aggregate_type_1 <- df_p_value_list %>%
  filter(delta==0,sigma==1,p<200) %>%
  pivot_longer(c(naive_vec,selective_vec,selective_vec_sample,selective_vec_MED),
               names_to = "p_type", values_to = "p_values") %>%
  mutate(p_type = fct_relevel(as.factor(p_type), "naive_vec", after = Inf)) %>%
  group_by(p,p_type) %>%
  mutate(theoretical = ecdf(p_values)(p_values)
         ) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(y = p_values, x = (theoretical), colour=p_type),size=0.6) +
  facet_wrap(.~p, nrow=1,scales = "free",
             labeller = as_labeller(q_label_names))+
  ylab('P-value quantiles')+
  xlab('Uniform(0,1) quantiles')+
  #ggtitle(TeX('Selective test p-values'))+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1)+
  scale_colour_brewer(name="Tests",palette="Dark2",
                      labels = unname(TeX(c('$p_{selective}$',
                                            '$\\hat{p}_{selective}(\\hat{\\sigma}_{MED})$',
                                            '$\\hat{p}_{selective}(\\hat{\\sigma}_{Sample})$',
                                            '$p_{Naive}$'))))+
  theme_bw(base_size=17) +
  theme(plot.title = element_text(hjust = 0.5,size=17),
        legend.position="none",
        legend.title = element_text(size=17),
        axis.text = element_text(size= 17),
        legend.text = element_text(size=17,hjust = 0),
        axis.title=element_text(size=17))+
  guides(colour = guide_legend(override.aes = list(size=5)))


png(paste0(plot_output_dir,'Figure_4.png'),
    width = 13 ,height=3.4,res=400,units='in')
p_aggregate_type_1
dev.off()


