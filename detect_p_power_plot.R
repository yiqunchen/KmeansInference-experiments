library(ggplot2)
library(tidyverse)
library(latex2exp)
library(ggpubr)

input_dir <- "~/Desktop/dissertation/k_means_project/kmeans_inference_experiment/cluster_output/"
setwd(input_dir)
plot_output_dir <- "~/Desktop/dissertation/k_means_project/KmeansInference-experiments/plot_output/"


sim_files <- list.files(path = input_dir,
                        pattern = glob2rx("validation_sim_k*delta*_p_10.*"),full.names = T)

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

sigma_label_names <- c(
  `0.5` = TeX("$\\sigma=0.5$"),
  `1` = TeX("$\\sigma=1$"),
  `0.25` = TeX("$\\sigma=0.25$")
)


q_label_names <- c(
  `2` = "q=2",
  `10` = "q=10",
  `50` = "q=50",
  `100` = "q=100")


detect_prob_df <- df_p_value_list %>%
  filter(delta>1,p==10,sigma<2) %>%
  group_by(p,sigma,delta) %>%
  summarise(mean_detect_p = mean(adjust_rand_index==1,na.rm = T),
            sd_detect_p = sqrt(mean_detect_p*(1-mean_detect_p)/n())) %>%
  ungroup()

p_detect_prob <- detect_prob_df %>%
  ggplot(aes(x=delta,y=mean_detect_p,linetype=as.factor(sigma),group=as.factor(sigma)))+
  geom_point()+
  geom_pointrange(aes(ymin = mean_detect_p+sd_detect_p,
                      ymax = mean_detect_p-sd_detect_p))+
  geom_line()+
  ylab(TeX('Detection probability'))+
  xlab(TeX('Distance between true clusters: $\\delta$'))+
  theme_bw(base_size=18) +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        legend.position="none",
        legend.title = element_text(size=18),
        legend.key.size = unit(2.5, "lines"),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18,hjust = 0),
        axis.title=element_text(size=18))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_linetype_discrete(name=TeX('$\\sigma$'))

png(paste0(plot_output_dir,'Figure_4_a.png'),
    width =6 ,height=6,res=400,units='in')
p_detect_prob
dev.off()


cond_power_df <- df_p_value_list %>%
  filter(delta>1,p==10,sigma<2) %>%
  select(-naive_vec) %>%
  pivot_longer(c(selective_vec,selective_vec_sample,selective_vec_MED),
               names_to = "p_type", values_to = "p_values") %>%
  group_by(p,p_type,sigma,delta) %>%
  summarise(mean_cond_power = mean((p_values<=alpha_thresh)&(adjust_rand_index==1),na.rm = T)/mean(adjust_rand_index==1,na.rm = T),
            sd_cond_power = sqrt(mean_cond_power*(1-mean_cond_power)/n())) %>%
  ungroup()

p_cond_power <- cond_power_df %>%
  ggplot(aes(x=delta,y=mean_cond_power,
             colour=p_type,linetype=as.factor(sigma)))+
  geom_point()+
  geom_pointrange(aes(ymin = mean_cond_power+sd_cond_power,
                      ymax = mean_cond_power-sd_cond_power))+
  geom_line()+
  ylab(TeX('Conditional power at $\\alpha$=$0.05$'))+
  xlab(TeX('Distance between true clusters: $\\delta$'))+
  theme_bw(base_size=18) +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        legend.position="none",
        legend.title = element_text(size=15),
        axis.text = element_text(size=18),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=18))+
  guides(colour = guide_legend(override.aes = list(size=1)))+
  scale_colour_brewer(name="Tests",palette="Dark2",
                      labels = unname(TeX(c('$p_{selective}$',
                  '$\\hat{p}_{selective}(\\hat{\\sigma}_{MED})$',
                  '$\\hat{p}_{selective}(\\hat{\\sigma}_{Sample})$'))))+
  scale_linetype(guide="none")

png(paste0(plot_output_dir,'Figure_4_b.png'),
    width =6 ,height=6,res=400,units='in')
p_cond_power
dev.off()


smoothed_power <- df_smooth_power %>%
  filter(sigma<2,sigma>0.25) %>%
  ggplot(aes(x=(true_mu_diff)*sigma,
             y=as.numeric(p_values<=alpha_thresh),
             colour = p_type)) +
  facet_wrap(.~as.factor(sigma), nrow=1)+
  geom_smooth()+
  theme_bw() +
  ylab(TeX('Power at $\\alpha = 0.05$'))+
  xlab(TeX('$ || \\mu^T\\nu ||_2$'))+
  theme_bw(base_size=18) +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        legend.position="none",
        legend.title = element_text(size=15),
        axis.text = element_text(size=18),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=18))+
  scale_colour_brewer(name="Tests",palette="Dark2",
                      labels = unname(TeX(c('$p_{selective}$',
                                            '$\\hat{p}_{selective}(\\hat{\\sigma}_{MED})$',
                                            '$\\hat{p}_{selective}(\\hat{\\sigma}_{Sample})$'))))


png(paste0(plot_output_dir,'Figure_8.png'),
    width =8 ,height=4.3,res=400,units='in')
smoothed_power
dev.off()

