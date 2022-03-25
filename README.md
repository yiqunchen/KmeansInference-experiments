# KmeansInference-experiments
Code and instructions to reproduce results and figures from the paper "Selective inference for k-means clustering"

Chen YT and Witten DM. (2022+) [Selective inference for k-means clustering](https://arxiv.org/abs/). arXiv:TBD [statME].

### Step-by-step instructions for reproducing figures
`~/Desktop/dissertation/gflasso_project/KmeansInference-experiment/`. Please double check your working directories and modify the relevant lines in the code accordingly (e.g., if you are running them on a computing cluster).

1. Download this repository and install relevant packages.
```
git clone https://github.com/yiqunchen/KmeansInference-experiments.git
cd KmeansInference-experiments
Rscript install_packages.R
```
2. Generate Figure 1 (an example to illustrate the inflated Type I error of the naive approach, and the selective Type I error control of our proposal).
```
Rscript ./Figure_1_Gen_Data.R
Rscript ./Figure_1_plot.R
```
2. Generate Figure 2 (an example to illustrate the inflated Type I error of the naive sample splitting approach).
```
Rscript ./Figure_2_Gen_Data.R
Rscript ./Figure_2_plot.R
```
3. Generate Figure 3 (intuition for the perturbation x'(phi)).
```
Rscript ./Figure_3_plot.R
```
4. Generate Figure 4 (demonstrate selective Type I error control). This step can be computationally intensive! The first R script runs the simulation study  described in Section 5.2 of our manuscript. This step might be computationally intensive!
```
Rscript ./sim_type_1_data.R
Rscript ./sim_type_1_plot.R
```
5. Generate Figures 5 and 8 (detection probability, conditional power, and the regression-spline-based estimate of power). The first R script runs the simulation study  described in Section 5.3 of our manuscript. This step might be computationally intensive!
```
Rscript ./sim_power_data.R 
Rscript ./detect_p_power_plot.R
```
6. Generate Figure 6 (analysis of the Palmer penguins data). You need to install the `palmerpenguins` package.
```
Rscript ./Figure_6.R
```
7. Generate Figure 7 (analysis of the MNIST data). You need the `keras` package installed.
```
Rscript ./Figure_7.R
```
8. Generate Table 1 and Figure 10 (analysis of the single-cell data). You need to install the relevant single-cell packages in `install_packages.R`, in addition to download the single cell datasets from Zheng et al. (2017).
In the code block below, we will make use of the data of the following five cell types: [monocytes](https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd14_monocytes/cd14_monocytes_filtered_gene_bc_matrices.tar.gz), [memory T cells](https://cf.10xgenomics.com/samples/cell-exp/1.1.0/memory_t/memory_t_filtered_gene_bc_matrices.tar.gz), [naive T cells](https://cf.10xgenomics.com/samples/cell-exp/1.1.0/naive_t/naive_t_filtered_gene_bc_matrices.tar.gz), [B cells](https://cf.10xgenomics.com/samples/cell-exp/1.1.0/b_cells/b_cells_filtered_gene_bc_matrices.tar.gz), and [natural killer cells](https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd56_nk/cd56_nk_filtered_gene_bc_matrices.tar.gz).

After downloading the data, you will need to need to rename the five downloaded folders as `filtered_matrices_mex_mono`,  `filtered_matrices_mex_memory`, `filtered_matrices_mex_naive_cytotoxic`,`filtered_matrices_mex_bcell`, and `filtered_matrices_mex_nk`. The following analysis also assumes that all five folders are listed under the directory `raw`.
```
Rscript ./Table_1_and_Figure_10.R
```
9. Generate Figure 9. The first R script runs the simulation study  described in Appendix A.8 of our manuscript and might be computationally intensive. 
```
Rscript ./sim_power_data_less_sparse.R
Rscript ./detect_p_power_plot_q_50.R
Rscript ./detect_p_power_plot_less_sparse.R
```

### Reference: 
- Chen YT and Witten DM. (2022+) Selective inference for k-means clustering. arXiv preprint. https://arxiv.org/abs/xxxx.xxxxx.
- Zheng, G. X. Y., Terry, J. M., Belgrader, P., Ryvkin, P., Bent, Z. W., Wilson, R., Ziraldo, S. B., Wheeler, T. D., McDermott, G. P., Zhu, J., Gregory, M. T., Shuga, J., Montesclaros, L., Underwood, J. G., Masquelier, D. A., Nishimura, S. Y., Schnall-Levin, M., Wyatt, P. W., Hindson, C. M., Bharadwaj, R., Wong, A., Ness, K. D., Beppu, L. W., Deeg, H. J., McFarland, C., Loeb, K. R., Valente, W. J., Ericson, N. G., Stevens, E. A., Radich, J. P., Mikkelsen, T. S., Hindson, B. J., and Bielas, J. H. (2017), “Massively parallel digital transcriptional profiling of single cells,” Nature Communications, 8, 14049. https://doi.org/10.1038/ncomms14049.




