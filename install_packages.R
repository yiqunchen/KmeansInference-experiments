packages <- c("devtools", "tidyverse", "latex2exp", "RColorBrewer",
              "shades","patchwork","palmerpenguins","keras",
              "umap","DropletUtils","scater", "ggfortify")
install.packages(setdiff(packages, rownames(installed.packages())))  
if(!("KmeansInference" %in% rownames(installed.packages()))){
  devtools::install_github("yiqunchen/KmeansInference")
}
