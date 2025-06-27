setwd("C:/Users/yh95l/Desktop/CMTE")

unlink(list.files("Data_1.1_1.2", full.names = TRUE), recursive = TRUE, force = TRUE)
file.remove(list.files("results_1.1_1.2", full.names = TRUE))

source("./parameters_1.1_1.2.R")
source("./DataGen_1.1_1.2.R")
source("./Evaluation.R")

library(foreach)
library(doParallel)
library("TRES")

source("./alg_CMTE/cmte_exp_1.1_1.2.R")
source("./alg_TRR/trr_exp_1.1_1.2.R")
source("./alg_TMDDM/tmddm_exp_1.1_1.2.R")

#generate table
source("./table_output_1.1_1.2.R")
