setwd("C:/Users/yh95l/Desktop/CMTE")

file.remove(list.files("Data_1.3", full.names = TRUE))
file.remove(list.files("results_1.3", full.names = TRUE))

source("./parameters_1.3.R")
source("./DataGen_1.3.R")
source("./Evaluation.R")

library(foreach)
library(doParallel)
library("TRES")

source("./alg_CMTE/cmte_exp_1.3.R")
source("./alg_TRR/trr_exp_1.3.R")
source("./alg_TMDDM/tmddm_exp_1.3.R")

#generate table
source("./table_output_1.3.R")
