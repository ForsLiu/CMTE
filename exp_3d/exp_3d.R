setwd("C:/Users/yh95l/Desktop/CMTE")

unlink(list.files("Data", full.names = TRUE), recursive = TRUE, force = TRUE)
file.remove(list.files("results", full.names = TRUE))

source("./exp_3d/parameters_3d.R")
print(nrow(param_grid))
set.seed(123)
source("./exp_3d/DataGen_3d.R")
source("./Evaluation.R")

library(foreach)
library(doParallel)
library("TRES")

source("./alg_CMTE/cmte_exp.R")
source("./alg_TRR/trr_exp.R")
source("./alg_TMDDM/tmddm_exp.R")

#generate table
source("./exp_3d/table_output_3d.R")
