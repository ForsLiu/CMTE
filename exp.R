setwd("C:/Users/yh95l/Desktop/CMTE")

unlink(list.files("Data", full.names = TRUE), recursive = TRUE, force = TRUE)
file.remove(list.files("results", full.names = TRUE))

source("./parameters.R")
source("./DataGen.R")
source("./Evaluation.R")

library(foreach)
library(doParallel)
library("TRES")

source("./alg_CMTE/cmte_exp.R")
source("./alg_TRR/trr_exp.R")
source("./alg_TMDDM/tmddm_exp.R")

#generate table
source("./table_output.R")