# one-line example to run GEN or GFL
# load the data set generated in the simulation section
#
# Jie Jian
# 2022 Aug
#
load("./data/simulationData_50_1.RData")

################## GEN ##################
source("./R/GEN_TVNet_utilities.R")
source("./R/GEN_TVNet.R")
library(MASS)
library(Matrix)

GEN_result=TVNet_GEN(YY,lambda1=0.1,lambda2=0.1,
                     alpha=10,tol_rho=1e-8,
                     tol_sigma=1e-8,tol_admm=1e-8)
# the result contains "est_rho", the estimated 
# partial correlation matrices at each time point

################## GFL ##################
source("./R/GFL_admm_zupdate.R")
source("./R/GFL_TVNet_utilities.R")
source("./R/GFL_TVNet.R")

GFL_result=TVNet_GFL_ADMM(YY,lambda1=0.1,lambda2=0.1,
                          alpha=10,tol_rho=1e-8,
                          tol_sigma=1e-8,tol_admm=1e-8)
# the result contains "est_rho", the estimated 
# partial correlation matrices at each time point