###################################################
#
# To simulate the scenario 1 with a single trial
#
# Jie Jian
# 2021 Feb
#
###################################################

source("GEN_TVNet_utilities.R")
source("GEN_TVNet.R")
library(MASS)
library(Matrix)
load("s1_block_xi_true_t30.RData")

se1=c(0,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,5)
se2=c(0,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,5,10,20,25,30,35,40,45,50,60,80,100)
  
m <- 50 # total participants in the experiment (number of subject)

# attribute of the network
k <- length(adj_i) # total time point
p <- ncol(xi_cov[[1]]) # total test point in each brain
bs.k <- ncol(B)

t.start <- Sys.time()

for (data.num in 1){
  
  set.seed(data.num)
  
  xi <- vector("list",length = bs.k)
  # where random comes from
  for (i in 1:bs.k) {
    xi[[i]] <- mvrnorm(m,rep(0,p),xi_cov[[i]]) # need to check
  }
  
  ##### generate the matrix at each time point ########
  YY <- vector("list",length = k)
  sample_concentration <- vector("list",length = k)
  sample_rho <- vector("list",length = k)
  sample_covariance <- vector("list",length = k)
  
  for (i in 1:k){
    random.fluc <- 0
    for (j in 1:bs.k){
      random.fluc <- random.fluc +  B[i,j]*xi[[j]]
    }
    YY[[i]] <- random.fluc + matrix(rep(1,m*p),m,p)*i
    sample_concentration[[i]] <- solve(cov(YY[[i]]))
    sample_rho[[i]] <- concentration_to_partialCor(sample_concentration[[i]])
  }
  
  for (i in 1:k){
    sample_covariance[[i]] <- cov(YY[[i]])
  }
  
  ##### end of data generation ####
  
  cat("\n","dataset",data.num," - time: ",as.character(Sys.time()), sep = "")
  save(YY,sample_concentration,sample_rho,sample_covariance,
       file = paste0("simulationData_",m,"_",data.num,".RData"))
  result=vector("list",length = length(se1)*length(se2))
  ind=0
  for (lam1 in 1:length(se1) ){
    for (lam2 in 1:length(se2) ){
      ind=ind+1
      GEN_result <- TVNet_GEN(YY,se1[lam1],se2[lam2],alpha=10,tol_rho=1e-8,tol_sigma=1e-8,tol_admm=1e-8)
      cat("\n","(l1,l2)=(",se1[lam1],",",se2[lam2],") - time: ",as.character(Sys.time()), sep = "")
      result[[ind]]=GEN_result
    }
  }
}
save(result,file = "s1_gen_50.RData")
t.end <- Sys.time()
print(t.end-t.start)