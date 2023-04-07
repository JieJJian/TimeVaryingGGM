###################################################
#
# To simulate the scenario 1 with replications via GFL
#
# Jie Jian
# 2021 March
#
###################################################

source("GFL_admm_zupdate.R")
source("GFL_TVNet_utilities.R")
source("GFL_TVNet.R")
library(MASS)
library(Matrix)
load("s1_block_xi_true_t30.RData")

se1=0.1
se2=0.5

NumRepeat=100

m <- 50 # total participants in the experiment (number of subject)

# attribute of the network
k <- length(adj_i) # total time point
p <- ncol(true_rho[[1]]) # total test point in each brain
bs.k <- ncol(B)

t.start <- Sys.time()
result=vector("list",length = NumRepeat)

for (data.num in 1:NumRepeat){
  
  set.seed(data.num)
  
  xi <- vector("list",length = bs.k)
  # where random comes from
  for (i in 1:bs.k) {
    xi[[i]] <- mvrnorm(m,rep(0,p),xi_cov[[i]],empirical = TRUE)
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
  
  GFL_result <- TVNet_GFL_ADMM(YY,se1,se2,alpha=10,tol_rho=1e-5,tol_sigma=1e-5,tol_admm=1e-5)
 
  GFL_data=list("YY"=YY,"sample_con"=sample_concentration,"sample_rho"=sample_rho,
            "sample_cov"=sample_covariance)
  
  result[[data.num]]=list("GFL_result"=GFL_result,"data"=GFL_data)
  
}

save(result,file = paste0("GFL_simulation1_repeat_",m,".RData"))

t.end <- Sys.time()
print(t.end-t.start)