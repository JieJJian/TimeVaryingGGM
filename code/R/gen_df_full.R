###################################################
#
# To calculate the degrees of freedom
#
# Jie Jian
# 2021 Feb
#
###################################################

se1=c(0,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,5)
se2=c(0,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,5,10,20,25,30,35,40,45,50,60,80,100)
n=50
eps=0.000001
load("s1_gen_50.RData")
load("simulationData_50_1.RData")
source("gen_df_full_utilities.R")

DF <- rep(0,length(result))

k = length(YY)
p = ncol(YY[[1]])

# center mean 
for (i in 1:k){
  YY[[i]] <- scale(YY[[i]],scale=FALSE)
}

for (ii in 1:length(result)) {
  cat("\n",ii," - time: ",as.character(Sys.time()), sep = "")
  GEN_result=result[[ii]]
  lambda2 <- se2[ifelse(ii%%length(se2)==0,length(se2),ii%%length(se2))]
  
  # term2
  if((lambda2!=0)&(GEN_result$non0!=0)){
    
    df=0
    
    for (i in 1:p) {
        
        matReform=reformMat(YY,GEN_result$sigma,GEN_result$est_rho,i,p)
        DDa=matReform$DD
        XXa=matReform$XX+eps*diag(ncol(matReform$XX))
        invMat=solve(XXa+n*lambda2*DDa)
        df=df+sum(diag(invMat%*%(XXa)))
    }
    
    DF[ii] = df
  } else {
    DF[ii] =  2*(GEN_result$non0)
  }
   
}
DF_full=matrix(DF,length(se1),length(se2),byrow = TRUE)
save(DF_full,file = "df_full.RData")

# Xa for each i (i=1,...,p) and t (t=1,...,T)


# rho for each i (i=1,...,p) and t (t=1,...,T)
