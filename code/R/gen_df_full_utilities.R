###################################################
#
# To calculate the degrees of freedom
#
# Jie Jian
# 
# created on 2021-Feb
# updated on Feb 6: output l2 diff and rss
# updated on Feb 23: add functions to calculate the cases when lambda1=0 or lambda2=0
# Generalized elastic net 
# 2021-Jun-12: delete the weight
#
###################################################

# Xa for each i (i=1,...,p) and t (t=1,...,T)
# rho for each i (i=1,...,p) and t (t=1,...,T)

# reformMat(YY,GEN_result$sigma,GEN_result$est_rho,i,p)
reformMat = function(y,sig,rhoMat,ip,p){
  k=length(y)
  
  XXalist=vector("list",k)
  Dl = matrix(0,(k-1)*(p-1),k*(p-1))
  for (ii in 1:((k-1)*(p-1))) {
    Dl[ii,ii] <- 1
    Dl[ii,(ii+p-1)] <- -1
  }
  deleteDl=c()
  
  Aind=c()
  #reform rho vector
  for (i in 1:k) {
    
    Y=y[[i]]
    sigma=sig[,i]
    X=matrix(0,n,p)
    
    for (j in 1:p) {
      X[,j]=Y[,j]*sqrt(sigma[j])/sqrt(sigma[ip])
    }
    
    rhoT=rhoMat[[i]][ip,]
    indRecord=(rhoT!=0)
    indRecord[ip]=FALSE
    
    XXalist[[i]]=t(X[,indRecord])%*%X[,indRecord]
    indRecord=indRecord[-ip]
    deleteDl=c(deleteDl,indRecord)
    Aind=c(Aind,sum(indRecord))
  }
  
  DD=t(Dl[,deleteDl])%*%Dl[,deleteDl]
  XX=matrix(0,sum(Aind),sum(Aind))
  for (j in 1:k) {
    if(Aind[j]!=0){
      if (j==1){
        XX[1:Aind[1],1:Aind[1]]=XXalist[[1]]
      } else {
        XX[(1+sum(Aind[1:(j-1)])):sum(Aind[1:j]),(1+sum(Aind[1:(j-1)])):sum(Aind[1:j])]=XXalist[[j]]
      }
    }
  }
  
  result=list("DD"=DD,"XX"=XX)
  return(result)
}