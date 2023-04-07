# solve the flsa problem
# input: response vector with the length of k * p(p-1)/2 
# solve the problem of min ||resY - z||^2_2 + \mu_1 ||z||_1 + \mu_2 ||Dz||_1

library(flsa)

admm_flsa_zupdate = function(resY,k,mu1,mu2){
  
  P = length(resY)/k
  
  matY = matrix(resY,k,P,byrow = TRUE)
  matZ = matrix(0,k,P)
  
  for (i in 1:P) {
    # matZ[,i]=ifelse(mu1==0,
    #                 as.vector(flsa(matY[,i], lambda1=mu1, lambda2 =mu2)),
    #                 flsa(matY[,i], lambda1=mu1, lambda2 =mu2)[1,1,])
    
    if(mu1==0){
      matZ[,i]=as.vector(flsa(matY[,i], lambda1=mu1, lambda2 =mu2))
    } else {
      matZ[,i]=flsa(matY[,i], lambda1=mu1, lambda2 =mu2)[1,1,]
    }
  }
  
  return(as.vector(t(matZ)))
}