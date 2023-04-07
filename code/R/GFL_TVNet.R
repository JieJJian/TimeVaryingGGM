###################################################
#
# Time-varying network via GFL
#
# Jie Jian
# 
# updated on Feb 6: output l2 diff and rss
# updated on Feb 23: add functions to calculate the cases when lambda1=0 or lambda2=0
# Generalized elastic net 
# update on March 27: modify ADMM algorithm to generate sparsity on rho and difference
# update on July 28: add case with no penalty (LR)
#
###################################################

TVNet_GFL <- function(YY,lambda1,lambda2,alpha,tol_rho=1e-5,tol_sigma=1e-5,tol_admm=1e-5){
  
  if (lambda2==0){
    if(lambda1==0){
      return(LR(YY,tol_rho,tol_sigma))
    }else{return(TVNet_GEN_l1(YY,lambda1,alpha,tol_rho,tol_sigma,tol_admm))}
  } else {
    return(TVNet_GFL_ADMM(YY,lambda1,lambda2,alpha,tol_rho,tol_sigma,tol_admm))
  }
}



# calculate the normal (regular) regression
# when lambda1=0 and lambda2=0
LR=function(YY,tol_rho=1e-5,tol_sigma=1e-5){
  # attribute of the network
  k <- length(YY) # total time point
  p <- ncol(YY[[1]]) # total test point in each brain
  m <- nrow(YY[[1]]) # total participants in the experiment
  
  # adjust the tolerances
  tol_rho <- tol_rho*k*p*(p-1)
  tol_sigma <- tol_sigma*p*k
  
  # center mean 
  for (i in 1:k){
    YY[[i]] <- scale(YY[[i]],scale=FALSE)
  }
  
  ############################# Step 1: Give initial and updated sigma ##########################
  sigma_old <- matrix(0,p,k)
  
  sigma_new <- matrix(0,p,k)
  for (i in 1:k){
    sigma_new[,i] <- sigma_initial(YY[[i]])
  }
  
  ##----------------------- Step 2 & 3: Update rho and sigma ----------------------##
  
  # initial rho
  rho_old <- rep(0,k*p*(p-1)/2)
  
  # updated rho
  reformX <- X_generation(sigma_new,YY,0,0,includeR=FALSE)
  Xt_list <- reformX$X
  rhoMatrix=matrix(0,p*(p-1)/2,k)
  for (i in 1:k) {
    rhoMatrix[,i]=as.vector(lm(as.vector(YY[[i]])~Xt_list[[i]])$coefficients)[2:(p*(p-1)/2+1)]
  }
  rho_new <- as.vector(rhoMatrix)
  
  # updated sigma
  
  count <- 0
  
  # record the lambda values
  iteration_ADMM_record <- vector()
  
  while  ( (sum((rho_new - rho_old)^2) > tol_rho) 
           | (norm((sigma_new - sigma_old),"F") > tol_sigma) ) {
    
    count <- count + 1 
    rho_old <- rho_new
    
    ############################# Step 3: Update sigma ########################
    sigma_old <- sigma_new
    for (i in 1:k){
      sigma_new[,i] <- sigma_update(YY[[i]],rho_new[((i-1)*p*(p-1)/2+1):(i*p*(p-1)/2)],sigma_old[,i])
    }
    
    ############################# Step 2: Update rho ##########################
    reformX <- X_generation(sigma_new,YY,0,0,includeR=FALSE)
    Xt_list <- reformX$X
    rhoMatrix=matrix(0,p*(p-1)/2,k)
    for (i in 1:k) {
      rhoMatrix[,i]=as.vector(lm(as.vector(YY[[i]])~Xt_list[[i]])$coefficients)[2:(p*(p-1)/2+1)]
    }
    rho_new <- as.vector(rhoMatrix)
    
   # iteration_ADMM_record <- c(iteration_ADMM_record,calc_rho$iteration_ADMM)
  }
  
  rhorho <- as_upper_triangular(rho_new,k,p) #list
  l1 <- norm1(rhorho)
  
  # number of fuse group
  rhoMat = matrix(rho_new,(p*(p-1)/2),k,byrow = FALSE)
  fusenon0=sum(rhoMat[,1]!=0)
  for (i in 2:k) {
    fusenon0=fusenon0+sum((rhoMat[,i]!=0)& (rhoMat[,i]!=rhoMat[,(i-1)]))
  }
  
  # number of non0
  non0 <- sum(rhoMat!=0)
  
  # rss
  # input: Xt_list,YY_tilde,rho_new
  RSS=rss_cal(Xt_list,YY,rho_new)
  
  est.rho <- vector("list",length = k)
  for (i in 1:k){
    est.rho[[i]] <- rhorho[[i]] + t(rhorho[[i]])
    diag(est.rho[[i]]) <- rep(1,p)
  }
  
  est.con <- vector("list",length = k)
  for (i in 1:k){
    est.con[[i]] <- partialCor_to_concentration(rhorho[[i]],sigma_new[,i])
  }
  
  result_List <- list("est_rho" = est.rho, "rho.vector"=rho_new,"sigma" = sigma_new,"est_con"=est.con,
                      "iteration_ADMM_record" = iteration_ADMM_record,
                      "number_iteration" = count,
                      "l1_origin"=l1$origin_l1,"l1_diff"=l1$diff_l1,
                      "l2_diff_square"=l1$diff_l2,"rss"=RSS,
                      "non0"=non0,"fusenon0"=fusenon0)
  
  return(result_List)
  
  
}


# when lambda1 and lambda2 are non-zero
TVNet_GFL_ADMM <- function(YY,lambda1,lambda2,alpha,tol_rho=1e-5,tol_sigma=1e-5,tol_admm=1e-5){
  
  # attribute of the network
  k <- length(YY) # total time point
  p <- ncol(YY[[1]]) # total test point in each brain
  m <- nrow(YY[[1]]) # total participants in the experiment
  
  # adjust the tolerances
  tol_rho <- tol_rho*k*p*(p-1)
  tol_sigma <- tol_sigma*p*k
  tol_admm <- tol_admm*k*p*(p-1)
  
  # center mean 
  for (i in 1:k){
    YY[[i]] <- scale(YY[[i]],scale=FALSE)
  }
  
  # # Weights of each test point for each participant
  # weight_mat <- matrix(0,p,k)
  # for (i in 1:k){
  #   weight_mat[,i] <- uniform_weight(p)
  # }
  
  # # YY tilde: original data multiplied by sqrt(weight)
  # YY_tilde <- vector(mode = "list", length = k)
  # for (i in 1:k){
  #   w <- weight_mat[,i]
  #   YY_tilde[[i]] <- YY[[i]] # initialize YY_tilde[[i]] as a matrix
  #   for (j in 1:p){
  #     YY_tilde[[i]][,j] <- YY[[i]][,j]*w[j]
  #   }
  #   #YY_tilde[[i]] <- YY_tilde[[i]]/sqrt(m) # 1/m
  # }
  
  ############################# Step 1: Give initial and updated sigma ##########################
  sigma_old <- matrix(0,p,k)
  
  sigma_new <- matrix(0,p,k)
  for (i in 1:k){
    sigma_new[,i] <- sigma_initial(YY[[i]])
  }
  
  ##----------------------- Step 2 & 3: Update rho and sigma ----------------------##
  
  # initial rho
  rho_old <- rep(0,k*p*(p-1)/2)
  
  # updated rho
  reformX <- X_generation(sigma_new,YY,lambda2,alpha,includeR=FALSE)
  Xt_list <- reformX$X
  calc_rho <- ADMM_3(Xt_list,YY,alpha,tol_admm,lambda1,lambda2)
  rho_new <- calc_rho$rho
  
  # updated sigma
  
  count <- 0
  
  # record the lambda values
  iteration_ADMM_record <- vector()
  
  while  ( (sum((rho_new - rho_old)^2) > tol_rho) 
           | (norm((sigma_new - sigma_old),"F") > tol_sigma) ) {
     
    count <- count + 1 
    rho_old <- rho_new
    
    ############################# Step 3: Update sigma ########################
    sigma_old <- sigma_new
    for (i in 1:k){
      sigma_new[,i] <- sigma_update(YY[[i]],rho_new[((i-1)*p*(p-1)/2+1):(i*p*(p-1)/2)],sigma_old[,i])
    }
    
    ############################# Step 2: Update rho ##########################
    reformX <- X_generation(sigma_new,YY,lambda2,alpha,includeR=FALSE)
    Xt_list <- reformX$X
    calc_rho <- ADMM_3(Xt_list,YY,alpha,tol_admm,lambda1,lambda2)
    rho_new <- calc_rho$rho
    
    iteration_ADMM_record <- c(iteration_ADMM_record,calc_rho$iteration_ADMM)
  }
  
  rhorho <- as_upper_triangular(rho_new,k,p) #list
  l1 <- norm1(rhorho)

  # number of fuse group
  rhoMat = matrix(rho_new,(p*(p-1)/2),k,byrow = FALSE)
  fusenon0=sum(rhoMat[,1]!=0)
  for (i in 2:k) {
    fusenon0=fusenon0+sum((rhoMat[,i]!=0)& (rhoMat[,i]!=rhoMat[,(i-1)]))
  }
  
  # number of non0
  non0 <- sum(rhoMat!=0)
  
  # rss
  # input: Xt_list,YY_tilde,rho_new
  RSS=rss_cal(Xt_list,YY,rho_new)
  
  est.rho <- vector("list",length = k)
  for (i in 1:k){
    est.rho[[i]] <- rhorho[[i]] + t(rhorho[[i]])
    diag(est.rho[[i]]) <- rep(1,p)
  }
  
  est.con <- vector("list",length = k)
  for (i in 1:k){
    est.con[[i]] <- partialCor_to_concentration(rhorho[[i]],sigma_new[,i])
  }
  
  result_List <- list("est_rho" = est.rho, "rho.vector"=rho_new,"sigma" = sigma_new,"est_con"=est.con,
                      "iteration_ADMM_record" = iteration_ADMM_record,
                      "number_iteration" = count,
                      "l1_origin"=l1$origin_l1,"l1_diff"=l1$diff_l1,
                      "l2_diff_square"=l1$diff_l2,"rss"=RSS,
                      "non0"=non0,"fusenon0"=fusenon0)
  
  return(result_List)
  
}

