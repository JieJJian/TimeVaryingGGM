###################################################
#
# GFL_TVNet_utilities
#
# Jie Jian
# 
# create on 2021-Jan-27
# updated on 2021-Feb-13 (add pml)
# update on March 27: modify ADMM algorithm to generate sparsity on rho and difference
#
###################################################



#
# "uniform_weight" function
# input: p
# output: vector of weights
uniform_weight <- function(p){
  w <- rep(1,p)
  return(w/sqrt(sum(w^2)))
}


#
# "sigma_initial" function
# input: matrix Y whose columns are observations of each covariate.
# output: initial estimate of vector sigma whose elements are the 
# diagnol of concentration matrix;
sigma_initial <- function(Y){
  return(1/diag(cov(Y)))
}



#
# "X_generation" function
# generate the 
# input: p (number of test position);
#        k (number of time point).
# output: blocks

X_generation <- function(sigma_new,YY_tilde,lambda,alpha,includeR=TRUE){
  
  # attribute of the network
  k <- length(YY_tilde) # total time point
  p <- ncol(YY_tilde[[1]]) # total test point in each brain
  m <- nrow(YY_tilde[[1]]) # total participants in the experiment
  
  X_list <- vector("list",length = k)
  R_list <- vector("list",length = k)  
  
  for (i in 1:k){
    
    Xt <- matrix(0,m*p,p*(p-1)/2)
    
    sigma <- sigma_new[,i]
    # weight <- weight_mat[,i]
    Y <- YY_tilde[[i]]
    
    for (j in 1:(p-1)){
      # Diagonal part
      diag_mat <- matrix(0,m,(p-j))
      for (h in 1:(p-j)){
        diag_mat[,h] <- Y[,j+h]*sqrt(sigma[j+h])/sqrt(sigma[j])
      }
      Xt[ ((j-1)*m+1):(j*m) , ( (2*p-j-1)*j/2 -p+j+1 ):((2*p-j-1)*j/2) ] <- diag_mat
      # Lower part
      for (h in (j+1):p){
        lower_mat <- matrix(0,m,p-j)
        lower_mat[,(h-j)] <- Y[,j]*sqrt(sigma[j])/sqrt(sigma[h])
        Xt[ ((h-1)*m+1):(h*m) , ( (2*p-j-1)*j/2 -p+j+1 ):((2*p-j-1)*j/2) ] <- lower_mat
      }
    }
    
    X_list[[i]] <- Xt
    
    if(includeR==TRUE){
      if (lambda!=0){
        if (i==1 || i==k ){
          A <- (t(Xt)%*%Xt/(m*lambda) + (1+alpha/(2*lambda) )*diag(p*(p-1)/2) )
        } else {
          A <- (t(Xt)%*%Xt/(m*lambda) + (2+alpha/(2*lambda) )*diag(p*(p-1)/2) )
        }
        
        if (i==1){
          R_list[[i]] <- solve(A)
        } else {
          R_list[[i]] <- solve(A-R_list[[(i-1)]])
        }
        
      } else {
        A <- (2*t(Xt)%*%Xt/m + alpha*diag(p*(p-1)/2) )
        R_list[[i]] <- solve(A)
      }
    } else {
      R_list=0
    }
    
    
  }
  result <- list("X"=X_list,"R"=R_list)
  return(result)
}

#
# ADMM algorithm
ADMM_3 <- function(Xlist,Ylist,alpha,tol,lamb1,lamb2){
  
  # attribute of the network
  k <- length(Ylist) # total time point
  p <- ncol(Ylist[[1]]) # total test point in each brain
  m <- nrow(Ylist[[1]]) # total participants in the experiment
  
  XY <- XY_cal(Xlist,Ylist) # calculate t(X)%*%Y
  Rlist <- vector("list",k) # calculate Rlist
  for (i in 1:k) {
    Rlist[[i]]=solve(2*t(Xlist[[i]])%*%Xlist[[i]]/m+alpha*diag(p*(p-1)/2))
  }
  
  count_ADMM <- 0
  
  # initial no warm start
  z_old <- rep(0,k*p*(p-1)/2)
  u <- rep(0,k*p*(p-1)/2)
  rho_old <- rep(0,k*p*(p-1)/2)
  
  # round1
  longVec <- (2/m)*XY+alpha*(z_old-u)
  rho_new <- admm3_x_update(Rlist,longVec,alpha,m)
  z_new <- admm_flsa_zupdate(rho_new+u,k,lamb1/alpha,lamb2/alpha)
  u <- u+rho_new-z_new
  
  while  ( (norm((rho_new - z_new),"2")> tol) || 
           (norm(alpha*(z_new-z_old),"2")> tol) ) {
    count_ADMM <- count_ADMM + 1
    rho_old <- rho_new
    z_old <- z_new
    
    longVec <- (2/m)*XY+alpha*(z_old-u)
    rho_new <-  admm3_x_update(Rlist,longVec,alpha,m)
    z_new <- admm_flsa_zupdate(rho_new+u,k,lamb1/alpha,lamb2/alpha)
    u <- u+rho_new-z_new
    
  }
  result_ADMM <- list("rho" = z_new[1:(k*p*(p-1)/2)], "iteration_ADMM" = count_ADMM)
  return(result_ADMM)
}


#admm3_x_update 
admm3_x_update <- function(Rlist,longVec,lamb2,m){
  l <- ncol(Rlist[[1]])
  k <- length(Rlist)
  
  for (i in 1:k){
    longVec[(1+l*(i-1)):(l*i)]<-Rlist[[i]]%*%longVec[(1+l*(i-1)):(l*i)]
  }
  
  return(longVec)
}



# long vector in the x-update
XY_cal <- function(XList,YList){
  p <- ncol(YList[[1]])
  k <- length(XList)
  l <- p*(p-1)/2
  
  vector_retrun <- rep(0,(k*l))
  
  for(i in 1:k){
    vector_retrun[(1+l*(i-1)):(l*i)] <- t(XList[[i]])%*%as.vector(YList[[i]])
  }
  return(vector_retrun)
}



# 
# "soft_th_vector" function
# soft thresholding operator
# Input: threshold and input in the bracket
# Output: a vector whose elements are the output of the soft thresholding
soft_th_vector <- function(threshold,x_vector){
  x <- ifelse(abs(x_vector)<=threshold,0,x_vector)
  x <- ifelse(x>threshold,(x-threshold),x)
  x <- ifelse(x<(-threshold),(x+threshold),x)
  return(x)
}


#
# "sigma_update" function
# input: matrix Y whose columns are observations of each covariate;
#        updated theta (vector);
#        old sigma (vector).
# output: updated sigma
sigma_update <- function(Y,theta_new,sigma0){
  # convert the vector theta to matrix theta
  p <- ncol(Y)
  t <- c()
  for (i in 1:(p-1)) {
    tt <- c()
    tt <- cbind(matrix(0,1,i),matrix(theta_new[(1+(i-1)*(2*p-i)/2):((i-1)*(2*p-i)/2+p-i)],1,))
    t <- rbind(t,tt)
  }
  t <- rbind(t,matrix(0,1,p))
  theta_new <- t
  
  B <- beta(theta_new,sigma0)
  diag(B) <- rep(-1,length(sigma0))  
  sigma_update <- vector()
  for (i in 1:length(sigma0)) {
    b <- B[i,]
    s <- nrow(Y)/sum((Y%*%b)^2)
    sigma_update <- append(sigma_update,s)
  }
  return(sigma_update)
}

#
# "beta" function
# input: theta (matrix);
#        sigma (vector).
# output: beta
beta <- function(theta,sigma){
  beta <- matrix(0,nrow(theta),ncol(theta))
  for (i in 1:nrow(theta)) {
    for(j in i:ncol(theta)){
      beta[i,j] <- theta[i,j]*sqrt(sigma[j])/sqrt(sigma[i])
      beta[j,i] <- theta[i,j]*sqrt(sigma[i])/sqrt(sigma[j])
    }
  }
  return(beta)
}


#
# turn the rho vector to k upper triangular matrix

as_upper_triangular <- function(rho,k,p){
  l <- p*(p-1)/2
  rhorho <- vector("list",length = k)
  for (i in 1:k){
    rho_t <- matrix(0,p,p)
    vec_t <- rho[((i-1)*l+1):(i*l)]
    for (j in 1:(p-1)) {
      rho_t[j,((j+1):p)] <- vec_t[((2*p-j)*(j-1)/2+1):((2*p-j-1)*j/2)]
    }
    rhorho[[i]] <- rho_t
  }
  return(rhorho)
}

# partial correlation and diagonal concentration to concentration matrix
# input: partial correlation matrix and diagonal entry of concentration matrix
# output: full concentration matrix
partialCor_to_concentration <- function(partialCor,diag.con){
  
  p <- nrow(partialCor)
  
  conMat <- matrix(0,p,p)
  for (ii in 1:(p-1)) {
    for (jj in (ii+1):p) {
      conMat[ii,jj] <- -partialCor[ii,jj]*sqrt(diag.con[ii]*diag.con[jj])
    }
  }
  
  conMat <- t(conMat) + conMat
  diag(conMat) <- diag.con
  
  return(conMat)
  
}



# concentration to partial correlation function
# input: concentration matrix
# output: full partial correlation matrix
concentration_to_partialCor <- function(conMat){
  
  p <- nrow(conMat)
  
  partialCor <- matrix(0,p,p)
  for (ii in 1:(p-1)) {
    for (jj in (ii+1):p) {
      partialCor[ii,jj] <- -conMat[ii,jj]/sqrt(conMat[ii,ii]*conMat[jj,jj])
    }
  }
  
  partialCor <- t(partialCor) + partialCor
  diag(partialCor) <- rep(1,p)
  
  return(partialCor)
  
}


######### 1 norm of rho and difference of rho#############
# input: rho: a list of k upper triangular matrice

norm1 <- function(rho){
  # attribute of the network
  k <- length(rho) # total time point
  
  norm1_origin <- 0
  norm1_diff <- 0
  norm2_diff <- 0
  
  # turn each matrix in the list rho into a upper triangular
  for (i in 1:k){
    rho[[i]][lower.tri(rho[[i]],diag = TRUE)]==0
  }
  
  for (i in 1:k) {
    norm1_origin <- norm1_origin + sum(abs(rho[[i]]))
  }
  if (k>1){
    for (i in 2:k) {
      norm1_diff <- norm1_diff + sum(abs(rho[[i]]-rho[[(i-1)]]))
      norm2_diff <- norm2_diff + sum((rho[[i]]-rho[[(i-1)]])^2)
    }
  } else {
    norm1_diff <- NA
    norm2_diff <- NA
  }
  
  result <- list("origin_l1"=norm1_origin,"diff_l1"=norm1_diff,"diff_l2"=norm2_diff)
  
}


# rss
rss_cal <- function(XList,YList,rhoVector){
  p <- ncol(YList[[1]])
  k <- length(XList)
  l <- p*(p-1)/2
  
  rss=0
  
  for(i in 1:k){
    rss <- rss+sum((as.vector(YList[[i]])-XList[[i]]%*%rhoVector[((i-1)*l+1):(i*l)])^2)
  }
  return(rss)
}


# calculate the first term in the BIC criterion
pmlT1 <- function(est_con,sample_cov){
  k <- length(est_con)
  p <- ncol(est_con[[1]])
  
  term1 <- 0
  term2 <- 0
  if(sum(sapply(1:k, FUN = function(x){det(est_con[[x]])})>0)==k){
    for (i in 1:k){
      term1 <- term1-log(det(est_con[[i]]))+sum(diag(est_con[[i]]%*%sample_cov[[i]]))
    }
  } else {
    term1 <- NA
  }
  
  return(term1)
}




