knots_eq3 <- function(x, k, m){
  #external knots are on boundary
  #return boundary with internal knots only
  #used in bs or bsplineS
  c(min(x), seq(from=min(x), to=max(x), length.out=m+2)[-c(1,m+2)], max(x))
}

create_adjacency <- function(V,n) {
  differences <- apply(V,2,FUN=function(x) {norm(as.matrix(x),'f')})
  connected_ix <- which(differences == 0);
  index = t(combn(n,2));
  i <- index[connected_ix,1]
  j <- index[connected_ix,2]
  A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  A[(j-1)*n + i] <- 1
  return(A)
}


#####x:vector，sigma:double#####


prox_L2<- function(x,sigma){
  n <- length(x)
  px <- rep(0,n)
  lv <- norm(x,type = "2")
  if(lv==0.){
    px <- x
  }
  else{
    px <- max(0.,1.-(sigma/lv))*x
  }
  return(px)
}


######UPDATA B#################

update_B <- function(diagD,V,Lambda,n,gamma1,index,theta,X,y,Hs,B_old,individual_time){
  m <- nrow(index)
  p <- nrow(V)
  B_new <- matrix(nrow=p,ncol=n)
  V_tilde <- V+Lambda/theta
  eyemat_n <- diag(n)
  eyemat_p <- diag(p)
  onemat_n <- ones(n)
  AtA <- kronecker(n*eyemat_n - onemat_n,eyemat_p)
  GHlist <- Grad_L(B_old,X,y,individual_time)
  Hs <- diag(GHlist[[2]])
  invmat <- solve(-Hs+gamma1*diagD+theta*AtA)
  lv <- rep(0,n*p)
  for (l in 1:m){
    i <- index[l,1]
    j <- index[l,2]
    e1 <- rep(0,n)
    e2 <- rep(0,n)
    e1[i] <- 1
    e2[j] <- 1
    Alt <- kronecker(e1-e2,eyemat_p)
    lv <- lv + Alt%*%V_tilde[,l]
  }
  Grad_L <- GHlist[[1]]
  rhs <- Grad_L-t(as.vector(B_old)%*%Hs)+theta*lv
  b_new <- invmat%*%rhs
  B_new <- matrix(b_new,nrow=p,ncol=n)
  return(B_new)
}

#####calculate the gradient of beta############################



Pupbeta <- function(B_old,X,y,individual_time,subject){
  trydata <- data.frame(exp(X%*%B_old),y,individual_time,X)
  n1<- ncol(exp(X%*%B_old))
  n2 <- ncol(y)
  n3 <- ncol(X)
  sub_trydata <- subset(trydata,individual==subject)
  pupbeta <- matrix(nrow=n1*n3,ncol=0)
  for(i in 1:length(sub_trydata$time)){
    yi0 <- sum(sub_trydata[i,(n1+1):(n1+n2)])
    alphasum <- rowSums(sub_trydata[i,1:n1])
    alphaij <- t(sub_trydata[i,1:n1])
    a<- alphaij/alphasum
    vai <- matrix(ai,ncol=1)
    tvai <- matrix(ai,nrow=1)
    Aij <- kronecker(yi0*(diag(a)-a%*%t(a)),t(sub_trydata[i,(n1+n2+3):(n1+n2+n3+2)]))
    pupbeta <- cbind(pupbeta,Aij)
    
    
  }
  return(pupbeta)
  
}

QIJ <- function(B_old,X,y,individual_time,subject){
  trydata <- data.frame(exp(X%*%B_old),y,individual_time,X)
  n1<- ncol(exp(X%*%B_old))
  n2 <- ncol(y)
  n3 <- ncol(X)
  qij <- matrix(nrow=n1*n3,ncol=0)
  sub_trydata <- subset(trydata,individual==subject)
  for(i in 1:length(sub_trydata$time)){
    yi0 <- sum(sub_trydata[i,(n1+1):(n1+n2)])
    alphasum <- rowSums(sub_trydata[i,1:n1])
    alphaij <- t(sub_trydata[i,1:n1])
    a<- alphaij/alphasum
    vai <- matrix(ai,ncol=1)
    tvai <- matrix(ai,nrow=1)
    Qij <- kronecker(yi0*(diag(rep(1,n1))-2*diag(a))%*%(diag(a)-a%*%t(a)),t((sub_trydata[i,(n1+n2+3):(n1+n2+n3+2)])^2))
    qij <- cbind(qij ,Qij)
  }
  return(qij)
}

sub_estimate_cov<- function(trydata,subject,time1,time2,n1,n2){
  s1 <- subset(trydata,time==time1)$individual
  s2 <- subset(trydata,time==time2)$individual
  a <- intersect(s1,s2)
  sub_codata1 <- subset(trydata,individual%in%a&time==time1)
  sub_codata2 <- subset(trydata,individual%in%a&time==time2)
  expectmean1 <- sub_codata1[,1:n1]/rowSums(sub_codata1[,1:n1])*rowSums(sub_codata1[,(1+n1):(n1+n2)])
  expectmean2 <- sub_codata2[,1:n1]/rowSums(sub_codata2[,1:n1])*rowSums(sub_codata2[,(1+n1):(n1+n2)])
  
  subv <- t(as.matrix(sub_codata1[,(1+n1):(n1+n2)]-expectmean1))%*%as.matrix(sub_codata2[,(1+n1):(n1+n2)]-expectmean1)/length(a)
  return(subv)
}


estimate_cov <- function(B_old,X,y,individual_time,subject){
  trydata <- data.frame(exp(X%*%B_old),y,individual_time,X)
  n1<- ncol(exp(X%*%B_old))
  n2 <- ncol(y)
  n3 <- ncol(X)
  sub_trydata <- subset(trydata,individual==subject)
  time <- sub_trydata$time
  es_cov_time1<- matrix(nrow=n1*length(time),ncol=0)

  for(time1 in time){
    es_cov_time2<- matrix(nrow=0,ncol=n1)
    for(time2 in time){
      es_cov_time2 <- rbind(es_cov_time2,sub_estimate_cov(trydata,subject,time1,time2,n1,n2))
    }
    
    es_cov_time1<- cbind(es_cov_time1,es_cov_time2)
  }
  return(ginv(es_cov_time1))
  
  
  
}



Grad_L<- function(B_old,X,y,individual_time){
  trydata <- data.frame(exp(X%*%B_old),y,individual_time,X)
  n1<- ncol(exp(X%*%B_old))
  n2 <- ncol(y)
  n3 <- ncol(X)
  Grad_l<- rep(0,n1*n3)
  hessian <- rep(0,n1*n3)
  for(subject in unique(trydata$individual)){
    subjectdata <- subset(trydata,individual==subject)
    pupbeta <- Pupbeta(B_old,X,y,individual_time,subject)
    inversecov <- estimate_cov(B_old,X,y,individual_time,subject)
    expecti <-  subjectdata[,1:n1]/rowSums(subjectdata[,1:n1])*rowSums( subjectdata[,(1+n1):(n1+n2)])
    yiminiui <- matrix(t(subjectdata[,(1+n1):(n1+n2)]-expecti),ncol=1)
    Grad_l<- Grad_l+pupbeta%*%inversecov%*%yiminiui
    Qij<- QIJ(B_old,X,y,individual_time,subject)
    hessian <-hessian- diag(pupbeta%*%inversecov%*%t(pupbeta))+as.vector(t(yiminiui)%*%inversecov%*%t(Qij))
    
    
  }
  hessian[which(hessian>-0.01)]<- -0.01
  Hs <- diag(hessian,nrow=p*n,ncol = p*n)
  ghlist <- list(Grad_l,hessian)
  return(ghlist)
}





#########calculate the hessian of the beta###########################
#######################################################################








##########calculate the initial beta###########################
####################################################################

B_ini0 <- function(X,y){
  B_0 <- matrix(nrow = ncol(X),ncol=ncol(y))
  logy <- log(y+1)
  for (j in 1 :12){
    B_0[,j] <- solve(t(X)%*%X)%*%t(X)%*%logy[,j]
  }
  return( B_0)
}
##############  Update V #######################
update_V <- function(B,Lambda,gamma2,theta,tau,index){
  m <- nrow(index)
  p <- nrow(B)
  V <- matrix(nrow=p,ncol=m)
  sigma <- gamma2/theta
  mcp <- tau*theta/(tau*theta -1)
  for(l in 1:m){
    i <- index[l,1]
    j <- index[l,2]
    x <- B[,i]-B[,j]-(1.0/(theta))*Lambda[,l]
    if(norm(x,type = "2")>=tau*gamma2){
      V[,l] <- x
    }
    else{
      V[,l] <- prox_L2(x,sigma)*mcp
    }
  }
  return(V)
}


##############  Update Lambda #######################

update_Lambda<- function(Lambda,B,V,theta,index){
  m <- nrow(index)
  p <- nrow(B) 
  Lambda_new <- matrix(nrow=p,ncol=m)
  for(l in 1:m){
    i <- index[l,1]
    j <- index[l,2]
    x <- V[,l]-B[,i]+B[,j]
    Lambda_new[,l] <- Lambda[,l] + theta*x
  }
  return(Lambda_new)
}


###########   Stopping criterion related #############
tolerance_primal<- function(B,V,eps_abs,eps_rel,index){
  m <- nrow(index)
  p <- nrow(V) 
  db <- 0.0
  dv <- 0.0
  output <- sqrt(p*m)*(eps_abs)
  for(l in 1:m){
    i <- index[l,1]
    j <- index[l,2]
    dbl <- (B[,i]-B[,j])^2
    db <- db + sum(dbl)
    dvl <- V[,l]^2
    dv <- dv +sum(dvl)
  } 
  db <- sqrt(db)
  dv <- sqrt(dv)
  output <- output +(eps_rel)*max(db,dv)
  return(output)
}


tolerance_dual <- function(Lambda,index,n,eps_abs,eps_rel){
  m <- nrow(index)
  p <- nrow(Lambda)
  output <- sqrt(p*n)*(eps_abs)
  AtLambda <- rep(0,p*n)
  eyemat_p <- diag(p)
  for(l in 1:m){
    i <- index[l,1]
    j <- index[l,2]
    e1 <- rep(0,n)
    e2 <- rep(0,n)
    e1[i] <- 1
    e2[j] <- 1
    Alt <- kronecker(e1-e2,eyemat_p)
    AtLambda = AtLambda + Alt%*%Lambda[,l]
  }
  s <- norm(AtLambda,type = "2")
  output <- output + s*(eps_rel)
  return(output)
}

########## This function computes the L2 norm of the primal residual.###########

residual_primal <- function(B,V,index){
  m <- nrow(index)
  p <- nrow(B)
  r <- 0.0
  for(l in 1:m){
    i <- index[l,1]
    j <- index[l,2]
    x <- (V[,l]-B[,i]+B[,j])^2
    r <- r +sum(x)
  }
  residual <- sqrt(r)
  return(residual)
}


##################This function the L2 norm of the dual residual############
residual_dual <- function(V,V_old,index,n,theta){
  m <- nrow(index)
  p <- nrow(V)
  s <- 0.0
  diff <- V - V_old
  for(i in 1 :n){
    v1 <- rep(0,p)
    v2 <- rep(0,p)
    for(ix in 1 :m){
      l1 <- index[ix,1]
      l2 <- index[ix,2]
      if(l1==i){
        v1 <- v1 + diff[,ix]
      }
      else{
        v1 <- v1
      }
      if(l2==i){
        v2 <- v2 + diff[,ix]
      }
      else{
        v2 <- v2
      }
    }
    v <- (v1-v2)^2
    s <- s + sum(v)
  }
  
  residual <- theta*sqrt(s)
  return(residual)
}

#######ADMM###########
prclust_admm <- function(X,y,diagD,B_0,index,gamma1,gamma2,theta,tau,n,p,max_iter,eps_abs,eps_rel,individual_time){
  m <- nrow(index)
  V <- zeros(p,m)
  Lambda<- zeros(p,m)
  V <- update_V(B=B_0, Lambda=Lambda, gamma2=gamma2, theta=theta, tau=tau, index=index)
  B_new <- B_0
  iter <- 1 
  for (its in 1:max_iter){
    V_old <- V
    Lambda_old <- Lambda
    B_old <- B_new
    B_new <- update_B(X=X, diagD=diagD, y=y, V=V_old, Lambda=Lambda_old, n=n, gamma1=gamma1, index=index, theta=theta,Hs=Hs,B_old=B_new,individual_time)
    V <- update_V(B=B_new, Lambda=Lambda_old, gamma2=gamma2, theta=theta, tau=tau, index=index)
    Lambda <- update_Lambda(Lambda=Lambda_old,B=B_new,V=V,theta=theta,index=index)
    rp <- residual_primal(B=B_new,V= V, index=index)
    #tp = tolerance_primal(B=B_new,V= V, eps_abs=eps_abs,eps_rel= eps_rel, index=index)
    rd = residual_dual(V=V, V_old=V_old, index=index, n=n, theta=theta)
    #td = tolerance_dual(Lambda=Lambda, index=index, n=n, eps_abs=eps_abs, eps_rel=eps_rel)
    if ((rp <= eps_abs) & (rd <= eps_rel)){
      break
    }
    iter <- iter +1
  }
  if (iter>max_iter){
    iter <- max_iter
  }
  finalresult <- list(B=B_new,Lambda=Lambda,V=V,rd=rd,rp=rp,iterations=iter)
  return( finalresult)
  
}




