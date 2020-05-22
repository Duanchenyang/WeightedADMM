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

update_B <- function(diagD,V,Lambda,n,gamma1,index,theta,X,y,Hs,B_old){
  m <- nrow(index)
  p <- nrow(V)
  B_new <- matrix(nrow=p,ncol=n)
  V_tilde <- V+Lambda/theta
  eyemat_n <- diag(n)
  eyemat_p <- diag(p)
  onemat_n <- ones(n)
  AtA <- kronecker(n*eyemat_n - onemat_n,eyemat_p)
  Hs <- Hs(B_old,X,y)
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
  Grad_L <- Grad_L(B_old,X,y)
  rhs <- Grad_L-t(as.vector(B_old)%*%Hs)+theta*lv
  b_new <- invmat%*%rhs
  B_new <- matrix(b_new,nrow=p,ncol=n)
  return(B_new)
}

#####calculate the gradient of beta############################
########B_old is beta(p,n) p=6,n=12,X is B-spline(256*6),y is (256*12)#####

Grad_L <- function(B_old,X,y){
  Grad_l <- matrix(nrow=nrow(B_old),ncol = ncol(B_old))
  for(j in 1 :12){
    original1 <-(digamma(rowSums(exp(X%*%B_old)))-digamma(rowSums(y)+rowSums(exp(X%*%B_old))))*exp(X%*%B_old[,j])[,1]*X*wij
    indexmat1 <- rowSums(is.na((digamma(rowSums(exp(X%*%B_old)))-digamma(rowSums(y)+rowSums(exp(X%*%B_old))))*exp(X%*%B_old[,j])[,1]*X))
    original1[which(indexmat1>0),]<- (-rowSums(y)*X)[which(indexmat1>0),]
    original2 <- (digamma(y[,j]+exp(X%*%B_old[,j]))[,1]-digamma(exp(X%*%B_old[,j]))[,1])*exp(X%*%B_old[,j])[,1]*X*wij
    indexmat2 <- rowSums(is.na((digamma(y[,j]+exp(X%*%B_old[,j]))[,1]-digamma(exp(X%*%B_old[,j]))[,1])*exp(X%*%B_old[,j])[,1]*X))
    original2[which(indexmat2>0),]<- (y[,j]*X)[which(indexmat2>0),]
    Grad_l[,j]<- colSums(original1+original2)
  }
  return(as.vector(Grad_l))
}


#########calculate the hessian of the beta###########################
#######################################################################


Hs <-function(B_old,X,y){
  p <- nrow(B_old)
  n <- ncol(B_old)
  Hessian_diag <- matrix(nrow=p,ncol = n)
  for(j in 1:12){
    Hessian_diag[,j]<-colSums((trigamma(rowSums(exp(X%*%B_old)))-trigamma(rowSums(y)+rowSums(exp(X%*%B_old)))+trigamma(y[,j]+exp(X%*%B_old[,j]))[,1]-trigamma(exp(X%*%B_old[,j]))[,1]  )*((exp(X%*%B_old[,j]))^2)[,1]*(X^2)*wij
                              +(digamma(rowSums(exp(X%*%B_old)))-digamma(rowSums(y)+rowSums(exp(X%*%B_old)))+digamma(y[,j]+exp(X%*%B_old[,j]))[,1]-digamma(exp(X%*%B_old[,j]))[,1])*exp(X%*%B_old[,j])[,1]*(X^2)*wij)
  }
  Hessian_diag<- as.vector(Hessian_diag)
  Hessian_diag[which(Hessian_diag>-1)]<- -1
  Hs <- diag(Hessian_diag,nrow=p*n,ncol = p*n)
  return(Hs)
}


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
prclust_admm <- function(X,y,diagD,B_0,index,gamma1,gamma2,theta,tau,n,p,max_iter,eps_abs,eps_rel){
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
    B_new <- update_B(X=X, diagD=diagD, y=y, V=V_old, Lambda=Lambda_old, n=n, gamma1=gamma1, index=index, theta=theta,Hs=Hs,B_old=B_new)
    V <- update_V(B=B_new, Lambda=Lambda_old, gamma2=gamma2, theta=theta, tau=tau, index=index)
    Lambda <- update_Lambda(Lambda=Lambda_old,B=B_new,V=V,theta=theta,index=index)
    rp <- residual_primal(B=B_new,V= V, index=index)
    tp = tolerance_primal(B=B_new,V= V, eps_abs=eps_abs,eps_rel= eps_rel, index=index)
    rd = residual_dual(V=V, V_old=V_old, index=index, n=n, theta=theta)
    td = tolerance_dual(Lambda=Lambda, index=index, n=n, eps_abs=eps_abs, eps_rel=eps_rel)
    if ((rp <= tp) & (rd <= td)){
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




