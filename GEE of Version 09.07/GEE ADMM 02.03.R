#####################################################
##This is alpha function######
#####################################################

####################
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

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
######################



Dispersion_parameter1 <- function(B_old,X,y,individual_time){
  L <- ncol(X)
  n <- ncol(y)
  N<- nrow(individual_time)
  element_data <- data.frame(X,y,individual_time)
  alpha <- exp(as.matrix(element_data[,1:L])%*%B_old)
  numerator <- as.vector(t(y-alpha/rowSums(alpha)*rowSums(y)))^2
  denominator<- as.vector(t(alpha/rowSums(alpha)*rowSums(y)))*as.vector(t((1-alpha/rowSums(alpha))*(rowSums(y)+rowSums(alpha))/(1+rowSums(alpha))))
  return(sum(numerator/denominator)/(n*N-1))
}

Dispersion_parameter2<- function(B_old,X,y,individual_time){
  L <- ncol(X)
  n <- ncol(y)
  N<- nrow(individual_time)
  element_data <- data.frame(X,y,individual_time)
  alpha <- exp(as.matrix(element_data[,1:L])%*%B_old)
  numerator_matrix <- y-alpha/rowSums(alpha)*rowSums(y)
  denominator_matrix1<- alpha/rowSums(alpha)*rowSums(y)
  denominator_matrix2<- (1-alpha/rowSums(alpha))*(rowSums(y)+rowSums(alpha))/(1+rowSums(alpha))
  phk_0<- 0
  for(i in 1:n){
    phik <- sum(numerator_matrix[,i]^2/denominator_matrix1[,i]/denominator_matrix2[,i])/(N-1/n)
    phk_0<-phk_0+phik 
  }
  return(phk_0/n)
}



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

update_B <- function(diagD,V,Lambda,n,gamma1,index,theta,X,y,Hs,B_old,individual_time){
  m <- nrow(index)
  p <- nrow(V)
  B_new <- matrix(nrow=p,ncol=n)
  V_tilde <- V+Lambda/theta
  eyemat_n <- diag(n)
  eyemat_p <- diag(p)
  onemat_n <- ones(n)
  AtA <- kronecker(n*eyemat_n - onemat_n,eyemat_p)
  Hs <- Hessian(B_old,X,y,individual_time)
  #print(diag(Hs))
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
  Grad_L <- Grad_L(B_old,X,y,individual_time)
  rhs <- Grad_L-t(as.vector(B_old)%*%Hs)+theta*lv
  b_new <- invmat%*%rhs
  B_new <- matrix(b_new,nrow=p,ncol=n)
  return(B_new)
}

#####calculate the gradient of beta############################

element_need <- function(B_old,X,y,individual_time,subject){
  L <- ncol(X)
  n <- ncol(y)
  element_data <- data.frame(X,y,individual_time)
  sub_element_data <- subset(element_data,individual==subject)
  alpha <- exp(as.matrix(sub_element_data[,1:L])%*%B_old)
  t_pupbeta <- matrix(nrow=ncol(alpha)*L,ncol=0)
  qij <- matrix(nrow=ncol(alpha)*L,ncol=0)
  yminusu <- NULL
  pearson_residual <- NULL
  V<- matrix(nrow=0,ncol=0)
  R<- matrix(nrow=0,ncol=0)
  diag_person_resudual<-matrix(nrow=0,ncol=0)
  for(j in 1:nrow(sub_element_data)){
    a <- alpha[j,]/sum(alpha[j,])
    A<- diag(a)-matrix(a,ncol=1)%*%matrix(a,nrow=1)
    B<- matrix(unlist(sub_element_data[j,1:L]),ncol=1)
    yij0 <- sub_element_data$total_n[j]
    yiminusui <- matrix(unlist(sub_element_data[j,(L+1):(L+n)]))[,1]-a*yij0
    Identity_n <- diag(rep(1,n))
    C <- (Identity_n-2*diag(a))%*%A*yij0
    
    t_pupbeta <- cbind(t_pupbeta,kronecker(A,B)*yij0)
    yminusu <- c(yminusu,yiminusui)
    qij<- cbind(qij,kronecker(C,B^2))
    
    V_block <- yij0*(yij0+sum(alpha[j,]))/(1+sum(alpha[j,]))*A
    V <- adiag(V,V_block)
    Rblock <- diag(1/sqrt(diag(V_block)))%*%V_block%*%diag(1/sqrt(diag(V_block)))
    R <- adiag(R,Rblock)
    phi <- Dispersion_parameter1(B_old,X,y,individual_time)
    pearson_residual_j <- yiminusui/sqrt(phi*yij0*(yij0+sum(alpha[j,]))/(1+sum(alpha[j,]))*a*(1-a))
    pearson_residual <- c(pearson_residual,pearson_residual_j)
    
  }
  
  return(list(t_pupbeta,yminusu,qij,V,pearson_residual,R))
  
}



cal_alpha_beta <- function(B_old,X,y,individual_time){
  cora <- 0
  ni <- length(unique(individual_time$individual))
  K <- ncol(y)
  regressiondata<- NULL
  for (subject in unique(individual_time$individual)){
    #change alpha##
    capturenumber <- subset(individual_time,individual==subject)$Capture.Number
    numtime<- length(unique( capturenumber))
    j1minusj2 <- ones(numtime*K,numtime*K)
    
    element<- element_need(B_old,X,y,individual_time,subject)
    pearson_residual<- element[[5]]
    matrix_pearson_residual <- matrix(pearson_residual,ncol=1)%*%matrix(pearson_residual,nrow=1)
    matrix_pearson_residual[lower.tri(matrix_pearson_residual)]=0
    Rijkl <- as.vector( matrix_pearson_residual)
    newRijkl <- Rijkl[Rijkl!=0]
    abs_j1_j2 <- abs(as.vector( j1minusj2))[Rijkl!=0]
    fordata<- data.frame(newRijkl ,abs_j1_j2) 
    regressiondata<- rbind(regressiondata, fordata)
  }
  
  ###change model#####
  #print(sd(regressiondata$newRijkl))
  mod <-summary(nls(newRijkl~Alpha^abs_j1_j2,data=regressiondata,start = list(Alpha = 0.2),algorithm="port", lower=-1, upper=1))
  #mod <-summary(nls(newRijkl~(Alpha^abs_j1_j2)*((B2/Alpha)^indicator_beta),data=regressiondata,start = list(Alpha = 0.2, B2 = 0.5),algorithm="port", 
  #lower=c(-1,-0.99), upper=c(1,0.99) ))
  
  cor_Alpha <- mod$coefficients[1,1]
  return( list(cor_Alpha))
  
}





Grad_L<- function(B_old,X,y,individual_time){
  Grad_l <- rep(0,nrow(B_old)*ncol(B_old))
  ni <- length(unique(individual_time$individual))
  #ntime <- length(unique(individual_time$time))
  K <- ncol(y)
  
  result<- cal_alpha_beta(B_old,X,y,individual_time)
  cor_Alpha<-  result[[1]]
  print(cor_Alpha)
  #cor_Alpha<-0.130994
  
  for(subject in unique(individual_time$individual)){
    element<- element_need(B_old,X,y,individual_time,subject)
    t_pupbeta<- element[[1]]
    yminusu<-element[[2]]
    qij<-element[[3]]
    V<- element[[4]]
    R <- element[[6]]
    
    capturenumber <- subset(individual_time,individual==subject)$Capture.Number
    ####might need change#########
    corR <-kronecker(ones(length(capturenumber),length(capturenumber))-diag(rep(1,length(capturenumber))),ones(K,K))*cor_Alpha
    V_inv <- diag(1/sqrt(diag(V)))%*%ginv(R+corR)%*%diag(1/sqrt(diag(V)))
    Grad_l <-Grad_l+as.vector(t_pupbeta%*%V_inv%*%matrix(yminusu,ncol=1))
  }
  return(Grad_l)
}


Hessian <- function(B_old,X,y,individual_time){
  p <- nrow(B_old)
  n <- ncol(B_old)
  #ntime <- length(unique(individual_time$time))
  K <- ncol(y)
  Hessian_diag <- rep(0,n*p)
  result<- cal_alpha_beta(B_old,X,y,individual_time)
  cor_Alpha<-  result[[1]]
  #cor_Alpha<-0.130994
  
  for(subject in unique(individual_time$individual)){
    element <- element_need(B_old,X,y,individual_time,subject)
    t_pupbeta <- element[[1]]
    yminusu <-element[[2]]
    qij <-element[[3]]
    V <- element[[4]]
    R <- element[[6]]
    capturenumber <- subset(individual_time,individual==subject)$Capture.Number
    ####might need change#########
    corR <-kronecker(ones(length(capturenumber),length(capturenumber))-diag(rep(1,length(capturenumber))),ones(K,K))*cor_Alpha
    V_inv <-  diag(1/sqrt(diag(V)))%*%ginv(R+corR)%*%diag(1/sqrt(diag(V)))
    
    Hessian_diag<- Hessian_diag-diag(t_pupbeta%*%V_inv%*%t(t_pupbeta))+
      as.vector(matrix(yminusu,nrow=1)%*%V_inv%*%t(qij))
    
  }
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
    tp = tolerance_primal(B=B_new,V= V, eps_abs=eps_abs,eps_rel= eps_rel, index=index)
    rd = residual_dual(V=V, V_old=V_old, index=index, n=n, theta=theta)
    td = tolerance_dual(Lambda=Lambda, index=index, n=n, eps_abs=eps_abs, eps_rel=eps_rel)
    if ((rp <= tp) & (rd <= td)){
      break
    }
    iter <- iter +1
    print(iter)
  }
  if (iter>max_iter){
    iter <- max_iter
  }
  finalresult <- list(B=B_new,Lambda=Lambda,V=V,rd=rd,rp=rp,iterations=iter)
  return( finalresult)
}