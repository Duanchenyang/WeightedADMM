library(dirmult)
library(plyr)
library(fda)#bspline basis
library(Matrix)
library(matrixcalc)
library(igraph)
library(tidyverse)
library(matlab)

########Initial value##########
time <- seq(0,1,0.07)
indiv_n<-50
miss_n <- sample(10:15,indiv_n,replace = TRUE)
b <- rexp(indiv_n, rate = 0.001)
a1 <- exp(5*cos(2*pi*time))
a2 <- exp(5*cos(2*pi*time))
a3 <- exp(5*cos(2*pi*time))
a4 <- exp(5*cos(2*pi*time))
a5 <- exp(1+exp(1.5*time))
a6 <- exp(1+exp(1.5*time))
a7 <- exp(1+exp(1.5*time))
a8 <- exp(1+exp(1.5*time))
a9 <- exp( exp(1.8*(1-time)))
a10 <- exp( exp(1.8*(1-time)))
a11 <- exp( exp(1.8*(1-time)))
a12 <- exp( exp(1.8*(1-time)))
curveparamenter <- data.frame(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)
ranges <- 5000:20000
nknots <- 3
order <- 3
n <- 12
gamma1 <- 10
gamma2 <- 0.1
p <- order + nknots
C <- matrix(0, nrow=nknots+order-2, ncol=nknots+order)
for (j in 1:(nknots+order-2)){
  d_j <- c(rep(0,j-1),1,-2,1,rep(0,(nknots+order)-3-(j-1)))
  e_j <- c(rep(0,j-1), 1 ,rep(0,(nknots+order)-3-(j-1)))
  C <- C + e_j%*%t(d_j)
}
D = t(C)%*%C
diagD <- kronecker(diag(1, n), D)

index = t(combn(n,2))



evaluation <- function(k){ 
  source('function of calculate weights.R')
  source('Function of generating data.R')
  source('function of weighted ADMM.R')
  
  
  ####generate data
  generate.data<- generate_data(indiv_n,time,miss_n,b,ranges,curveparamenter)
  simulate_data_new <- combine_data(generate.data)
  names(simulate_data_new) <- c("time","Taxa.1","Taxa.2","Taxa.3","Taxa.4","Taxa.5","Taxa.6",
                                "Taxa.7","Taxa.8","Taxa.9","Taxa.10","Taxa.11","Taxa.12","total_n","individual")
  
  simulate_data_new$Capture.Number <- as.factor(simulate_data_new$time)
  simulate_data_new$individual<- as.factor(simulate_data_new$individual)
  data <- simulate_data_new
  ####finish generate data
  
  ####calculate the weights####
  
  
  y <- as.matrix(simulate_data_new[,2:13])
  X <- bsplineS(simulate_data_new$time,knots_eq3(simulate_data_new$time, k = order, m = nknots), norder = order)
  B_0 <- B_ini0(X,y)
  Rpoint <- doubleR0get(simulate_data_new,B_0,nknots,order)
  wij1 <- colSums(solve(Rpoint))/sum(colSums(solve(Rpoint)))*indiv_n*length(time)
  ####finish calculate the weights####
  
  ######Start calculating the parameter#####
  try.sample <- prclust_admm(X=X, y=y, diagD=diagD, B_0=B_0, index=index,
                             gamma1 = gamma1, gamma2 = gamma2, theta=1100, tau = 2, n=n, p=p,  max_iter=10000,
                             eps_abs=1e-3, eps_rel=1e-3,wij1)
  B_coefficient<- try.sample$B
  B1 <- data.frame(as.numeric(B_coefficient))
  wij2 <- rep(1,length(wij1))
  try.sample2 <- prclust_admm(X=X, y=y, diagD=diagD, B_0=B_0, index=index,
                              gamma1 = gamma1, gamma2 = gamma2, theta=1100, tau = 2, n=n, p=p,  max_iter=10000,
                              eps_abs=1e-3, eps_rel=1e-3,wij2)
  B_coefficient2<- try.sample2$B
  B2 <- data.frame(as.numeric(B_coefficient2))
  write.csv(B1,file= paste('weighted_ADMM',k,'.csv',sep=""))
  write.csv(B2,file= paste('regular_ADMM',k,'.csv',sep=""))
  
}

args = commandArgs(trailingOnly = TRUE)
k = as.numeric(args[1])
evaluation(k)


