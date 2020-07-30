library(dirmult)
library(plyr)
library(fda)#bspline basis
library(Matrix)
library(matrixcalc)
library(igraph)
library(tidyverse)
library(matlab)

time <- seq(0,1,0.05)
indiv_n<-10
miss_n <- sample(10:15,indiv_n,replace = TRUE)
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
curveparamenterold <- data.frame(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)
atotal <- rowSums(curveparamenterold)
curveparamenter<- data.frame(matrix(0,ncol=12,nrow=length(a1)))

for (i in 1:length(a1)){
  curveparamenter[i,] <- as.vector(curveparamenterold[i,]/atotal[i]*50)
}

ranges <- 20000
nknots <- 3
order <- 3
n <- 12
gamma1 <- 0.5
gamma2 <- 12
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

b <- rexp(indiv_n, rate = 0.7)
generate.data<- generate_data(indiv_n,time,miss_n,b,ranges,curveparamenter)
simulate_data_new <- combine_data(generate.data)
names(simulate_data_new) <- c("time","Taxa.1","Taxa.2","Taxa.3","Taxa.4","Taxa.5","Taxa.6",
                              "Taxa.7","Taxa.8","Taxa.9","Taxa.10","Taxa.11","Taxa.12","total_n","individual")

simulate_data_new$Capture.Number <- as.factor(simulate_data_new$time)
simulate_data_new$individual<- as.factor(simulate_data_new$individual)
data <- simulate_data_new

####finish generate data
y <- as.matrix(simulate_data_new[,2:13])
X <- bsplineS(simulate_data_new$time,knots_eq3(simulate_data_new$time, k = order, m = nknots), norder = order)
individual_time <- simulate_data_new[,c(1,15)]
B_0 <- B_ini0(X,y)
matrix(sample(5:7,72,replace=TRUE),ncol=12)


try.sample <- prclust_admm(X=X, y=y, diagD=diagD, B_0=B_0, index=index,
                           gamma1 = gamma1, gamma2 = gamma2, theta=20, tau = 2, n=n, p=p,  max_iter=1000,
                           eps_abs=1e-3, eps_rel=1e-3,individual_time=individual_time)
Ad_final <- create_adjacency(try.sample$V, n);
G_final <- graph.adjacency(Ad_final, mode = 'upper')
#clustering membership
cls_final <- components(G_final)
#number of clusters
k_final <- cls_final$no

