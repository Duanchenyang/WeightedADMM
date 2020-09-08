library(dirmult)
library(plyr)
library(fda)#bspline basis
library(Matrix)
library(matrixcalc)
library(igraph)
library(tidyverse)
library(matlab)
library(magic)


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

ranges <- 5000:20000
nknots <- 3
order <- 3
n <- 12
gamma1 <- 1
gamma2 <- 3000
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
miss_n <- seq(0.6,1,by=0.1)
generate.data<- generate_data(indiv_n,time,miss_n,b,ranges,curveparamenter)
simulate_data_new <- combine_data(generate.data)
names(simulate_data_new) <- c("time","Taxa.1","Taxa.2","Taxa.3","Taxa.4","Taxa.5","Taxa.6",
                              "Taxa.7","Taxa.8","Taxa.9","Taxa.10","Taxa.11","Taxa.12","total_n","individual")

simulate_data_new$Capture.Number <- simulate_data_new$time*(length(unique(simulate_data_new$time))-1)+1
simulate_data_new$individual<- as.factor(simulate_data_new$individual)
data <- simulate_data_new

####finish generate data
y <- as.matrix(simulate_data_new[,2:13])
X <- bsplineS(simulate_data_new$time,knots_eq3(simulate_data_new$time, k = order, m = nknots), norder = order)
individual_time <- simulate_data_new[,c(1,14,15,16)]
B_0 <- B_ini0(X,y)


try.sample <- prclust_admm(X=X, y=y, diagD=diagD, B_0=B_0, index=index,
                           gamma1 = gamma1, gamma2 = gamma2, theta=2000, tau = 2, n=n, p=p,  max_iter=10000,
                           eps_abs=1e-3, eps_rel=1e-3,individual_time=individual_time)

Ad_final <- create_adjacency(try.sample$V, n)
G_final <- graph.adjacency(Ad_final, mode = 'upper')
#clustering membership
cls_final <- components(G_final)
#number of clusters
k_final <- cls_final$no
cls_final










#####generate the plot
finalresult <- cls_final
B_ini1 <- function(X,y){
  B_0 <- matrix(nrow = ncol(X),ncol=ncol(y))
  logy <- log(y+1)
  for(i in 1 :length(finalresult$csize)){
    choose_class <-unique(finalresult$membership)[i]
    sameg <- which(finalresult$membership==choose_class)
    X_sameg<- NULL
    logy_sameg<-NULL
    for(j in 1:length(sameg)){
      X_sameg<- rbind(X_sameg,X)
      logy_sameg<- c(logy_sameg,logy[,sameg[j]])
    }
    B_0solution <- solve(t(X_sameg)%*%X_sameg)%*%t(X_sameg)%*%logy_sameg
    B_0[,sameg]<- matrix(rep(B_0solution,length(sameg)),nrow=ncol(X))
  }
  
  return( B_0)
}

B_coefficient<- try.sample$B
B_coefficient2<- B_ini1(X,y)

##########generate the corresponding plot


B00 <- bsplineS(simulate_data_new$time ,knots_eq3(simulate_data_new$time, k = order, m = nknots), norder = order)

estimate_data_weights <- gather(data.frame(exp(B00%*%B_coefficient)/rowSums(exp(B00%*%B_coefficient))),X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,key="par",value="value_weights")
estimate_data_regular <- gather(data.frame(exp(B00%*%B_coefficient2)/rowSums(exp(B00%*%B_coefficient2))),X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,key="par",value="value_regular")
true_data <- gather(simulate_data_new[1:14],Taxa.1,Taxa.2,Taxa.3,Taxa.4,Taxa.5,Taxa.6,Taxa.7,Taxa.8,Taxa.9,Taxa.10,Taxa.11,Taxa.12,key="par",value="true_value")

curve <- c(exp(5*cos(2*pi*simulate_data_new$time)),exp(5*cos(2*pi*simulate_data_new$time))
           ,exp(5*cos(2*pi*simulate_data_new$time)),exp(5*cos(2*pi*simulate_data_new$time))
           ,exp(1+exp(1.5*simulate_data_new$time))
           ,exp(1+exp(1.5*simulate_data_new$time))
           ,exp(1+exp(1.5*simulate_data_new$time))
           ,exp(1+exp(1.5*simulate_data_new$time))
           ,exp( exp(1.8*(1-simulate_data_new$time)))
           ,exp( exp(1.8*(1-simulate_data_new$time)))
           ,exp( exp(1.8*(1-simulate_data_new$time)))
           ,exp( exp(1.8*(1-simulate_data_new$time))))/(exp(5*cos(2*pi*simulate_data_new$time))+exp(1+exp(1.5*simulate_data_new$time))+exp( exp(1.8*(1-simulate_data_new$time))))/4
dataforplot <- data.frame(true_data[1:4],estimate_data_weights[2],estimate_data_regular[2],curve)




ggplot(dataforplot)+
  geom_point(aes(x=time,y=true_value/total_n),color="gray50",size=0.8)+
  geom_line(aes(x=time,y=value_regular),size=0.75)+
  geom_line(aes(x=time,y=value_weights),size=0.75,col="blue")+
  geom_line(aes(x=time,y=curve),size=0.75,col="red")+
  facet_wrap( ~ par , ncol=4)+
  labs(x='Age (month)',y='Relative abundance')
