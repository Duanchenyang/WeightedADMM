library(plyr)
library(fda)#bspline basis
library(Matrix)
library(matrixcalc)
library(igraph)
library(Rcpp)
library(tidyverse)
library(matlab)
data <-read.csv('GEE of Version 09.07/juvenile_growth.csv', header = TRUE)
timerange <- data$age.years*12
y <- as.matrix(round(data[,12:23]*data[,2],0))
X <- bsplineS(timerange,knots_eq3(timerange, k = order, m = nknots), norder = order)
individual <- data$Numeric.Animal.ID
time <- data$age.years*12
total_n<- data$read.number
Capture.Number<- data$Capture.Number
individual_time <- data.frame(time,total_n,individual,Capture.Number)



nknots <- 3
order <- 3
n <- 12
gamma1 <- 1
gamma2 <- 30000
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
B_0 <- B_ini0(X,y)


try.sample <- prclust_admm(X=X, y=y, diagD=diagD, B_0=B_0, index=index,
                           gamma1 = gamma1, gamma2 = gamma2, theta=20000, tau = 2, n=n, p=p,  max_iter=10000,
                           eps_abs=1e-3, eps_rel=1e-3,individual_time=individual_time)

Ad_final <- create_adjacency(try.sample$V, n)
G_final <- graph.adjacency(Ad_final, mode = 'upper')
#clustering membership
cls_final <- components(G_final)
#number of clusters
k_final <- cls_final$no
cls_final









#####generate the plot############

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



B_coeffient <- try.sample$B
B_coeffient<- B_ini1(X,y)
trydata <- data.frame(exp(X%*%B_coeffient)/rowSums(exp(X%*%B_coeffient)),timerange,y,data[,2])

trydata1<- gather(trydata,Tvelifera.sc,TveliferaB,TveliferaUD,Tmutans1,Tmutans2,Tmutans3,
                  TmutansMSD,Tmutans.sc, TmutansUD,TspBg,TspBf,Tparva,key="type2",value="true_response2")
trydata1$type2 <- factor(trydata1$type2,level=c("Tvelifera.sc", "TveliferaB" ,       
                                                "TveliferaUD", "Tmutans1",   
                                                "Tmutans2"  , "Tmutans3"  ,           
                                                "TmutansMSD"    , "Tmutans.sc"   ,     
                                                "TmutansUD"       ,     "TspBg" ,       
                                                "TspBf"   , "Tparva"))

trydata2 <- gather(trydata,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,key="type",value="true_response")

trydata2$type <- factor(trydata2$type,level=c("X1","X2","X3","X4","X5","X6",
                                              "X7","X8","X9","X10","X11","X12"))
try.labs<-levels(trydata1$type2)

names(try.labs)<- c("X1","X2","X3","X4","X5","X6",
                    "X7","X8","X9","X10","X11","X12")

data.new <- data.frame(trydata2[,c("type","timerange","true_response")],trydata1[,c("true_response2","type2")],data[,2])
data.new$group <- rep(0,length(data.new$type))
for(i in 1:12){
  data.new$group[data.new$type==paste("X",i,sep = "")] <- cls_final$membership[i]}
data.new$group<- as.factor(data.new$group)
ggplot(data.new)+
  geom_point(aes(x=timerange,y=true_response2/`data...2.`),color="gray50",size=0.8)+
  geom_line(aes(x=timerange,y=true_response),size=0.75)+
  facet_wrap( ~ type, ncol=4,labeller = as_labeller(try.labs))+
  labs(x='Age (month)',y='Relative abundance')



