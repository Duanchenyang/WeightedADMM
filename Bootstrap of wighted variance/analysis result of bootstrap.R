
#####calculate the variance of relative abundance of  weighted data 

varweightq <- seq(0,0.98,0.01)

for(q in 1:length(seq(0,0.98,0.01))){

  varianceweight<-rep(1,500)
  for (i in 1:500){
    print(i)
    try_csv <- try(read.csv(paste('~/Desktop/data saver/weighted_ADMM',i,'.csv',sep="")),silent = TRUE)
    try_answer<-'try-error' %in% class(try_csv)
    if (try_answer==FALSE){
      weighted_ADMM <- read.csv(paste('~/Desktop/data saver/weighted_ADMM',i,'.csv',sep=""))
      weighted_ADMM1 <- matrix(weighted_ADMM[,2],ncol=12)
      B00 <- bsplineS(seq(0,0.98,0.01)[q] ,knots_eq3(simulate_data_new$time, k = order, m = nknots), norder = order)
      varianceweight[i]<- (B00%*%weighted_ADMM1)[1,5]
    }
    else{varianceweight[i]<-NA}
    
  }
  varweightq[q]<- var(varianceweight,na.rm = TRUE)
}

#####calculate the variance of relative abundance of  unweighted data 
varregularq <- seq(0,0.98,0.01)
for(q in 1:length(seq(0,0.98,0.01))){

  varianceregular<-rep(1,500)

for (i in 1:500){
  print(i)
  try_csv <- try(read.csv(paste('~/Desktop/data saver/regular_ADMM',i,'.csv',sep="")),silent = TRUE)
  try_answer<-'try-error' %in% class(try_csv)
  if (try_answer==FALSE){
    weighted_ADMM <- read.csv(paste('~/Desktop/data saver/regular_ADMM',i,'.csv',sep=""))
    weighted_ADMM1 <- matrix(weighted_ADMM[,2],ncol=12)
    B00 <- bsplineS(seq(0,0.98,0.01)[q],knots_eq3(simulate_data_new$time, k = order, m = nknots), norder = order)
    varianceregular[i]<- (B00%*%weighted_ADMM1)[1,5]
  }
  else{varianceregular[i]<-NA}
  
}
  varregularq[q]<- var(varianceregular,na.rm = TRUE)
}

######plot#######
varq <- seq(0,0.98,0.01)
data.frame(varq,varweightq,varregularq)%>%
ggplot()+
  geom_point(aes(x=varq,y=varweightq),col="blue")+
  geom_point(aes(x=varq,y=varregularq))

####compare beta variance#######
beta_matrix_wighted <- matrix(nrow=72)
beta_matrix_regular <- matrix(nrow=72)

for (i in 1:500){
  try_csv <- try(read.csv(paste('~/Desktop/data saver/regular_ADMM',i,'.csv',sep="")),silent = TRUE)
  try_answer<-'try-error' %in% class(try_csv)
  if (try_answer==FALSE){
    regular_ADMM <- read.csv(paste('~/Desktop/data saver/regular_ADMM',i,'.csv',sep=""))
    weighted_ADMM <- read.csv(paste('~/Desktop/data saver/weighted_ADMM',i,'.csv',sep=""))
    beta_matrix_wighted<-cbind(beta_matrix_wighted,matrix(weighted_ADMM[,2]))
    beta_matrix_regular<- cbind(beta_matrix_regular,matrix(regular_ADMM[,2]))
  }
}
###variance of weighted beta
apply(beta_matrix_wighted[,c(-1,-70)],1,var)
###variance of unweighted beta
apply(beta_matrix_regular[,-1],1,var)


