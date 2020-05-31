rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

doubletry_simulate_Rtrace_ij1_ij2 <- function(subject,j1,j2,data,B_0,nknots,order,simulate_data_new){
  a1 <- subset(data,Capture.Number==j1)$individual
  a2 <- subset(data,Capture.Number==j2)$individual
  time1 <- subset(data,Capture.Number==j1&individual==subject)$time
  time2 <- subset(data,Capture.Number==j2&individual==subject)$time
  y10 <-subset(data,Capture.Number==j1&individual==subject)$total_n
  y20 <-subset(data,Capture.Number==j2&individual==subject)$total_n
  a <- intersect(a1,a2)
  sub_bafflodata1 <- subset(data,individual%in%a&Capture.Number==j1)
  sub_bafflodata2 <- subset(data,individual%in%a&Capture.Number==j2)
  y1 <- as.matrix(sub_bafflodata1[,2:13]/sub_bafflodata1$total_n)
  
  dmData <- matrix(0, nrow(y1), ncol(y1))
  for(i in 1:nrow(y1)){
    dmData[i,] <- stats::rmultinom(1, y10, y1[i,])
  }
  y2 <- as.matrix(sub_bafflodata2[,2:13]/sub_bafflodata1$total_n)
  dmData2 <- matrix(0, nrow(y2), ncol(y2))
  for(i in 1:nrow(y1)){
    dmData2[i,] <- stats::rmultinom(1, y20, y2[i,])
  }
  B1 <- bsplineS(time1,knots_eq3(simulate_data_new$time, k = order, m = nknots), norder = order)
  B2 <- bsplineS(time2,knots_eq3(simulate_data_new$time, k = order, m = nknots), norder = order)
  getalpha1 <-exp(B1%*%B_0)
  getalpha2 <- exp(B2%*%B_0)
  result1 <- digamma(rep.row(getalpha1,nrow(dmData))+dmData)
  result2 <- digamma(rep.row(getalpha2,nrow(dmData2))+dmData2)
  Bsum <- sum(as.vector(B1)*as.vector(B2))
  if(j1==j2){
    final_result <- sum(diag(var(result1))*as.vector(getalpha1)*as.vector(getalpha2))*Bsum
  }
  else{
    final_result <- sum(diag(cov(result1,result2))*as.vector(getalpha1)*as.vector(getalpha2))*Bsum
  }
  return(final_result)
}


doubleR0get <- function(simulate_data_new,B_0,nknots,order){
  
  nl <- nrow(simulate_data_new)
  R0 <- matrix(0,nl,nl)
  for (indexi in 1:nl){
    
    for(indexj in 1:nl){
      IDi <- simulate_data_new[indexi,]$individual
      IDj <- simulate_data_new[indexj,]$individual
      capturei <- simulate_data_new[indexi,]$Capture.Number
      capturej <- simulate_data_new[indexj,]$Capture.Number
      if(IDi==IDj){
        R0[indexi,indexj]<- doubletry_simulate_Rtrace_ij1_ij2(IDi,capturei,capturej,simulate_data_new,B_0,nknots,order,simulate_data_new)
      }
      
      else{
        R0[indexi,indexj] <- 0
      }
      
    }
  }
  print("finish")
  return(R0)
}
