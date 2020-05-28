######This file is aimed at investigating whether our simulated weights 
####is close to the theoretical weights
#####we only need to compare the cov(digamma(y_imk+a_imk),digamma(y_ink+a_ink)).







###################################################
################function of generate the data #####
###################################################


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

generate_data <-function(indiv_n,time,miss_n,b,ranges,curveparamenter){
  simulate_data_new<- list()
  for (i in 1:indiv_n){
    Get_alpha <- curveparamenter+matrix(b[i],ncol=ncol(curveparamenter),nrow=nrow(curveparamenter))
    simulate_data <- NULL
    for (j in 1:length(time)) {
      if(length(ranges)==1){
        total_n=ranges
      }
      else{
        total_n <- sample(ranges,1)
      }
      dm <- Dirichlet.multinomial(total_n,as.numeric(Get_alpha[j,]))
      dm_new <- c(dm,total_n,i)
      simulate_data <- rbind(simulate_data,dm_new )
    } 
    # simulate_data_new[[i]] <- sample_n(data.frame(time,simulate_data),size=miss_n[i])
    simulate_data_new[[i]] <- data.frame(time,simulate_data)
    rownames(simulate_data_new[[i]]) <- NULL
    
  }
  return(simulate_data_new)
}


combine_data <- function(generate.data){
  simulate_data <- NULL
  for (i in 1:indiv_n){
    simulate_data<- rbind(simulate_data,generate.data[[i]])
  }
  return(simulate_data)
}

Dirichlet.multinomial <- function(Nrs, shape){
  if(missing(Nrs) || missing(shape))
    stop("Nrs and/or shape missing.")
  
  # Create the data from the rmultinom
  dmData <- matrix(0, length(Nrs), length(shape))
  for(i in 1:length(Nrs))	
    dmData[i,] <- stats::rmultinom(1, Nrs[i], dirmult::rdirichlet(1, shape))
  
  # Label the created data
  colnames(dmData) <- paste("Taxa", 1:ncol(dmData))
  
  return(dmData)
}

generate.data_function <- function(time,indiv_n,b,ranges){
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
  
  generate.data<- generate_data(indiv_n,time,miss_n,b,ranges,curveparamenter)
  simulate_data <- combine_data(generate.data)
  names(simulate_data) <- c("time","Taxa.1","Taxa.2","Taxa.3","Taxa.4","Taxa.5","Taxa.6",
                            "Taxa.7","Taxa.8","Taxa.9","Taxa.10","Taxa.11","Taxa.12","total_n","individual")
  return(simulate_data)
  
}

time <- c(0,1)
indiv_n<-100
miss_n <- sample(10:15,indiv_n,replace = TRUE)
b <- rexp(indiv_n, rate = 0.001)

simulate_data_new <- generate.data_function(time,indiv_n,b,20000)
simulate_data_new$Capture.Number <- as.factor(simulate_data_new$time)
simulate_data_new$individual<- as.factor(simulate_data_new$individual)
data <- simulate_data_new

###################################################
##########finish function of generate the data ####
###################################################












################new version of calculate weights#################

wights_trysimulate_var_cov_ij1_ij2 <- function(j1,j2,data,k,subject){
  a1 <- subset(data,Capture.Number==j1)$individual
  a2 <- subset(data,Capture.Number==j2)$individual
  a <- intersect(a1,a2)
  time1 <- subset(data,Capture.Number==j1&individual==subject)$time
  time2 <- subset(data,Capture.Number==j2&individual==subject)$time
  y10 <- subset(data,Capture.Number==j1&individual==subject)$total_n
  y20 <- subset(data,Capture.Number==j2&individual==subject)$total_n
  getalpha1 <- c(rep(exp(5*cos(2*pi*time1)),4),rep(exp(1+exp(1.5*time1)),4),rep(exp( exp(1.8*(1-time1))),4))+rep(b[subject],12)
  getalpha2 <- c(rep(exp(5*cos(2*pi*time2)),4),rep(exp(1+exp(1.5*time2)),4),rep(exp( exp(1.8*(1-time2))),4))+rep(b[subject],12)
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
  cov(digamma(dmData[,k-1]+getalpha1[k-1]),digamma(dmData2[,k-1]+getalpha2[k-1]))
}

trysimulate_var_cov_matrix <- function(data,k){
  #nl <- nrow(data)
  nl<-10
  R0 <- matrix(0,nrow=nl,ncol=nl)
  for (indexi in 1:nl){
    for (indexj in 1:nl){
      
      IDi <- data[indexi,]$individual 
      IDj <- data[indexj,]$individual
      capturei <- data[indexi,]$Capture.Number
      capturej <- data[indexj,]$Capture.Number
      if(IDi==IDj){
        R0[indexi,indexj] <- wights_trysimulate_var_cov_ij1_ij2 (capturei,capturej,data,k,IDi)
      }
      else{
        R0[indexi,indexj] <- 0
      }
    }
  }
  return(R0) 
}




wights_theory_simulate_cov_ij1_ij2 <- function(j1,j2,data,k,subject){
  y10 <- subset(data,Capture.Number==j1&individual==subject)$total_n
  y20 <- subset(data,Capture.Number==j2&individual==subject)$total_n
  time1 <- subset(data,Capture.Number==j1&individual==subject)$time
  time2 <- subset(data,Capture.Number==j2&individual==subject)$time
  getalpha1 <- c(rep(exp(5*cos(2*pi*time1)),4),rep(exp(1+exp(1.5*time1)),4),rep(exp( exp(1.8*(1-time1))),4))+rep(b[subject],12)
  getalpha2 <- c(rep(exp(5*cos(2*pi*time2)),4),rep(exp(1+exp(1.5*time2)),4),rep(exp( exp(1.8*(1-time2))),4))+rep(b[subject],12)
  indiv_n<-1000
  curveparamenter1 <- t(data.frame(c(rep(exp(5*cos(2*pi*time1)),4),rep(exp(1+exp(1.5*time1)),4),rep(exp( exp(1.8*(1-time1))),4))))
  curveparamenter2 <- t(data.frame(c(rep(exp(5*cos(2*pi*time2)),4),rep(exp(1+exp(1.5*time2)),4),rep(exp( exp(1.8*(1-time2))),4))))
  miss_n <- sample(10:15,indiv_n,replace = TRUE)
  b <- rexp(indiv_n, rate = 0.001)
  dmData <- combine_data(generate_data(indiv_n,time1,miss_n,b,y10,curveparamenter1))
  dmData2 <- combine_data(generate_data(indiv_n,time2,miss_n,b,y20,curveparamenter2))
  cov(digamma(dmData[,k]+getalpha1[k-1]),digamma(dmData2[,k]+getalpha2[k-1]))
  
}



wights_theory_simulate_var_ij1_ij2 <- function(j,data,k,subject){
  j = capturei
  k=2
  subject=IDi
  y10 <- subset(data,Capture.Number==j&individual==subject)$total_n
  time1 <- subset(data,Capture.Number==j&individual==subject)$time
  getalpha1 <- c(rep(exp(5*cos(2*pi*time1)),4),rep(exp(1+exp(1.5*time1)),4),rep(exp( exp(1.8*(1-time1))),4))+rep(b[subject],12)
  curveparamenter1 <- t(data.frame(c(rep(exp(5*cos(2*pi*time1)),4),rep(exp(1+exp(1.5*time1)),4),rep(exp( exp(1.8*(1-time1))),4))))
  indiv_n<-1000
  miss_n <- sample(10:15,indiv_n,replace = TRUE)
  b <- rexp(indiv_n, rate = 0.001)
 
  dmData <- combine_data(generate_data(indiv_n,time1,miss_n ,b,y10,curveparamenter1))
  var(digamma(dmData[,k]+getalpha1[k-1]))
  
}

wights_empirical_var_cov_ij1_ij2 <- function(data,k){
  nl<-10
  R0 <- matrix(0,nrow=nl,ncol=nl)
  
  for (indexi in 1:nl){
    for (indexj in 1:nl){

      IDi <- data[indexi,]$individual 
      IDj <- data[indexj,]$individual
      capturei <- data[indexi,]$Capture.Number
      capturej <- data[indexj,]$Capture.Number
      if(IDi==IDj){

        aaaa <- wights_theory_simulate_cov_ij1_ij2(capturei,capturej,data,k,IDi)
        bbbb<- wights_theory_simulate_var_ij1_ij2(capturei,data,k,IDi)
        if(capturei==capturej){
          R0[indexi,indexj] <- bbbb}
        else{
          R0[indexi,indexj] <- aaaa}
      }
      
      
      else{
        R0[indexi,indexj] <- 0
      }
    }
  }
  return(R0)
}
####simulation value of weights
trysimulate_var_cov_matrix(data,2)
#####theoretical value of weights
wights_empirical_var_cov_ij1_ij2(data,2)
