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

