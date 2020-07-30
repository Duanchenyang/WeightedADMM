
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
  #geom_line(aes(x=time,y=value_regular),size=0.75)+
  geom_line(aes(x=time,y=value_weights),size=0.75,col="blue")+
  geom_line(aes(x=time,y=curve),size=0.75,col="red")+
  facet_wrap( ~ par , ncol=4)+
  labs(x='Age (month)',y='Relative abundance')


