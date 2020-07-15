rm(list=ls())
library(plyr)
library(fitR)
###本程序可同时适用发病人数(本例)和发病人数比例两种情形，互换时需对line 11的SEIR_ode、line 144、146、154进行修改
##########PTR1(time-dependent transmission rate)
data(SEITL_deter)
PTR1 <- SEITL_deter
###################################################
PTR1$simulate <- function(init.theta,init.state,times) {
  
  SEIR_ode <- function(time, state, theta) {
    
    # param
    beta <- theta[[1]][["beta"]]
    epsilon <- 1/theta[[1]][["D_lat"]]
    nu <- 1/theta[[1]][["D_inf"]]
    
    # states
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    R <- state[["R"]]
    Inc <- state[["Inc"]]
    
    N <- S + E +I + R
    #N <- S + I + R
    
    dS <- -beta*S*I/N 
    #dS <- -beta*S*I
    dE <- beta*S*I/N - epsilon*E
    #dE <- beta*S*I - epsilon*E
    dI <- epsilon*E - nu*I
    #dI <- beta*S*I - nu*I
    dR <- nu*I
    dInc <- epsilon*E
    #dInc <- beta*S*I
    
    return(list(c(dS,dE,dI,dR,dInc)))
    #return(list(c(dS,dI,dR,dInc)))
  }
  
  init.Inc <- init.state[["Inc"]]
  BetaTT <- init.theta[[2]][["beta"]]
  traj <- data.frame('time'=1,'S'=init.state[["S"]],'E'=init.state[["E"]],'I'=init.state[["I"]],'R'=init.state[["R"]],'Inc'=init.state[["Inc"]])
  #traj <- data.frame('time'=0,'S'=init.state[["S"]],'I'=init.state[["I"]],'R'=init.state[["R"]],'Inc'=init.state[["Inc"]])
  Beta <- data.frame('beta'= matrix(1,nrow = length(times), ncol = 1))
  for(t in times){
    betaT <- BetaTT[t-min(times)+1]
    Beta[t-min(times)+1,1] <- betaT
    init.theta[[1]][["beta"]] <- betaT
    if((t+1)<=max(times)){
      
      index = t-min(times)+1
      init.state[["S"]] = traj[index,2]
      init.state[["E"]] = traj[index,3]
      init.state[["I"]] = traj[index,4]
      init.state[["R"]] = traj[index,5]
      init.state[["Inc"]] = traj[index,6]
      
      trajUnit <- as.data.frame(ode(init.state, c(t-min(times)+1,t+1-min(times)+1), SEIR_ode, init.theta, method = "ode45"))
      trajUnit <- trajUnit[-1,]
      traj <- rbind(traj, trajUnit)
    }
    
  }
  
  traj <- cbind(traj, Beta)
  
  #traj <- as.data.frame(ode(init.state, times, SEIR_ode, init.theta, method = "ode45"))
  
  # compute incidence of each time interval
  traj$Inc <- c(init.Inc, diff(traj$Inc))
  
  return(traj)
}

PTR1$rPointObs <- function(model.point, theta){
  
  obs.point <- rpois(n=1, lambda=theta[[1]][["RR"]]*model.point[["Inc"]])
  #obs.point <- rnorm(n=1, mean = theta[[1]][["RR"]]*model.point[["Inc"]], sd=0.1)
  
  return(c(obs=obs.point))
}

################################################################################
################################################################################
setwd("E:/科研/paper/bi-modal/mfiidd")
data <- read.csv(header = TRUE,"USA_2009.csv")
#plot(x = data$Week.Number.2009, y = data$Confirmed.Cases.2009, type = "l", lty=1 ,col = "red", lwd=2, xlab = "Week", ylab = "Number of cases")
init.theta <- c(D_inf = 3, D_lat = 1)
S0 <- 0.8 #按理应该与init.state1的0.857保持一致，但0.945的效果最佳
I0 <- 0.01

#calculate β(t)
nu <- 1/init.theta[["D_inf"]]
alpha <- 1/init.theta[["D_lat"]]
data$Confirmed.Cases.2009 <- data$Confirmed.Cases.2009*700/287000000
#plot(x = data$Week.Number.2009, y = data$Confirmed.Cases.2009, type = "l", lty=1 ,col = "red", lwd=2, xlab = "Week", ylab = "Number of cases")
y = data$Confirmed.Cases.2009 # the number of new cases per day
#x = data$Week.Number.2009
x = c(1:length(data$Week.Number.2009))
#plot(x, y, type = "l" ,col = "red", lwd=2, xlab = "Time", ylab = "Fraction of population")
#lines(spline(x, y))
#lines(spline(x, y, n = 201), col = 2)

Meth <- "monoH.FC" #method of splinefun
#f is the function of the time series data on the number of infected individuals
f <- splinefun(x = x, y = y, method = Meth) 
Inc <- f(x) #values of function f
#lines(x,Inc)
#deriv=1对应一阶导函数，2对应二阶导函数，3对应三阶导函数。derive只能在0,1,2,3中取值
#Inc1 <- f(x,deriv = 1) #values of the first derivative of function f
#f1 is the first derivative function of function f
#f1 <- splinefun(x = x, y = Inc1, method = Meth)
G <- Inc/alpha # values of function g
g <- splinefun(x = x, y = G, method = Meth) # g is the function of function g
G1 <- g(x,deriv = 1) #values of the first derivative of function g
g1 <- splinefun(x=x, y=G1, method = Meth) # g1 is the first derivative function of function g
H <- G1 + alpha*G # values of function h
h <- splinefun(x=x, y=H, method = Meth) # h is the function of function h
#H[which(H<=0)]<- 0 #H中的负数值都与0接近，因此将它们都置零，目的是使Inc非负
eInc1 <- exp(nu*x)*Inc # values of function EInC1
#eInc1[which(eInc1<0)]
EInC1 <- splinefun(x=x, y=eInc1, method = Meth) # function of function EInC1
SInH <- H
It <- H
BetaT <- data.frame('beta'= matrix(1,nrow = length(x), ncol = 1))

for(t in x){
  T <- t-min(x)+1
  SInH[T] <- S0 - integrate(h,0,T)$value
  #SInH <- abs(SInH)
  It[T] <- exp(-nu*T)*(I0 + integrate(EInC1,0,T)$value)
  #if(G[t]<=0) { G[t] <- -G[t] }
  BetaT[T,1] <- H[T]/(It[T]*SInH[T])
  
}
#plot(x = x, y = BetaT$beta,type = "l", lty=1 ,col = "red", lwd=2, xlab = "Time", ylab = "Transmission rate")
which(G<=0)
which(H<0)
which(SInH<0)
which(It<0)
#which(BetaT<0)
min(SInH)

#############################################################################################
#####use the generated β(t) to simulate the time series of the number of infected individuals

POP <- 27700
S0 <- 0.8
E0 <- 0.002
#init.state1 <- c('S' = 0.857, 'E'=0.002, 'I' = 0.004, 'R' = 0.137,'Inc' = data[1,2])
#init.state1 <- c('S' = S0*POP, 'E'=0.002*POP, 'I' = I0*POP, 'R' = 0.137*POP,'Inc' = data[1,2]*POP)
init.state1 <- c('S' = S0*POP, 'E'=E0*POP, 'I' = I0*POP, 'R' = (1-S0-I0-E0)*POP,'Inc' = data[1,2]*POP)
init.theta1 <- list(c(D_inf = init.theta[["D_inf"]], D_lat = init.theta[["D_lat"]], RR = 0.75),BetaT)
x=data$Week.Number.2009
epi_time = x
obs.traj1 <- rTrajObs(fitmodel = PTR1, theta = init.theta1, init.state = init.state1, times = epi_time)
#traj <- PTR1$simulate(init.theta=init.theta1, init.state=init.state1, times=epi_time)
plot(x = epi_time, y = obs.traj1$Inc, type = "l", lty=1 ,col = "red", lwd=2, xlab = "Week", ylab = "Number of cases")
lines(x = epi_time, y = obs.traj1$obs, type = "l",lty=1, col = "black",lwd=1)
#legend(20,600, legend = c("Incidence","Incidence observation","Transmission rate"),bty = "n", lty = c(1,1,2), lwd = c(2,1,2),col = c("red","black","dark gray") )
legend("topleft", legend = c("Incidence","Incidence observation","Transmission rate"),bty = "n", lty = c(1,1,2), lwd = c(2,1,2),col = c("red","black","dark gray") )
par(new=TRUE)
plot(x=epi_time, y=obs.traj1$beta, type = "l", lty=2, col="dark gray", lwd=2, yaxt="n", ylab = "",ylim = c(min(obs.traj1$beta),max(obs.traj1$beta)),xaxt="n", xlab = "" )
axis(side = 4)

# add a title for the right axis 
mtext("Transmission rate", side=4, line=2.5, cex.lab=1, las=0, col="black")

AT = sum(obs.traj1$Inc)/POP*S0

##################################################################################################
##从该行起至204行，目的是使出的图的横坐标轴范围显示在1~100之间。可暂时不用，直接跳过。
x = seq(1, 99, by=2.8) #time range[0,100]
x1 = seq(1,99,by=0.01)
x2 = seq(1,99,by=2)
x3 = seq(1,99,by=5.6)
IncF <- splinefun(x=x, y=obs.traj1$Inc, method = "monoH.FC")
#obsF <- splinefun(x=x, y=obs.traj1$obs, method = "monoH.FC")
obsF <- splinefun(x=x, y=y*POP, method = "monoH.FC")
betaF <- splinefun(x=x, y=obs.traj1$beta, method = "monoH.FC")
IncFV <- IncF(x1, deriv = 0)
obsFV <- obsF(x2, deriv = 0)
betaFV <- betaF(x3,deriv = 0)

par(mar=c(5, 4, 2, 4))
#ylim的下限取0是因为考虑到obs.traj1$Inc的显示
#plot(x = x, y = y, type = "p", lty=1 ,col = "red", lwd=2, xlab = "Time", ylab = "Fraction of population")
#plot(x = x, y = y*POP, type = "p", lty=1 ,col = "red", lwd=2, xlab = "Time", ylab = "Fraction of population")
#lines(x = x, y = obs.traj1$I, type = "l", lty=1 ,col = "blue", lwd=2, xlab = "Time", ylab = "Fraction of population",ylim = c(0,max(obs.traj1$I)))
##plot(x = x1, y = IncFV, type = "l", lty=1 ,col = "red", lwd=2, xlab = "Day", ylab = "Number of cases", ylim = c(min(y*POP), max(y*POP)), xlim = c(0,max(x)+4))
plot(x = x1, y = IncFV, type = "l", lty=1 ,col = "red", lwd=2, xlab = "Day", ylab = "Number of cases", xlim = c(0,max(x)+4))

lines(x = x2, y = obsFV, type = "l",lty=1, col = "black",lwd=1)
#lines(x = x, y = y*POP, type = "p",lty=1, col = "black",lwd=2)
#legend("topright", legend = c("infected(simulated)","infected(obs)","incidence(simulated)"),lty = c(1,1,1), lwd = c(2,2,2),col = c("blue","red","green") )

legend(20,600, legend = c("Incidence","Incidence observation","Transmission rate"),bty = "n", lty = c(1,1,2), lwd = c(2,1,2),col = c("red","black","dark gray") )

par(new=TRUE)
plot(x=x3, y=betaFV, type = "l", lty=2, col="dark gray", lwd=2, yaxt="n", ylab = "",ylim = c(min(obs.traj1$beta),max(obs.traj1$beta)),xaxt="n", xlab = "" )
axis(side = 4)

# add a title for the right axis 
mtext("Transmission rate", side=4, line=2.5, cex.lab=1, las=0, col="black")
##################################################################################################

####################################################################################################################
##判断two wave
##########################方案2
GT = init.theta[["D_lat"]] + 0.5*init.theta[["D_inf"]]#generation time
final_size = sum(obs.traj1$Inc)
i = which.max(obs.traj1$Inc) #最大值位于从左往右的第几个（不一定是time）
Si = obs.traj1$Inc[i]
#towards right
if( (i+2) <= length(obs.traj1$Inc) ) {
  len = length(c((i+2):length(obs.traj1$Inc)))
  TWr = data.frame(matrix(999,nrow = len, ncol = 11))
  TWr[,5] = i
  colnames(TWr)=c('TW','Si','Sj','vij','i','j','|j-i|','P2','TPc','GT','W2')#'TW' represents a new metric of mine inspired from the reference paper, '|j-i|' represents time gap between two peaks. 'P2' represents 
  #(Sj-vij)/Sj, Sj is the value of the second peak. 'TPc' represents another metric suggested by the author(Thomas J Hladish). 'GT' represents generation time.'W2' represents (Sj-vij)/(Si-vij)
  # The title of the paper from which I reference is "Epidemic Wave Dynamics Attributable to Urban Community Structure:A Theoretical Characterization of Disease Transmission in a Large Network"
  t = 1
  for(j in (i+2):length(obs.traj1$Inc)){
    Sij = obs.traj1$Inc[c(i:j)]
    vij = min(Sij)
    Sj = obs.traj1$Inc[j]
    TWr[t,2] = Si
    TWr[t,3] = Sj
    TWr[t,4] = vij
    TWr[t,6] = j
    TWr[t,7] = j-i
    if((Si>0)&&(Sj>0)){ # &&是标量的逻辑与运算，&是向量的逻辑与运算
      TWr[t,1] = (((Si-vij)/Si)*((Sj-vij)/Sj))^(1/2)
    }
    else{
      TWr[t,1] = 0
    }
    if(Sj>0){
      TWr[t,8] = (Sj-vij)/Sj
    }
    else{
      TWr[t,8] = 0
    }
    TWr[t,9] = ((Si-vij)*(Sj-vij)/final_size)^(1/2)
    TWr[t,10] = GT
    if((Si-vij)>0){
      TWr[t,11] = (Sj-vij)/(Si-vij)
    }
    t=t+1
  }
}

#towards left
if(i>2){
  len = i-2
  TWl = data.frame(matrix(999,nrow = len, ncol = 11))
  TWl[,5] = i
  colnames(TWl)=c('TW','Si','Sj','vij','i','j','|j-i|','P2','TPc','GT','W2')
  t = 1
  for(j in 1:(i-2)){
    Sij = obs.traj1$Inc[c(j:i)]
    vij = min(Sij)
    Sj = obs.traj1$Inc[j]
    TWl[t,2] = Si
    TWl[t,3] = Sj
    TWl[t,4] = vij
    TWl[t,6] = j
    TWl[t,7] = i-j
    if((Si>0)&&(Sj>0)){
      TWl[t,1] = (((Si-vij)/Si)*((Sj-vij)/Sj))^(1/2)
    }
    else{
      TWl[t,1] = 0
    }
    if(Sj>0){
      TWl[t,8] = (Sj-vij)/Sj
    }
    else{
      TWl[t,8] = 0
    }
    TWl[t,9] = ((Si-vij)*(Sj-vij)/final_size)^(1/2)
    TWl[t,10] = GT
    if((Si-vij)>0){
      TWl[t,11] = (Sj-vij)/(Si-vij)
    }
    t=t+1
  }
}

#combine
TW = rbind(TWl,TWr)
ind = which.max(TW[,8])
ind
P2 = TW$P2[ind] #(Sj-vij)/Sj, Sj is the value of the second peak
W2 = TW$W2[ind] #(Sj-vij)/(Si-vij)

thre1 = 0.5
thre2 = 10
thre3 = 1
#TW$TW[ind] #2-peak metric
#TW$`|j-i|`[ind] #time gap between two peaks
if(P2 >= thre1) {
  if((TW$`|j-i|`[ind] >= thre3*GT) && (W2 >= (1/thre2))){peakNum = 2}
  else {
    if((TW$`|j-i|`[ind] < thre3*GT) && (W2 >= (1/thre2))) {peakNum = 1.8} #time gap between two peaks is not long enough to generate a second generation
    else {peakNum = 1} #只要W2 < (1/thre2)，那么就认为是单峰
    
  }
}

if((P2 < thre1)&(P2>0)){
  if((TW$`|j-i|`[ind] >= thre3*GT) && (W2 >= (1/thre2))) {peakNum = 1.5} #the valley is not deep enough to be looked as an interwave period, but the time gap is enough
  else {
    if((TW$`|j-i|`[ind] < thre3*GT) && (W2 >= (1/thre2))) {peakNum = 1.2}
    else {peakNum = 1} #只要W2 < (1/thre2)，那么就认为是单峰
  }
}
if(P2==0) {peakNum = 1}#只要P2==0，那么也认为是单峰
peakNum
