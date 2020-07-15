library(plyr)
library(fitR)
###本程序可同时适用发病人数和发病人数比例两种情形，互换时需对line 10的SEIR_ode、line 135、137、145进行修改
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
    
    #dS <- -beta*S*I/N 
    dS <- -beta*S*I
    #dE <- beta*S*I/N - epsilon*E
    dE <- beta*S*I - epsilon*E
    dI <- epsilon*E - nu*I
    #dI <- beta*S*I - nu*I
    dR <- nu*I
    dInc <- epsilon*E
    #dInc <- beta*S*I
    
    return(list(c(dS,dE,dI,dR,dInc)))
    #return(list(c(dS,dI,dR,dInc)))
  }
  
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
  traj$Inc <- c(0, diff(traj$Inc))
  
  return(traj)
}

PTR1$rPointObs <- function(model.point, theta){
  
  #obs.point <- rpois(n=1, lambda=theta[["RR"]]*model.point[["Inc"]])
  obs.point <- rnorm(n=1, mean = theta[[1]][["RR"]]*model.point[["Inc"]], sd=0.1)
  
  return(c(obs=obs.point))
}

################################################################################
################################################################################
data <- read.csv(header = TRUE,"USA_2009.csv")
init.theta <- c(beta0 = 1.2, D_inf = 3, D_lat = 1, RR = 0.75)
#calculate β(t)
nu <- 1/init.theta[["D_inf"]]
alpha <- 1/init.theta[["D_lat"]]
data$Confirmed.Cases.2009 <- data$Confirmed.Cases.2009*700/287000000
y = data$Confirmed.Cases.2009 # the number of new cases per day
#x = data$Week.Number.2009
x = c(1:length(data$Week.Number.2009))
#plot(x, y, type = "l" ,col = "red", lwd=2, xlab = "Time", ylab = "Fraction of population")
#lines(spline(x, y))
#lines(spline(x, y, n = 201), col = 2)

Meth <- "monoH.FC" #method of splinefun
#f is the function of the time series data on the number of infected individuals
f <- splinefun(x = x, y = y, method = Meth) 
I <- f(x) #values of function f
#lines(x,I)
#deriv=1对应一阶导函数，2对应二阶导函数，3对应三阶导函数。derive只能在0,1,2,3中取值
I1 <- f(x,deriv = 1) #values of the first derivative of function f
#f1 is the first derivative function of function f
f1 <- splinefun(x = x, y = I1, method = Meth)
G <- (I1 + nu * I)/alpha # values of function g
g <- splinefun(x = x, y = G, method = Meth) # g is the function of function g
G1 <- g(x,deriv = 1) #values of the first derivative of function g
g1 <- splinefun(x=x, y=G1, method = Meth) # g1 is the first derivative function of function g
H <- G1 + alpha*G # values of function h
h <- splinefun(x=x, y=H, method = Meth) # h is the function of function h
H[which(H<=0)]<- 0 #H中的负数值都与0接近，因此将它们都置零，目的是使Inc非负
SInH <- H
BetaT <- data.frame('beta'= matrix(1,nrow = length(x), ncol = 1))
S0 <- 0.857 #按理应该与init.state1的0.857保持一致，但0.945的效果最佳
for(t in x){
  T <- t-min(x)+1
  SInH[T] <- S0 - integrate(h,0,T)$value
  #SInH <- abs(SInH)
  #if(G[t]<=0) { G[t] <- -G[t] }
  BetaT[T,1] <- H[T]/(I[T]*SInH[T])
  
}
#plot(x = x, y = BetaT$beta,type = "l", lty=1 ,col = "red", lwd=2, xlab = "Time", ylab = "Transmission rate")
which(G<=0)
which(H<0)
which(SInH<=0)
#which(BetaT<=0)
min(SInH)

#############################################################################################
#####use the generated β(t) to simulate the time series of the number of infected individuals

POP <- 27700
init.state1 <- c('S' = 0.857, 'E'=0.002, 'I' = 0.004, 'R' = 0.137,'Inc' = 0)
#init.state1 <- c('S' = 0.857*POP, 'E'=0.002*POP, 'I' = 0.004*POP, 'R' = 0.137*POP,'Inc' = 0*POP)
init.theta1 <- list(c(beta0 = 1.2, D_inf = 3, D_lat = 1, RR = 0.75),BetaT)
x=data$Week.Number.2009
epi_time = x
obs.traj1 <- rTrajObs(fitmodel = PTR1, theta = init.theta1, init.state = init.state1, times = epi_time)
#traj <- PTR1$simulate(init.theta=init.theta1, init.state=init.state1, times=epi_time)

#ylim的下限取0是因为考虑到obs.traj1$Inc的显示
plot(x = x, y = y, type = "p", lty=1 ,col = "red", lwd=2, xlab = "Time", ylab = "Fraction of population")
lines(x = x, y = obs.traj1$I, type = "l", lty=1 ,col = "blue", lwd=2, xlab = "Time", ylab = "Fraction of population",ylim = c(0,max(obs.traj1$I)))
lines(x = x, y = obs.traj1$Inc, type = "l", lty=1 ,col = "green", lwd=2, xlab = "Time", ylab = "Fraction of population")
#legend("topright", legend = c("infected(simulated)","infected(obs)","incidence(simulated)"),lty = c(1,1,1), lwd = c(2,2,2),col = c("blue","red","green") )
