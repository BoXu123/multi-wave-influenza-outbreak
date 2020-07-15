rm(list = ls())
library(fitR)
library(plyr)
##########age.Linear
##################model constrution
data(SEITL_deter)
age.Linear <- SEITL_deter
age.Linear$simulate <- function(theta,init.state,times) {
  
  age.Linear_ode <- function(time, state, parameters) {
    
    ## parameters
    beta11 <- parameters[["beta11"]]
    beta12 <- parameters[["beta12"]]
    beta21 <- parameters[["beta21"]]
    beta22 <- parameters[["beta22"]]
    l <- parameters[["l"]]
    epsilon <- 1/parameters[["D_lat"]]
    nu <- 1/parameters[["D_inf"]]
    
    ## states
    S1 <- state[["S1"]]
    E1 <- state[["E1"]]
    I1 <- state[["I1"]]
    R1 <- state[["R1"]]
    Inc1 <- state[["Inc1"]]
    S2 <- state[["S2"]]
    E2 <- state[["E2"]]
    I2 <- state[["I2"]]
    R2 <- state[["R2"]]
    Inc2 <- state[["Inc2"]]
    Inc <- state[["Inc"]]
    
    N <- S1 + E1 + I1 + R1 + S2 + E2 + I2 + R2
    
    dS1 <- -S1*(beta11*I1 + beta12*I2)/N - l*S1
    dE1 <- S1*(beta11*I1 + beta12*I2)/N - epsilon*E1 - l*E1
    dI1 <- epsilon*E1 - nu * I1 - l*I1
    dR1 <- nu * I1 - l*R1
    dS2 <- l*S1 - S2*(beta21*I1 + beta22*I2)/N
    dE2 <- l*E1 + S2*(beta21*I1 + beta22*I2)/N - epsilon*E2
    dI2 <- l*I1 + epsilon*E2 - nu*I2
    dR2 <- nu*I2 + l*R1
    dInc1 <- epsilon*E1
    dInc2 <- l*I1 + epsilon*E2
    dInc <- epsilon*E1 + (l*I1 + epsilon*E2)
    
    
    return(list(c(dS1, dE1, dI1, dR1, dS2, dE2, dI2, dR2, dInc1, dInc2, dInc)))
  }
  
  T <- theta[["T"]]
  a1 <- theta[["a1"]]
  b1 <- theta[["b1"]]
  a2 <- theta[["a2"]]
  b2 <- theta[["b2"]]
  Tint <- theta[["Tint"]]
  ####wave 1 (before time T)
  times1 <- c(1:T)
  N=init.state[["N"]]
  init.state1 <- c('S1'=init.state[["S1a"]]*N,'E1'=init.state[["E1a"]]*N,'I1'=init.state[["I1a"]]*N,'R1'=init.state[["R1a"]]*N,
                   'S2'=init.state[["S2a"]]*N,'E2'=init.state[["E2a"]]*N,'I2'=init.state[["I2a"]]*N,'R2'=init.state[["R2a"]]*N,
                   'Inc1'=init.state[["Inc1"]]*N,'Inc2'=init.state[["Inc2"]]*N,'Inc'=init.state[["Inc"]]*N )
  traj1 <- data.frame('time'=min(times1),'S1'=init.state[["S1a"]]*N,'E1'=init.state[["E1a"]]*N,'I1'=init.state[["I1a"]]*N,'R1'=init.state[["R1a"]]*N,
                     'S2'=init.state[["S2a"]]*N,'E2'=init.state[["E2a"]]*N,'I2'=init.state[["I2a"]]*N,'R2'=init.state[["R2a"]]*N,
                     'Inc1'=init.state[["Inc1"]]*N,'Inc2'=init.state[["Inc2"]]*N,'Inc'=init.state[["Inc"]]*N )
  N = sum(traj1[1, 2:9])
  s.chil <- data.frame('s.chil'= matrix(1,nrow = length(times), ncol = 1))
  for(t in times1){
    s.prop <- a1*t + b1 #线性下降
    s.chil[t-min(times1)+1,1] <- s.prop
    if((t+1)<=max(times1)){
      
      index = t-min(times1)+1
      init.state1[["S1"]] = s.prop*N #此处可以考虑采用"s.prop*traj1[index,2]"
      init.state1[["E1"]] = traj1[index,3]
      init.state1[["I1"]] = traj1[index,4]
      init.state1[["R1"]] = N - sum(traj1[index, 6:9]) - sum(traj1[index, 3:4]) - init.state1[["S1"]] #因为S1变少的部分都变成了R1
      #init.state1[["R1"]] = traj1[index,5]
      init.state1[["S2"]] = traj1[index,6]
      init.state1[["E2"]] = traj1[index,7]
      init.state1[["I2"]] = traj1[index,8]
      init.state1[["R2"]] = traj1[index,9]
      init.state1[["Inc1"]] = traj1[index,10]
      init.state1[["Inc2"]] = traj1[index,11]
      init.state1[["Inc"]] = traj1[index,12]
      
      trajUnit <- as.data.frame(ode(init.state1, c(t,t+1), age.Linear_ode, theta, method = "ode45"))
      trajUnit <- trajUnit[-1,]
      traj1 <- rbind(traj1, trajUnit)
      #N = sum(traj1[index, 2:9])
    }
  }
##############  
  ####wave 2 (between time T and time T+Tint)
  times2 <- c((T+1):(T+1+Tint))
  init.state2 <- c('S1'=s.chil[max(times1),1]*N,'E1'=traj1$E1[T],'I1'=traj1$I1[T],'R1'=traj1$R1[T],
                   'S2'=traj1$S2[T],'E2'=traj1$E2[T],'I2'=traj1$I2[T],'R2'=traj1$R2[T],
                   'Inc1'=traj1$Inc1[T],'Inc2'=traj1$Inc2[T],'Inc'=traj1$Inc[T] )
  traj2 <- data.frame('time'=min(times2),'S1'=s.chil[max(times1),1]*N,'E1'=traj1$E1[T],'I1'=traj1$I1[T],'R1'=traj1$R1[T],
                      'S2'=traj1$S2[T],'E2'=traj1$E2[T],'I2'=traj1$I2[T],'R2'=traj1$R2[T],
                      'Inc1'=traj1$Inc1[T],'Inc2'=traj1$Inc2[T],'Inc'=traj1$Inc[T] )
  for(t in times2){
    s.prop <- a2*t + b2 #线性下降
    s.chil[t-min(times1)+1,1] <- s.prop
    if((t+1)<=max(times2)){
      
      index = t-min(times2)+1
      init.state2[["S1"]] = s.prop*N #此处可以考虑采用"s.prop*traj2[index,2]"
      init.state2[["E1"]] = traj2[index,3]
      init.state2[["I1"]] = traj2[index,4]
      init.state2[["R1"]] = N - sum(traj2[index, 6:9]) - sum(traj2[index, 3:4]) - init.state2[["S1"]]
      #init.state2[["R1"]] = traj2[index,5]
      init.state2[["S2"]] = traj2[index,6]
      init.state2[["E2"]] = traj2[index,7]
      init.state2[["I2"]] = traj2[index,8]
      init.state2[["R2"]] = traj2[index,9]
      init.state2[["Inc1"]] = traj2[index,10]
      init.state2[["Inc2"]] = traj2[index,11]
      init.state2[["Inc"]] = traj2[index,12]
      
      trajUnit <- as.data.frame(ode(init.state2, c(t,t+1), age.Linear_ode, theta, method = "ode45"))
      trajUnit <- trajUnit[-1,]
      traj2 <- rbind(traj2, trajUnit)
      #N = sum(traj2[index, 2:9])
    }
  }
  #traj2 <- traj2[-1,]
  
  ####wave 3 (after time T+Tint)
  times3 <- c((T+1+Tint+1):max(times))
  for(t in times3) {
    s.prop <- s.chil[max(times2),1] #不变
    s.chil[t-min(times1)+1,1] <- s.prop
  }
  
  #init.state[["S1b"]]*N
  init.state3 <- c('S1'=s.chil[max(times2),1]*N,'E1'=traj2$E1[1+Tint],'I1'=traj2$I1[1+Tint],'R1'=traj2$R1[1+Tint],
                  'S2'=traj2$S2[1+Tint],'E2'=traj2$E2[1+Tint],'I2'=traj2$I2[1+Tint],'R2'=traj2$R2[1+Tint],
                  'Inc1'=traj2$Inc1[1+Tint],'Inc2'=traj2$Inc2[1+Tint],'Inc'=traj2$Inc[1+Tint] )
  
  traj3 <- as.data.frame(ode(init.state3, times3, age.Linear_ode, theta, method = "ode45"))
  traj3 <- traj3[-1,]
  s.chil <- s.chil[-c(min(times3)),]
#############
  #combine the two data frames
  traj <- rbind(traj1, traj2, traj3)
  
  traj <- cbind(traj, s.chil)
  traj <- traj[-min(times2),]
  # compute incidence of each time interval
  traj$Inc1 <- c(0, diff(traj$Inc1))
  traj$Inc2 <- c(0, diff(traj$Inc2))
  traj$Inc <- c(0, diff(traj$Inc))
  
  ####wave 3 (after time T+Tint)
  #times3 <- c((T+1+Tint+1):max(times))
  #init.state3 <- c('S1'=init.state[["S1c"]]*N,'E1'=traj2$E1[1+Tint],'I1'=traj2$I1[1+Tint],'R1'=traj2$R1[1+Tint],
  #                 'S2'=traj2$S2[1+Tint],'E2'=traj2$E2[1+Tint],'I2'=traj2$I2[1+Tint],'R2'=traj2$R2[1+Tint],
  #                 'Inc1'=traj2$Inc1[1+Tint],'Inc2'=traj2$Inc2[1+Tint],'Inc'=traj2$Inc[1+Tint] )
  #traj3 <- as.data.frame(ode(init.state3, times3, age.Linear_ode, theta, method = "ode45"))
  # compute incidence of each time interval
  #traj3$Inc1 <- c(init.state3["Inc1"], diff(traj3$Inc1))
  #traj3$Inc2 <- c(init.state3["Inc2"], diff(traj3$Inc2))
  #traj3$Inc <- c(init.state3["Inc"], diff(traj3$Inc))
  
  
  #combine the two data frames
  #traj <- rbind(traj1, traj2, traj3)
  
  return(traj)
}

age.Linear$rPointObs <- function(model.point, theta){
  
  ## the prevalence is observed through a Poisson process
  #obs.point <- rpois(n = 1, lambda = model.point[["Inc"]] * theta[["RR"]])
  obs.point <- rnorm(n=1, mean = model.point[["Inc"]] * theta[["RR"]], sd=0.2)
  
  return(c(obs = obs.point))
}

####
epi_time = c(1:200)
#epi_time = c(1:365)
#init.state <- c('N'=27700,'S1a' = 0.3, 'I1a' = 0.001, 'R1a' = 0,'S2a' = 0.35, 'I2a' = 0.001, 'R2a' = 0,'Inc1' = 0,'Inc2' = 0,'Inc' = 0, 'S1b'=0.24,'S1c'=0.075)
init.state <- c('N'=27700,'S1a' = 0.3, 'E1a'=0, 'I1a' = 0.0012, 'R1a' = 0,'S2a' = 0.22, 'E2a'=0, 'I2a' = 0.001, 'R2a' = 0.02,'Inc1' = 0,'Inc2' = 0,'Inc' = 0)
#init.theta <- c(beta11 = 100/365, beta12 = 10/365, beta21 = 10/365, beta22 = 40/365, l = 1/(15*365), nu = 10/365,RR=0.75, T=27, Tint=14)
init.theta <- c(beta11 = 170/365, beta12 = 10/365, beta21 = 10/365, beta22 = 53/365, l = 0, D_lat = 365/150, D_inf = 365/10, RR=0.75, T=27, Tint=18, 
                a1=-0.0023, b1=0.3023, a2=-0.01439, b2=0.665)
#obs.traj <- age.Linear$simulate(init.theta, init.state, times = epi_time)
obs.traj <- rTrajObs(fitmodel = age.Linear, theta = init.theta, init.state = init.state, times = epi_time)
#obs.traj <- obs.traj[-c(init.theta[["T"]]+init.theta[["Tint"]]+1),]
###########################################################
####plot designed by Bo Xu
#epi_time <- seq(0.5,100,by=0.5)
par(mar=c(5, 4, 2, 4))
#lines(x = epi_time[-c(length(epi_time), length(epi_time)-1)], y = obs.traj$Inc, type = "l", lty=1, col="black" )
plot(x = epi_time[-c(length(epi_time), length(epi_time)-1)], y = obs.traj$Inc, type = "l", lty=1 ,col = "red", lwd=2, xlab = "Day", ylab = "Number of cases")#Inc is incidence of children and adults calculated by model
#plot(x = epi_time[-c(length(epi_time))], y = obs.traj$s.chil, type = "l", lty=1 ,col = "red", lwd=2, xlab = "Day", ylab = "s.chil")
lines(x=epi_time[-c(length(epi_time), length(epi_time)-1)], y=obs.traj$Inc1, type = "l",lty=2, col="dark gray", lwd=2 )#Inc1 is incidence of children calculated by model
lines(x=epi_time[-c(length(epi_time), length(epi_time)-1)], y=obs.traj$Inc2, type = "l", lty=3, col="dark gray", lwd=2 )#Inc2 is incidence of adults calculated by model
#abline(v = init.theta[["T"]], type = "l", lty=2, col="blue",lwd=0.5)
#abline(v = init.theta[["T"]]+init.theta[["Tint"]], type = "l", lty=2, col="blue",lwd=0.5)
lines(x = epi_time[-c(length(epi_time), length(epi_time)-1)], y = obs.traj$obs, type = "l", lty=1,col = "black",lwd=2)#obs is from Inc considering reporting rate and poisson observation process
#legend("topright", legend = c("Incidence","Incidence observation","Changing time","Reporting rate"),bty = "n", lty = c(1,1,2,2), lwd = c(2,2,0.5,1),col = c("red","black","blue","gray") )
legend("topright", legend = c("Incidence","Incidence of children","Incidence of adults","Incidence observation","susceptibility of children"),bty = "n", lty = c(1,2,3,1,2), lwd = c(2,2,2,2,2),col = c("red","dark gray","dark gray","black","blue") )

par(new=TRUE)
plot(x=epi_time[-c(length(epi_time), length(epi_time)-1)], y=obs.traj$s.chil, type = "l", lty=2, col="blue", lwd=2, yaxt="n", ylab = "",ylim = c(min(obs.traj$s.chil),max(obs.traj$s.chil)+0.2),xaxt="n", xlab = "" )
axis(side = 4)

# add a title for the right axis 
mtext("susceptibility of children", side=4, line=2.5, cex.lab=1, las=0, col="black")

##################################################################################################################
##判定双峰
##########################方案3
GT = init.theta[["D_lat"]] + 0.5*init.theta[["D_inf"]]#generation time
#GT = init.theta[["D_lat"]] + init.theta[["D_inf"]] + init.theta[["D_imm"]]
#GT = init.theta[["D_lat"]] + init.theta[["D_inf"]] + init.theta[["D_imm"]] + init.theta[["D_win"]]
#GT = (1/init.theta[["k1"]]) + (1/init.theta[["gamma11"]]) #针对TCC_RRC_deter.R
final_size = sum(obs.traj$Inc)
i = which.max(obs.traj$Inc) #最大值位于从左往右的第几个（不一定是time）
Si = obs.traj$Inc[i]
#towards right
if( (i+2) <= length(obs.traj$Inc) ) {
  len = length(c((i+2):length(obs.traj$Inc)))
  TWr = data.frame(matrix(999,nrow = len, ncol = 11))
  TWr[,5] = i
  colnames(TWr)=c('TW','Si','Sj','vij','i','j','|j-i|','P2','TPc','GT','W2')#'TW' represents a new metric of mine inspired from the reference paper, '|j-i|' represents time gap between two peaks. 'P2' represents 
  #(Sj-vij)/Sj, Sj is the value of the second peak. 'TPc' represents another metric suggested by the author(Thomas J Hladish). 'GT' represents generation time.'W2' represents (Sj-vij)/(Si-vij)
  # The title of the paper from which I reference is "Epidemic Wave Dynamics Attributable to Urban Community Structure:A Theoretical Characterization of Disease Transmission in a Large Network"
  t = 1
  for(j in (i+2):length(obs.traj$Inc)){
    Sij = obs.traj$Inc[c(i:j)]
    vij = min(Sij)
    Sj = obs.traj$Inc[j]
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
    Sij = obs.traj$Inc[c(j:i)]
    vij = min(Sij)
    Sj = obs.traj$Inc[j]
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
  if( (i+2) <= length(obs.traj$Inc) ){#若条件满足，则表明TWr存在
    #combine
    TW = rbind(TWl,TWr)
  }else{
    TW = TWl
  }
}else{
  TW = TWr
}

#ind1 = which.max(TW[,8]) #若TW[,8]的最大值不止一个（即多值重复），则which.max(TW[,8])只会返回一个位置值，而which(TW[,8] == max(TW[,8]))会全部返回所有位置值
###以下9行的修改是在读取Tristan da cunha的实际发病数据后做出的
ind1 = which(TW[,8] == max(TW[,8])) #TW[,8]是P2; ind1所记录的都是TW[,8]或P2的最大值所在的位置
if(length(ind1) == 1){ #如果TW[,8]或P2的最大值只出现在一个位置上且该最大值不是1,则直接选取它
  if(TW[ind1,8] != 1) {ind = ind1}
  if(TW[ind1,8] == 1) {
    TW[ind1,11]
    #将其设为零，更新之后从头再来
    TW[ind1,8] = 0
    ind66 = which(TW[,8] == max(TW[,8])) #TW[,8]是P2; ind66所记录的都是更新之后TW[,8]或P2的最大值所在的位置
    if(length(ind66) > 1){#如果TW[,8]或P2有不止一个相等的最大值
      ind77 = which(TW[ind66,11] == max(TW[ind66,11])) #ind66所记录的都是更新之后TW[,8]或P2的最大值所在的位置,在ind66中找TW[,11]或W2的最大值在ind66中的位置
      if(length(ind77) == 1){ #若在ind66中，TW[,11]或W2的最大值只在一个位置出现，则选取它
        ind00 = ind66[ind77]
      }
      if(length(ind77) > 1){ #若在ind66中，TW[,11]或W2的最大值在多个位置出现，则在这些位置中找到TW[,7]或|i-j|最大的那个位置
        ind88 = which.max(TW[ind66[ind77],7])
        ind00 = ind66[ind77][ind88]
      }
    }
    if(length(ind66) == 1){
      ind00 = ind66
    }
    #比较ind1和ind00各自对应的TW[,11]或W2的大小
    if(TW[ind1,11] >= TW[ind00,11]){
      ind = ind1
      TW[ind1,8] = 1 #将之前设为0的TW[ind,8]设置回原值1
    }
    if(TW[ind1,11] < TW[ind00,11]){
      ind = ind00
    }
    #在TW[,8]或P2的多个相等的最大值等于1的前提下,将之前设为0的TW[ind1,8]设置回原值1
    TW[ind1,8] = 1
  }
}

if(length(ind1) > 1){ #如果TW[,8]或P2的最大值出现在不止一个位置上
  #####################################################
  #若TW[,8]或P2的多个相等的最大值等于1
  if(TW[ind1,8][1] == 1){ #P2等于1表明vij等于0
    ind2 = which(TW[ind1,11] == max(TW[ind1,11])) #ind1所记录的都是TW[,8]或P2的最大值所在的位置,在ind1中找TW[,11]或W2的最大值在ind1中的位置
    if(length(ind2) == 1){ #若在ind1中，TW[,11]或W2的最大值只在一个位置出现，则选取它
      ind = ind1[ind2]
    }
    if(length(ind2) > 1){ #若在ind1中，TW[,11]或W2的最大值在多个位置出现，则在这些位置中找到TW[,7]或|i-j|最大的那个位置
      ind3 = which.max(TW[ind1[ind2],7])
      ind = ind1[ind2][ind3]
    }
    
    #将其全部设为零，更新之后从头再来
    TW[ind1,8] = 0
    ind11 = which(TW[,8] == max(TW[,8])) #TW[,8]是P2; ind11所记录的都是更新之后TW[,8]或P2的最大值所在的位置
    if(length(ind11) > 1){#如果TW[,8]或P2有不止一个相等的最大值
      ind22 = which(TW[ind11,11] == max(TW[ind11,11])) #ind11所记录的都是更新之后TW[,8]或P2的最大值所在的位置,在ind11中找TW[,11]或W2的最大值在ind11中的位置
      if(length(ind22) == 1){ #若在ind11中，TW[,11]或W2的最大值只在一个位置出现，则选取它
        ind0 = ind11[ind22]
      }
      if(length(ind22) > 1){ #若在ind11中，TW[,11]或W2的最大值在多个位置出现，则在这些位置中找到TW[,7]或|i-j|最大的那个位置
        ind33 = which.max(TW[ind11[ind22],7])
        ind0 = ind11[ind22][ind33]
      }
    }
    if(length(ind11) == 1){
      ind0 = ind11
    }
    
    #比较ind和ind0各自对应的TW[,11]或W2的大小
    if(TW[ind,11] >= TW[ind0,11]){
      ind = ind
      TW[ind,8] = 1 #将之前设为0的TW[ind,8]设置回原值1
    }
    if(TW[ind,11] < TW[ind0,11]){
      ind = ind0
    }
    #在TW[,8]或P2的多个相等的最大值等于1的前提下,将之前设为0的TW[ind1,8]设置回原值1
    TW[ind1,8] = 1
  }
  
  #####################################################
  #若TW[,8]或P2的多个相等的最大值不等于1
  if(TW[ind1,8][1] != 1){
    ind2 = which(TW[ind1,11] == max(TW[ind1,11])) #ind1所记录的都是TW[,8]或P2的最大值所在的位置,在ind1中找TW[,11]或W2的最大值在ind1中的位置
    if(length(ind2) == 1){ #若在ind1中，TW[,11]或W2的最大值只在一个位置出现，则选取它
      ind = ind1[ind2]
    }
    if(length(ind2) > 1){ #若在ind1中，TW[,11]或W2的最大值在多个位置出现，则在这些位置中找到TW[,7]或|i-j|最大的那个位置
      ind3 = which.max(TW[ind1[ind2],7])
      ind = ind1[ind2][ind3]
    }
  }
  #ind2 = which(TW[ind1,11] == max(TW[ind1,11])) #TW[,11]是W2
  ##ind2 = which.max(TW[ind1,11]) #若TW[ind1,11]的最大值不止一个，则只取其中一个进入后续的运算
  ##ind = ind1[ind2]
}

ind
P2 = TW$P2[ind] #(Sj-vij)/Sj, Sj is the value of the second peak
W2 = TW$W2[ind] #(Sj-vij)/(Si-vij)

thre1 = 0.5
thre2 = 10
thre3 = 3 #在调查了16个模型的17条曲线（roadNetwork模型有country和city两条曲线）之后确定的
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

if(peakNum != 2){
  tk=0
  tpN=0
  for(k in 1:nrow(TW)){
    P2 = TW$P2[k] #(Sj-vij)/Sj, Sj is the value of the second peak
    W2 = TW$W2[k] #(Sj-vij)/(Si-vij)
    if(P2 >= thre1) {
      if((TW$`|j-i|`[k] >= thre3*GT) && (W2 >= (1/thre2))){pN = 2}
      else {
        if((TW$`|j-i|`[k] < thre3*GT) && (W2 >= (1/thre2))) {pN = 1.8} #time gap between two peaks is not long enough to generate a second generation
        else {pN = 1} #只要W2 < (1/thre2)，那么就认为是单峰
      }
    }
    if((P2 < thre1)&(P2>0)){
      if((TW$`|j-i|`[k] >= thre3*GT) && (W2 >= (1/thre2))) {pN = 1.5} #the valley is not deep enough to be looked as an interwave period, but the time gap is enough
      else {
        if((TW$`|j-i|`[k] < thre3*GT) && (W2 >= (1/thre2))) {pN = 1.2}
        else {pN = 1} #只要W2 < (1/thre2)，那么就认为是单峰
      }
    }
    if(P2==0) {pN = 1}#只要P2==0，那么也认为是单峰
    
    if(pN > peakNum) {
      tk = cbind(tk,k)
      tpN = cbind(tpN,pN)
    }
  }
  
  if((length(tk) > 1) && (length(tpN) > 1)){
    #tk = tk[,-1]
    #tpN = tpN[,-1]
    peakNum = max(tpN)
    ind44=which(tpN == max(tpN)) #将tpN最大值出现的所有位置都放在ind44中
    if(length(ind44)==1) { #若tpN最大值只在一个位置出现
      ind = tk[ind44]
    }
    if(length(ind44) > 1) { #若tpN最大值在多个位置出现
      ind55 = which.max(TW$W2[ tk[ind44] ])
      ind = tk[ind44][ind55]
    }
    #ind = tk[which.max(tpN)]
  }
}
print("ind=")
ind
print("peakNum=")
peakNum
