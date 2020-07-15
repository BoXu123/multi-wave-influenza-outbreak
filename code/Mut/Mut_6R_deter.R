rm(list = ls())
library(fitR)
##########Mut_6R_deter
##################model constrution
data("SEIT4L_deter")
Mut_6R_deter <- SEIT4L_deter
Mut_6R_deter$name <- c("deterministic Mut_6R model with daily incidence and constant population size")
Mut_6R_deter$theta.names <- c("R0","sigma","D_lat","D_inf","D_imm","rho","Tmut")
Mut_6R_deter$state.names  <- c("S","E1","I1","T11","T12","T13","T14","T15","T16","L1","E2","I2","T21","T22","T23","T24",
                               "T25","T26","L12","Inc")
Mut_6R_deter$simulate <- function(theta,init.state,times) {
  
  Tmut <- theta[["Tmut"]]
  
  twoVi_2R_ode_1 <- function(time, state, theta) {
    
    # param
    beta <- theta[["R0"]]/theta[["D_inf"]]
    sigma <- theta[["sigma"]]
    epsilon <- 1/theta[["D_lat"]]
    nu <- 1/theta[["D_inf"]]
    tau <- 1/theta[["D_imm"]]
    
    # states
    S <- state[["S"]]
    E1 <- state[["E1"]]
    I1 <- state[["I1"]]
    T11 <- state[["T11"]]
    T12 <- state[["T12"]]
    T13 <- state[["T13"]]
    T14 <- state[["T14"]]
    T15 <- state[["T15"]]
    T16 <- state[["T16"]]
    L1 <- state[["L1"]]
    E2 <- state[["E2"]]
    I2 <- state[["I2"]]
    T21 <- state[["T21"]]
    T22 <- state[["T22"]]
    T23 <- state[["T23"]]
    T24 <- state[["T24"]]
    T25 <- state[["T25"]]
    T26 <- state[["T26"]]
    L12 <- state[["L12"]]
    Inc <- state[["Inc"]]
    
    N <- S + E1 + I1 + T11 + T12 + T13 + T14 + T15 + T16 + L1 + E2 + I2 + T21 + T22 + T23 + T24 + T25 + T26 + L12
    
    dS <- -beta*S*I1/N
    dE1 <- beta*S*I1/N - epsilon*E1
    dI1 <- epsilon*E1 - nu*I1
    dT11 <- nu*I1 - 6*tau*T11
    dT12 <- 6*tau*T11 - 6*tau*T12
    dT13 <- 6*tau*T12 - 6*tau*T13
    dT14 <- 6*tau*T13 - 6*tau*T14
    dT15 <- 6*tau*T14 - 6*tau*T15
    dT16 <- 6*tau*T15 - 6*tau*T16
    dL1 <- 6*tau*T16
    dE2 <- 0
    dI2 <- 0
    dT21 <- 0
    dT22 <- 0
    dT23 <- 0
    dT24 <- 0
    dT25 <- 0
    dT26 <- 0
    dL12 <- 0
    dInc <- epsilon*E1
    
    return(list(c(dS,dE1,dI1,dT11,dT12,dT13,dT14,dT15,dT16,dL1,dE2,dI2,dT21,dT22,dT23,dT24,dT25,dT26,dL12,dInc)))
  }
  
  twoVi_2R_ode_2 <- function(time, state, theta) {
    
    # param
    beta <- theta[["R0"]]/theta[["D_inf"]]
    sigma <- theta[["sigma"]]
    epsilon <- 1/theta[["D_lat"]]
    nu <- 1/theta[["D_inf"]]
    tau <- 1/theta[["D_imm"]]
    
    # states
    S <- state[["S"]]
    E1 <- state[["E1"]]
    I1 <- state[["I1"]]
    T11 <- state[["T11"]]
    T12 <- state[["T12"]]
    T13 <- state[["T13"]]
    T14 <- state[["T14"]]
    T15 <- state[["T15"]]
    T16 <- state[["T16"]]
    L1 <- state[["L1"]]
    E2 <- state[["E2"]]
    I2 <- state[["I2"]]
    T21 <- state[["T21"]]
    T22 <- state[["T22"]]
    T23 <- state[["T23"]]
    T24 <- state[["T24"]]
    T25 <- state[["T25"]]
    T26 <- state[["T26"]]
    L12 <- state[["L12"]]
    Inc <- state[["Inc"]]
    
    N <- S + E1 + I1 + T11 + T12 + T13 + T14 + T15 + T16 + L1 + E2 + I2 + T21 + T22 + T23 + T24 + T25 + T26 + L12
    
    dS <- 0
    dE1 <- 0
    dI1 <- 0
    dT11 <- 0
    dT12 <- 0
    dT13 <- 0
    dT14 <- 0
    dT15 <- 0
    dT16 <- 0
    dL1 <- -sigma*beta*L1*I2/N
    dE2 <- sigma*beta*L1*I2/N - epsilon*E2
    dI2 <- epsilon*E2 - nu*I2
    dT21 <- nu*I2 - 6*tau*T21
    dT22 <- 6*tau*T21 - 6*tau*T22
    dT23 <- 6*tau*T22 - 6*tau*T23
    dT24 <- 6*tau*T23 - 6*tau*T24
    dT25 <- 6*tau*T24 - 6*tau*T25
    dT26 <- 6*tau*T25 - 6*tau*T26
    dL12 <- 6*tau*T26
    dInc <- epsilon*E2
    
    return(list(c(dS,dE1,dI1,dT11,dT12,dT13,dT14,dT15,dT16,dL1,dE2,dI2,dT21,dT22,dT23,dT24,dT25,dT26,dL12,dInc)))
    
  }
  
  #when the time before Tmut
  if(Tmut < max(times)) {
    TmutBefore = times[1:Tmut]
  }
  else { print("Wrong, Tmut is larger than max(times).")  }
  
  init.state1 <- init.state
  # put incidence at 0 in init.state
  init.state1["Inc"] <- 0
  trajBefore <- as.data.frame(ode(init.state1, TmutBefore, twoVi_2R_ode_1, theta, method = "ode45"))
  # compute incidence of each time interval
  trajBefore$Inc <- c(0, diff(trajBefore$Inc))
  
  #when the time after Tmut
  if(Tmut < max(times)){
    TmutAfter = times[(Tmut+1):max(times)]
  }
  else { print("Wrong, Tmut is larger than max(times).")  }
  
  #when the time comes to Tmut, all the compartments of virus 1 come to zero, the number of people susceptible('S') to and 
  #in the progression of temporary protection('T') against virus 1 at time of Tmut become a part of 'L1', in which the 
  #people become susceptible to virus 2. The number of people getting infection but not infectious with virus 1('E1') 
  #becomes a part of 'E2'. The initial value of 'I2' and 'Inc' both equal to the number of infectious people with 
  #virus 1 ('I1') at the time of Tmut.
  init.state2 <- c('S' = 0, 'E1' = 0, 'I1' = 0, 'T11' = 0, 'T12' = 0,'T13' = 0, 'T14' = 0,'T15' = 0, 'T16' = 0, 'L1'= trajBefore$L1[Tmut]+trajBefore$S[Tmut]+trajBefore$T11[Tmut]+trajBefore$T12[Tmut]+trajBefore$T13[Tmut]+trajBefore$T14[Tmut]+trajBefore$T15[Tmut]+trajBefore$T16[Tmut],
                  'E2'=trajBefore$E2[Tmut]+trajBefore$E1[Tmut], 'I2'=trajBefore$I1[Tmut], 'T21' = trajBefore$T21[Tmut], 'T22' = trajBefore$T22[Tmut],'T23' = trajBefore$T23[Tmut], 'T24' = trajBefore$T24[Tmut],'T25' = trajBefore$T25[Tmut], 'T26' = trajBefore$T26[Tmut], 'L12' = trajBefore$L12[Tmut], 'Inc' = trajBefore$I1[Tmut]) 
  # put incidence at 0 in init.state
  #init.state2["Inc"] <- 0
  trajAfter <- as.data.frame(ode(init.state2, TmutAfter, twoVi_2R_ode_2, theta, method = "ode45"))
  # compute incidence of each time interval
  trajAfter$Inc <- c(init.state2["Inc"], diff(trajAfter$Inc))
  
  #combine the two data frames
  traj <- rbind(trajBefore, trajAfter)
  return(traj)
  
}

Mut_6R_deter$rPointObs

#Mut_6R_deter$state.names  <- c("S","E1","I1","T11","T12","T13","T14","T15","T16","L1","E2","I2","T21","T22","T23","T24",
                              # "T25","T26","L12","Inc")
#Mut_6R_deter$theta.names <- c("R0","sigma","D_lat","D_inf","D_imm","rho","Tmut")
init.state <- c('S' = 27700, 'E1' = 0, 'I1' = 1, 'T11' = 0, 'T12' = 0,'T13' = 0, 'T14' = 0,'T15' = 0, 'T16' = 0, 'L1'= 5,
                'E2'=0, 'I2'=0, 'T21' = 0, 'T22' = 0,'T23' = 0, 'T24' = 0,'T25' = 0, 'T26' = 0, 'L12' = 0, 'Inc' = 0)
init.theta <- c(R0 = 9.2, sigma=0.18, D_lat = 2.4, D_inf = 1, D_imm = 1, rho = 0.69, Tmut =20)
#plotFit(Mut_6R_deter, theta = init.theta, init.state = init.state, data = FluTdC1971, n.replicates = 1, summary =TRUE, all.vars = TRUE)
#plotFit(Mut_6R_deter, theta = init.theta, init.state = init.state, data = FluTdC1971, n.replicates = 100, summary =TRUE, all.vars = FALSE)
#generate a whole trajectory with simulated observation. The observations are results of random draws from the model trajectory
epi_time = c(1:100)
obs.traj <- rTrajObs(fitmodel = Mut_6R_deter, theta = init.theta, init.state = init.state, times = epi_time)
#simulation <- Mut_6R_deter$simulate(theta = init.theta, init.state = init.state, times = epi_time)
#head(obs.traj)
Nreal = init.state[["S"]] + (obs.traj$L1[init.theta[["Tmut"]]] + obs.traj$S[init.theta[["Tmut"]]] + obs.traj$T11[init.theta[["Tmut"]]] + 
                               obs.traj$T12[init.theta[["Tmut"]]] + obs.traj$T13[init.theta[["Tmut"]]] + obs.traj$T14[init.theta[["Tmut"]]] + 
                               obs.traj$T15[init.theta[["Tmut"]]] + obs.traj$T16[init.theta[["Tmut"]]]) # 关注init.state2中的‘L1’
AT = (sum(obs.traj$Inc) - obs.traj$I2[init.theta[["Tmut"]] + 1])/Nreal #在sum(obs.traj$Inc)中，Tmut时刻的I1和(Tmut+1)时刻的I2重复计算了


#####################################################
### plot designed by Jun Cai
#install.packages("tidyverse")
library(tidyverse)

pdata <- obs.traj %>%
  select(time, Inc:obs) %>%
  gather(var, val, Inc:obs) %>%
  mutate(var = factor(var))

library(ggplot2)
#install.packages("ggsci")
library(ggsci)
#install.packages("ggthemes")
library(ggthemes)

pdf(file = "Mut_6R.pdf", width = 7.5, height = 5)

p <- ggplot(pdata, aes(x = time, y = val, group = var), size = 2)
p + geom_line(aes(color = var, linetype = var)) + 
  # geom_point(aes(shape = var)) + 
  scale_color_manual(values = c("red",  "black"), 
                     labels = c("Incidence", "Incidence observation")) + 
  # scale_color_npg(values = c("red", "blue", "green", "black"), 
  #                 labels = c("Incidence", "Incidence of virus 1", "Incidence of virus 2", "Incidence observation")) + 
  scale_linetype_manual(values = c("solid", "solid"), 
                        labels = c("Incidence", "Incidence observation")) + 
  labs(x = "Day", y = "Number of cases") + 
  theme_classic() + 
  theme(legend.position = c(0.8, 0.8), 
        legend.title = element_blank())

dev.off()

##################################################

#plotTraj(obs.traj)
plot(x = epi_time, y = obs.traj$Inc, type = "l", lty=2 ,col = "black", lwd=3, xlab = "time (days)", ylab = "Incidence from model")#Inc is incidence of virus 1 and 2 calculated by model
#lines(x=epi_time, y=obs.traj$Inc1, type = "l",lty=1, col="blue", lwd=1 )#Inc1 is incidence of virus 1 calculated by model
#lines(x=epi_time, y=obs.traj$Inc2, type = "l", lty=1, col="green", lwd=1 )#Inc2 is incidence of virus 2 calculated by model
lines(x = epi_time, y = obs.traj$obs, type = "l", lty=1,col = "red",lwd=1)#obs is from Inc considering reporting rate and poisson observation process
abline(v = init.theta[["Tmut"]], lty=2, col="blue",lwd=0.5)
legend("topright", legend = c("incidence","observation of incidence","time when mutation occured"),lty = c(2,1,2), lwd = c(3,1,0.5),col = c("black","red","blue") )
##################################################

##################################################################################################################
##判定双峰
##########################方案3
#GT = init.theta[["D_lat"]] + init.theta[["D_inf"]]#generation time
GT = init.theta[["D_lat"]] + 0.5*init.theta[["D_inf"]]#generation time
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
  #combine
  TW = rbind(TWl,TWr)
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
thre2 = 10 #max(W2)=0.1311535247，所以thre2在这里要大于1/0.1311535247=7.62才能保证peakNum为2
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
