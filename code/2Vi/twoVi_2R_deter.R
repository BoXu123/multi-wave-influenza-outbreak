rm(list = ls())
library(fitR)
###twoVi_2R_deter
######################################model construction
data("SEITL_deter")
twoVi_2R_deter <- SEITL_deter
twoVi_2R_deter$state.names <- c("S","E1","I1","T11","T12","L1","E2","I2","T21","T22","L12","Inc1","Inc2","Inc","IncL1")
twoVi_2R_deter$theta.names <- c("R01","R02","D_lat","D_inf","D_imm","rho" )
twoVi_2R_deter$name <- c("deterministic 2-virus model with daily incidence and constant population size")
twoVi_2R_deter$simulate <- function(theta,init.state,times) {
  
  twoVi_2R_ode <- function(time, state, theta) {
    
    # param
    beta1 <- theta[["R01"]]/theta[["D_inf"]]
    beta2 <- theta[["R02"]]/theta[["D_inf"]]
    epsilon <- 1/theta[["D_lat"]]
    nu <- 1/theta[["D_inf"]]
    tau <- 1/theta[["D_imm"]]
    
    # states
    S <- state[["S"]]
    E1 <- state[["E1"]]
    I1 <- state[["I1"]]
    T11 <- state[["T11"]]
    T12 <- state[["T12"]]
    L1 <- state[["L1"]]
    E2 <- state[["E2"]]
    I2 <- state[["I2"]]
    T21 <- state[["T21"]]
    T22 <- state[["T22"]]
    L12 <- state[["L12"]]
    Inc1 <- state[["Inc1"]]
    Inc2 <- state[["Inc2"]]
    Inc <- state[["Inc"]]
    IncL1 <- state[["IncL1"]]
    
    N <- S + E1 + I1 + T11 + T12 + L1 + E2 + I2 + T21 + T22 + L12
    
    dS <- -beta1*S*I1/N
    dE1 <- beta1*S*I1/N - epsilon*E1
    dI1 <- epsilon*E1 - nu*I1
    dT11 <- nu*I1 - 2*tau*T11
    dT12 <- 2*tau*T11 - 2*tau*T12
    dL1 <- 2*tau*T12 - beta2*L1*I2/N
    dE2 <- beta2*L1*I2/N - epsilon*E2
    dI2 <- epsilon*E2 - nu*I2
    dT21 <- nu*I2 - 2*tau*T21
    dT22 <- 2*tau*T21 - 2*tau*T22
    dL12 <- 2*tau*T22
    dInc1 <- epsilon*E1
    dInc2 <- epsilon*E2
    dInc <- epsilon*E1 + epsilon*E2
    dIncL1 <- 2*tau*T12
    
    return(list(c(dS,dE1,dI1,dT11,dT12,dL1,dE2,dI2,dT21,dT22,dL12,dInc1,dInc2,dInc,dIncL1)))
  }
  
  
  # put incidence at 0 in init.state
  init.state["Inc1"] <- 0
  init.state["Inc2"] <- 0
  init.state["Inc"] <- 0
  init.state["IncL1"] <- 0
  
  traj <- as.data.frame(ode(init.state, times, twoVi_2R_ode, theta, method = "ode45"))
  
  # compute incidence of each time interval
  traj$Inc1 <- c(0, diff(traj$Inc1))
  traj$Inc2 <- c(0, diff(traj$Inc2))
  traj$Inc <- c(0, diff(traj$Inc))
  traj$IncL1 <- c(0, diff(traj$IncL1))
  
  return(traj)
  
}

twoVi_2R_deter$rPointObs <- function(model.point, theta){
  
  obs.point <- rpois(n=1, lambda=theta[["rho"]]*model.point[["Inc"]])
  
  return(c(obs=obs.point))
}

twoVi_2R_deter$dprior <- function(theta, log = FALSE) {
  
  # package with truncated normal distribution
  library(truncnorm)
  
  log.prior.R01 <- dunif(theta[["R01"]], min = 1, max = 50, log = TRUE)
  log.prior.R02 <- dunif(theta[["R02"]], min = 1, max = 50, log = TRUE)
  # normal distribution with mean = 2 and sd = 1 and truncated at 0
  log.prior.latent.period <- log(dtruncnorm(theta[["D_lat"]], a = 0, b = Inf, 
                                            mean = 2, sd = 1))
  # normal distribution with mean = 2 and sd = 1 and truncated at 0
  log.prior.infectious.period <- log(dtruncnorm(theta[["D_inf"]], a = 0, b = Inf, 
                                                mean = 2, sd = 1))
  log.prior.temporary.immune.period <- dunif(theta[["D_imm"]], min = 0, max = 50, log = TRUE)
  log.prior.reporting.rate <- dunif(theta[["rho"]], min = 0, max = 1, log = TRUE)
  
  log.sum = log.prior.R01 + log.prior.R02 + log.prior.latent.period + log.prior.infectious.period + log.prior.temporary.immune.period  + log.prior.reporting.rate
  
  return(ifelse(log, log.sum, exp(log.sum)))
  
}

twoVi_2R_deter$dPointObs <- function(data.point, model.point, theta, log = FALSE){
  
  return(dpois(x=data.point[["obs"]],lambda=theta[["rho"]]*model.point[["Inc"]],log=log))
  
}
##################################################
################################plot curve
init.state <- c('S' = 27700, 'E1' = 0, 'I1' = 20, 'T11' = 0, 'T12' = 0, 'L1'= 3, 'E2'=0, 'I2'=10, 'T21' = 0, 'T22' = 0, 'L12' = 0, 'Inc1' = 0, 'Inc2' = 0, 'Inc' = 0, 'IncL1' = 0)
#init.theta <- c(R01 = 9.4015, R02=37.8601, D_lat = 2.0412, D_inf = 3.9359, D_imm = 20.5047, rho = 0.5846)
init.theta <- c(R01 = 4.2, R02=4, D_lat = 1.5, D_inf = 2, D_imm = 0.3, rho = 0.65)
#init.theta <- c(R01 = 1.9, R02=1.8, D_lat = 2.2, D_inf = 2.3, D_imm = 5.5, rho = 0.65)
#plotFit(twoVi_2R_deter, theta = init.theta, init.state = init.state, data = FluTdC1971, n.replicates = 1, summary =TRUE)
#generate a whole trajectory with simulated observation. The observations are results of random draws from the model trajectory
epi_time = c(1:100)
#epi_time = c(1:365)
obs.traj <- rTrajObs(fitmodel = twoVi_2R_deter, theta = init.theta, init.state = init.state, times = epi_time)
#head(obs.traj)
#Nreal = init.state[["S"]]+init.state[["E1"]]+init.state[["I1"]]+init.state[["T11"]]+init.state[["T12"]]+init.state[["L1"]]+init.state[["E2"]]+init.state[["I2"]]+init.state[["T21"]]+init.state[["T22"]]+init.state[["L12"]]
Nreal = init.state[["S"]]+sum(obs.traj$IncL1) #IncL1表示L1的增量，这部分人又成为virus 2的易感者
AT = sum(obs.traj$Inc)/Nreal

#####################################################
### plot designed by Jun Cai
#install.packages("tidyverse")
library(tidyverse)

pdata <- obs.traj %>%
  select(time, Inc1:obs) %>%
  gather(var, val, Inc1:obs) %>%
  mutate(var = factor(var))

library(ggplot2)
#install.packages("ggsci")
library(ggsci)
#install.packages("ggthemes")
library(ggthemes)

pdf(file = "twoVi_2R.pdf", width = 7.5, height = 5)

p <- ggplot(pdata, aes(x = time, y = val, group = var), size = 2)
p + geom_line(aes(color = var, linetype = var)) + 
  # geom_point(aes(shape = var)) + 
  scale_color_manual(values = c("red", "gray", "gray", "black"), 
                     labels = c("Incidence", "Incidence of virus 1", "Incidence of virus 2", "Incidence observation")) + 
  # scale_color_npg(values = c("red", "blue", "green", "black"), 
  #                 labels = c("Incidence", "Incidence of virus 1", "Incidence of virus 2", "Incidence observation")) + 
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "solid"), 
                        labels = c("Incidence", "Incidence of virus 1", "Incidence of virus 2", "Incidence observation")) + 
  labs(x = "Day", y = "Number of cases") + 
  theme_classic() + 
  theme(legend.position = c(0.8, 0.8), 
        legend.title = element_blank())

dev.off()

##################################################

### plot designed by Bo Xu
# plotTraj(obs.traj)
plot(x = epi_time, y = obs.traj$Inc, type = "l", lty = 2, col = "black", lwd = 3, 
     xlab = "time (days)", ylab = "Incidence from model") # Inc is incidence of virus 1 and 2 calculated by model
lines(x = epi_time, y = obs.traj$Inc1, type = "l", lty = 1, col = "blue", lwd = 1 ) # Inc1 is incidence of virus 1 calculated by model
lines(x = epi_time, y = obs.traj$Inc2, type = "l", lty = 1, col = "green", lwd = 1 ) # Inc2 is incidence of virus 2 calculated by model
lines(x = epi_time, y = obs.traj$obs, type = "l", lty = 1, col = "red", lwd = 1) # obs is from Inc considering reporting rate and poisson observation process
legend("topright", legend = c("incidence of both", "incidence of virus 1", "incidence of virus 2","observation of both"), lty = c(2, 1, 1, 1), lwd = c(3, 1, 1, 1), 
       col = c("black", "blue", "green", "red") )
##################################################

##################################################################################################################
##判定双峰
##########################方案3
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
}

#combine
TW = rbind(TWl,TWr)
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

##############################MCMC parameters inference
my_dLogPosterior <- function(fitmodel, theta, init.state, data) {
  
  # calculate the fitmodel prior for parameter vector theta using
  # fitmodel$dprior, and assign to variable log.prior
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  
  # calculate the log-likelihood of `theta`
  # and `init.state` with respect to the data using `dTrajObs`
  # and assign to a variable `log.likelihood`    
  log.likelihood <- dTrajObs(fitmodel, theta, init.state, data, log = TRUE)
  
  # calulate the log-posterior using the log-prior and log-likelihood
  log.posterior <- log.prior + log.likelihood
  
  return(log.posterior)
  
}

my_dLogPosterior_epi3 <- function(theta) {
  
  return(my_dLogPosterior(fitmodel = SIR,
                          theta = theta,
                          init.state = c(S = 999, I = 1, R = 0),
                          data = epi3))
}

# wrapper for posterior
my_posteriorTdC <- function(theta) {
  
  my_fitmodel <- twoVi_2R_deter
  #my_init.state <- c(S = 279, E = 0, I = 2, T = 3, L = 0, Inc = 0)
  # note that for the SEIT4L model there are 4 state variables for the Tcompartment 
  my_init.state <- c('S' = 277, 'E1' = 0, 'I1' = 3, 'T11' = 3, 'T12' = 0, 'L1'= 3, 'E2'=0, 'I2'=1, 'T21' = 0, 'T22' = 0, 'L12' = 0, 'Inc' = 0)
  
  return(logPosterior(fitmodel = my_fitmodel, theta = theta, init.state = my_init.state, 
                      data = FluTdC1971, margLogLike = dTrajObs, log = TRUE))#margLogLike expects a function that returns
  #the log-likelihood, which can be obtained by passing log = TRUE to dTrajObs. To do so, we can use the dot-dot-dot 
  #argument of logPosterior, which allows you to pass any extra argument to the function assigned to margLogLike.
  
}
# theta to initialise the MCMC
#twoVi_2R_deter$theta.names <- c("R01","R02","D_lat","D_inf","D_imm","rho" )
init.theta <- c(R01 = 2, R02=13, D_lat = 2, D_inf = 2, D_imm = 1, rho = 0.5)

# diagonal elements of the covariance matrix for the Gaussian proposal
proposal.sd <- c(R01 = 0.5, R02=0.5, D_lat = 0.5, D_inf = 0.5,  D_imm = 0.5, rho = 0.1)

# lower and upper limits of each parameter
lower <- c(R01 = 1, R02=1, D_lat = 0.1, D_inf = 0.1, D_imm = 0.1, rho = 0)
upper <- c(R01 = 50, R02=50, D_lat = Inf, D_inf = Inf, D_imm = Inf, rho = 1)

# number of iterations for the MCMC
n.iterations <- 20000

# additional parameters for the adaptive MCMC, see ?mcmcMH for more details
adapt.size.start <- 100
adapt.size.cooling <- 0.999
adapt.shape.start <- 200
trace <- mcmcMH(target = my_posteriorTdC,
                init.theta = init.theta,
                proposal.sd = proposal.sd,
                n.iterations = n.iterations,
                adapt.size.start = adapt.size.start,
                adapt.shape.start = adapt.shape.start,
                adapt.size.cooling= adapt.size.cooling,
                limits = list(lower = lower, upper= upper))
########################################
###############################Diagnostic
my_mcmc.TdC <- trace
library('coda')
# convert to a mcmc object for coda
my_trace <- mcmc(my_mcmc.TdC$trace)
# compute the acceptance rate
1 - rejectionRate(my_trace)
# between 0.1 and 0.6: looks good!

# plot the trace
library("lattice")  ## for xyplot
xyplot(my_trace)
# Let's find a suitable burning:
plotESSBurn(my_trace)
effectiveSize(my_trace)
trace.burn <- burnAndThin(my_trace, burn = 1250)
xyplot(x = trace.burn)
effectiveSize(trace.burn)
acfplot(x = trace.burn, lag.max = 60)
# Let's create a thinned trace
trace.burn.thin <- burnAndThin(trace.burn, thin = 20)
xyplot(x = trace.burn.thin)
# Let's check the ESS
effectiveSize(trace.burn.thin)
#Although the thinned trace has 20 times less fewer than the unthinned trace, it has a similar ESS. This is because the autocorrelation has been reduced.
# new autocorrelation plot
acfplot(x = trace.burn.thin, lag.max = 60)
# density plot
plotPosteriorDensity(list(trace.burn.thin))
#or
densityplot(x = trace.burn.thin)

#compare the posterior estimates of the thinned and unthinned traces
summary(trace.burn)
summary(trace.burn.thin)

#They are very similar. So why thin? Because autocorrelation produces clumpy samples that are unrepresentative, 
#in the short run, of the true underlying posterior distribution. We can check this by comparing the thinned and 
#unthinned distributions using the function plotPosteriorDensity of the fitR package
plotPosteriorDensity(list(unthinned = trace.burn, thinned = trace.burn.thin))

######################################plot curves with parameters inferenced
# the same init.state as for the fit
init.state <- c('S' = 277, 'E1' = 0, 'I1' = 3, 'T11' = 3, 'T12' = 0, 'L1'= 3, 'E2'=0, 'I2'=1, 'T21' = 0, 'T22' = 0, 'L12' = 0, 'Inc' = 0)

# by default plotPosteriorFit summarize the fit of 100 thetas sampled from
# the posterior
plotPosteriorFit(trace = trace.burn.thin, fitmodel = twoVi_2R_deter, init.state = init.state, 
                 data = FluTdC1971)
# alternatively, one can plot the fit of the mean of the posterior (in this
# case the observation is replicated 100 times)
plotPosteriorFit(trace = trace.burn.thin, fitmodel = twoVi_2R_deter, init.state = init.state, 
                 data = FluTdC1971, posterior.summary = "mean")
# or using the maximum a posteriori (MAP) estimate
plotPosteriorFit(trace = trace.burn.thin, fitmodel = twoVi_2R_deter, init.state = init.state, 
                 data = FluTdC1971, posterior.summary = "max")
##########################################################################

######################################calculate DIC
trace.combined <- as.data.frame(trace.burn.thin)
# take the mean of theta
theta.bar <- colMeans(trace.combined[twoVi_2R_deter$theta.names])#加了[SEITL_deter$theta.names]就只有前6个参数出来
#print(colMeans(trace.combined))#不加[SEITL_deter$theta.names]就会有9个参数全出来
print(theta.bar)
##         R01        R02      D_lat      D_inf      D_imm        rho 
##9.4015169 37.8600603  2.0411560  3.9359134 20.5046531  0.5846181

# compute its log-likelihood
init.state <- c('S' = 277, 'E1' = 0, 'I1' = 3, 'T11' = 3, 'T12' = 0, 'L1'= 3, 'E2'=0, 'I2'=1, 'T21' = 0, 'T22' = 0, 'L12' = 0, 'Inc' = 0)
log.like.theta.bar <- dTrajObs(twoVi_2R_deter, theta.bar, init.state, data = FluTdC1971, 
                               log = TRUE)
print(log.like.theta.bar)
## [1] -146.2264

# and its deviance
D.theta.bar <- -2 * log.like.theta.bar
print(D.theta.bar)
## [1] 292.4528

# the effective number of parameters
p.D <- var(-2 * trace.combined$log.likelihood)/2
print(p.D)
## [1] 7.867571

# and finally the DIC
DIC <- D.theta.bar + 2 * p.D
print(DIC)
## [1] 308.1879