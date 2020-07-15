rm(list=ls())
library(fitR)
##########TCC_RRC_deter (SEIR model with transmission coefficient and reporting rate changing along time)
##################model constrution
data("SEITL_deter")
TCC_RRC_deter <- SEITL_deter
# nu is birth and natural death rate. beta is transmission rate. q is a reduction factor in the transmissibility of the
# asymptomatic class. 1/k is latent period. rho is the proportion of latent individuals who progress to the clinically
# infectious class. gamma1 is the recovery rate of clinically infectious individuals and asymptomatic individuals. 
# alpha is the diagnostic rate. gamma2 is the recovery rate for hospitalized class. sigma is the mortality rate. RR is
# the reporting rate. T is the time point when the second wave started.
TCC_RRC_deter$theta.names <- c("nu","beta1","q1","k1","rho1","gamma11","alpha1","gamma21","sigma1","RR1",
                               "beta2","q2","k2","rho2","gamma12","alpha2","gamma22","sigma2","RR2","T")
TCC_RRC_deter$state.names <- c("S1","E1","I1","A1","J1","R1","D1","Hos1","S2","E2","I2","A2","J2","R2","D2","Hos2")
TCC_RRC_deter$simulate <- function(theta,init.state,times) {
  
  TCC_RRC_ode <- function(time, state, theta) {
    
    # param
    nu <- theta[["nu"]] 
    beta <- theta[["beta"]]
    q <- theta[["q"]]
    k <- theta[["k"]]
    rho <- theta[["rho"]]
    gamma1 <- theta[["gamma1"]]
    alpha <- theta[["alpha"]]
    gamma2 <- theta[["gamma2"]]
    sigma <- theta[["sigma"]]
    RR <- theta[["RR"]]
    
    # states
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    A <- state[["A"]]
    J <- state[["J"]]
    R <- state[["R"]]
    D <- state[["D"]]
    Hos <- state[["Hos"]]
    
    N <- S + E +I + A + J + R
    
    dS <- nu*N - beta*S*(I+J+q*A)/N - nu*S
    dE <- beta*S*(I+J+q*A)/N - (k+nu)*E
    dI <- k*rho*E - (alpha+gamma1+nu)*I
    dA <- k*(1-rho)*E - (gamma1+nu)*A
    dJ <- alpha*I - (gamma2+sigma+nu)*J
    dR <- gamma1*(A+I) + gamma2*J - nu*R
    dD <- sigma*J
    dHos <- alpha*I
    
    return(list(c(dS,dE,dI,dA,dJ,dR,dD,dHos)))
  }
  
  T <- theta[["T"]]
  ####wave 1 (before time T)
  times1 <- c(1:T)
  theta1 <- c(nu=theta[["nu"]], beta=theta[["beta1"]], q=theta[["q1"]], k=theta[["k1"]], rho=theta[["rho1"]],
              gamma1=theta[["gamma11"]], alpha=theta[["alpha1"]], gamma2=theta[["gamma21"]], sigma=theta[["sigma1"]],
              RR=theta[["RR1"]])
  init.state1 <- c('S'=init.state[["S1"]], 'E'=init.state[["E1"]], 'I'=init.state[["I1"]], 'A'=init.state[["A1"]], 'J'=
                     init.state[["J1"]], 'R'=init.state[["R1"]], 'D'=init.state[["D1"]], 'Hos'=init.state[["Hos1"]])
  # put incidence at 0 in init.state
  init.state1["Hos"] <- 0
  traj1 <- as.data.frame(ode(init.state1, times1, TCC_RRC_ode, theta1, method = "ode45"))
  
  # compute incidence of each time interval
  traj1$Hos <- c(0, diff(traj1$Hos))
  
  ####wave 2 (after time T)
  times2 <- c((T+1):max(times))
  theta2 <- c(nu=theta[["nu"]], beta=theta[["beta2"]], q=theta[["q2"]], k=theta[["k2"]], rho=theta[["rho2"]],
              gamma1=theta[["gamma12"]], alpha=theta[["alpha2"]], gamma2=theta[["gamma22"]], sigma=theta[["sigma2"]],
              RR=theta[["RR2"]])
   
  init.state2 <- c('S'=traj1$S[T], 'E'=9, 'I'=7, 'A'=traj1$A[T], 'J'=
                     traj1$J[T], 'R'=traj1$R[T], 'D'=traj1$D[T], 'Hos'=traj1$Hos[T])
  
  init.state2["Hos"] <- init.state2["Hos"]
  traj2 <- as.data.frame(ode(init.state2, times2, TCC_RRC_ode, theta2, method = "ode45"))
  # compute incidence of each time interval
  traj2$Hos <- c(init.state2["Hos"], diff(traj2$Hos))
  
  #combine the two data frames
  traj <- rbind(traj1, traj2)
  
  return(traj)
}

TCC_RRC_deter$rPointObs <- function(model.point, theta){
  
  ## the prevalence is observed through a Poisson process
  #obs.point <- rpois(n = 1, lambda = model.point[["Hos"]] * theta[["RR"]])
  ##the prevalence is observed through a Poisson distribution process
  #if(model.point[["time"]] <= theta[["T"]]) {
  #  obs.point <- rpois(n = 1, lambda = model.point[["Hos"]] * theta[["RR1"]])
  #}
  #else {
  #  obs.point <- rpois(n = 1, lambda = model.point[["Hos"]] * theta[["RR2"]])
  #}
  
  ##the prevalence is observed through a normal distribution process
  if(model.point[["time"]] <= theta[["T"]]) {
    obs.point <- rnorm(n = 1, mean = model.point[["Hos"]] * theta[["RR1"]], sd=0.1)
  }
  else {
    obs.point <- rnorm(n = 1, mean = model.point[["Hos"]] * theta[["RR2"]], sd=0.1)
  }
  return(c(obs = obs.point))
}

TCC_RRC_deter$dPointObs <- function(data.point, model.point, theta, log = FALSE){
  
  ## the prevalence is observed through a Poisson process with a reporting rate
  #return(dpois(x = data.point[["obs"]], lambda = model.point[["Hos"]] * theta[["RR"]], log = log))
  #the prevalence is observed through a normal distribution process
  if(model.point[["time"]] <= theta[["T"]]) {
    return(dnorm(x=data.point[["obs"]], mean = model.point[["Hos"]] * theta[["RR1"]], sd=0.1, log = log))
  }
  else {
    
    return(dnorm(x=data.point[["obs"]], mean = model.point[["Hos"]] * theta[["RR2"]], sd=0.1, log = log))
    
  }
}


epi_time = c(1:180)
# simulation when time interval between onsets equals to nonzero.
#"S1","E1","I1","A1","J1","R1","D1","Hos1","S2","E2","I2","A2","J2","R2","D2","Hos2"
init.state <- c('S1' = 174334,'E1'=157,'I1' = 80,'A1'=0,'J1'=0,'R1' = 0,'D1'=0,'Hos1' = 0)
#"nu","beta1","q1","k1","rho1","gamma11","alpha1","gamma21","sigma1","RR1","beta2","q2","k2","rho2","gamma12","alpha2","gamma22","sigma2","RR2","T"
init.theta <- c(nu = 1/(60*365), beta1 = 8.0, q1 = 0.003, k1=0.53, rho1=0.10, gamma11=0.34, alpha1=0.51, gamma21=1.10,
                sigma1=0.01, RR1=0.597, beta2 = 5.75, q2 = 0.014, k2=0.53, rho2=0.36, gamma12=0.45, alpha2=2.14, 
                gamma22=0.58, sigma2=0.02, RR2=0.83, T=72) # The 1st day is Jul 01, the 72nd day is Sep 10
#simulation <- TCC_RRC_deter$simulate(theta = init.theta, init.state = init.state, times = epi_time)
library(plyr)
obs.traj <- rTrajObs(fitmodel = TCC_RRC_deter, theta = init.theta, init.state = init.state, times = epi_time)
#####################################################
### plot designed by Jun Cai
#install.packages("tidyverse")
library(tidyverse)

pdata <- obs.traj %>%
  select(time, Hos:obs) %>%
  gather(var, val, Hos:obs) %>%
  mutate(var = factor(var))

library(ggplot2)
#install.packages("ggsci")
library(ggsci)
#install.packages("ggthemes")
library(ggthemes)

pdf(file = "TCC_RRC_1918-Flu.pdf", width = 7.5, height = 5)

p <- ggplot(pdata, aes(x = time, y = val, group = var), size = 2)
p + geom_line(aes(color = var, linetype = var)) + 
  # geom_point(aes(shape = var)) + 
  scale_color_manual(values = c("red","black"), 
                     labels = c("Hospitalized", "Hospitalized observation")) + 
  # scale_color_npg(values = c("red", "blue", "green", "black"), 
  #                 labels = c("Incidence", "Incidence of virus 1", "Incidence of virus 2", "Incidence observation")) + 
  scale_linetype_manual(values = c("solid","solid"), 
                        labels = c("Hospitalized", "Hospitalized observation")) + 
  labs(x = "Day", y = "Number of inpatients") + 
  theme_classic() + 
  theme(legend.position = c(0.8, 0.8), 
        legend.title = element_blank())

dev.off()

##################################################
### plot designed by Bo Xu
par(mar=c(5, 4, 2, 4))
plot(x = epi_time, y = obs.traj$Hos,type = "l", lty=2 ,col = "black", lwd=3, xlab = "time (days)", ylab = "Number of hospitalized people per day from model")
#lines(x=epi_time, y=obs.traj$inc1, type = "l",lty=1, col="blue", lwd=1 )#Inc1 is incidence of subregion 1 calculated by model
#lines(x=c(init.theta[["T"]]:max(epi_time)), y=obs.traj$inc2[init.theta[["T"]]:max(epi_time)], type = "l", lty=1, col="green", lwd=1 )#Inc2 is incidence of subregion 2 calculated by model
lines(x = epi_time, y = obs.traj$obs, lty=1,col = "red",lwd=1)#obs is from Hos considering reporting rate and normal distribution observation process
abline(v = init.theta[["T"]], lty=2, col="blue",lwd=0.5)
legend("topright", legend = c("hospitalized people","observation of hospitalized","onset of wave 2"),bty = "n",lty = c(2,1,2), lwd = c(3,1,0.5),col = c("black","red","blue") )

##################################################################################################################
##判定双峰
##########################方案3
GT = (1/init.theta[["k1"]]) + 0.5*(1/init.theta[["gamma11"]])
obs.traj$Inc = obs.traj$Hos #这一步替换很重要
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
#ind1 = which.max(TW[,8]) #若TW[,8]的最大值不止一个（即多值重复），则which.max(TW[,8])只会返回一个位置值，
#而which(TW[,8] == max(TW[,8]))会全部返回所有位置值
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
thre2 = 5 #W2最大值为0.2423335448，thre2要大于1/0.2423335448=4.13才能使peakNum为2
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