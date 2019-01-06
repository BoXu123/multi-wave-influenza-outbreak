rm(list = ls())
library(fitR)
##########AoN_4R_deter
##################model constrution
data("SEIT4L_deter")
AoN_4R_deter <- SEIT4L_deter
AoN_4R_deter$state.names <- c("S","E","I","R1","R2","R3","R4","L","RS","Inc")
AoN_4R_deter$theta.names <- c("R0","D_lat","D_inf","D_imm","alpha","rho")
AoN_4R_deter$simulate <- function(theta,init.state,times) {
  
  AoN_4R_ode <- function(time, state, theta) {
    
    # param
    beta <- theta[["R0"]]/theta[["D_inf"]]
    epsilon <- 1/theta[["D_lat"]]
    nu <- 1/theta[["D_inf"]]
    alpha <- theta[["alpha"]]
    gamma <- 1/theta[["D_imm"]]
    
    # states
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    R1 <- state[["R1"]]
    R2 <- state[["R2"]]
    R3 <- state[["R3"]]
    R4 <- state[["R4"]]
    L <- state[["L"]]
    Inc <- state[["Inc"]]
    RS <- state[["RS"]] #表示由R状态转变为S状态的人数
    
    N <- S + E + I + R1 + R2 + R3 + R4 + L
    
    dS <- -beta*S*I/N + (1-alpha)*4*gamma*R4
    dE <- beta*S*I/N - epsilon*E
    dI <- epsilon*E - nu*I
    dR1 <- nu*I - 4*gamma*R1
    dR2 <- 4*gamma*R1 - 4*gamma*R2
    dR3 <- 4*gamma*R2 - 4*gamma*R3
    dR4 <- 4*gamma*R3 - 4*gamma*R4
    dL <- alpha*4*gamma*R4
    dInc <- epsilon*E
    dRS <- (1-alpha)*4*gamma*R4
    
    return(list(c(dS,dE,dI,dR1,dR2,dR3,dR4,dL,dRS,dInc)))
  }
  
  
  # put incidence at 0 in init.state
  init.state["Inc"] <- 0
  init.state["RS"] <- 0
  
  traj <- as.data.frame(ode(init.state, times, AoN_4R_ode, theta, method = "ode45"))
  
  # compute incidence of each time interval
  traj$Inc <- c(0, diff(traj$Inc))
  traj$RS <- c(0, diff(traj$RS))
  
  return(traj)
}

#generates a single random observation from a single point
#in a model trajectory
AoN_4R_deter$rPointObs <- function(model.point, theta){
  
  obs.point<-rpois(n=1,lambda=theta[["rho"]]*model.point[["Inc"]])
  
  return(c(obs=obs.point))
}

# model simulation
init.state <- c('S' = 27700, 'E' = 0, 'I' = 1, 
                'R1' = 0, 'R2' = 0,'R3' = 0, 
                'R4' = 0, 'L'= 6,'RS'=0,'Inc' = 0)
init.theta <- c(R0 = 11.2, D_lat = 2.1, D_inf = 2.4, 
                D_imm = 11.6, alpha = 0.42, rho = 0.71)
epi_time = c(1:100)
#epi_time = c(1:365)
obs.traj <- rTrajObs(fitmodel = AoN_4R_deter, theta = init.theta, 
                     init.state = init.state, times = epi_time)

##clinical attack rate
#RS表示由R状态转变为S状态的人数，这部分人再次成为病毒的易感者
Nreal = init.state[["S"]]+sum(obs.traj$RS)
AT = sum(obs.traj$Inc)/Nreal

#head(obs.traj)
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

pdf(file = "AoN_4R.pdf", width = 7.5, height = 5)

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
##plot model-derived epidemic curve
#plotTraj(obs.traj)
#Inc is incidence of virus calculated by model
plot(x = epi_time, y = obs.traj$Inc, type = "l", 
     lty=1 ,col = "black", lwd=3, xlab = "Day", 
     ylab = "Incidence from model")
#obs is from Inc considering reporting rate and 
#poisson observation process
lines(x = epi_time, y = obs.traj$obs, type = "l", 
      lty=1,col = "red",lwd=1)
legend("topright", legend = c("incidence","observation of incidence"),
       lty = c(1,1), lwd = c(3,1),col = c("black","red") )

##################################################################################################################

##判定双峰
##########################方案3
#GT = init.theta[["D_lat"]] + init.theta[["D_inf"]]#generation time
#GT = init.theta[["D_lat"]] + init.theta[["D_inf"]] + init.theta[["D_imm"]]#generation time
GT = init.theta[["D_lat"]] + 0.5*init.theta[["D_inf"]]
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
thre2 = 10 #W2最大值是0.1032586261，则thre2要大于1/0.1032586261=9.68才能保证peakNum为2
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
