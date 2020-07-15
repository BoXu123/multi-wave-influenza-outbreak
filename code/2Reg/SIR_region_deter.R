rm(list = ls())
library(fitR)
##########SIR_region_deter (perform two SIR models corresponding to two sub-regions respectively, and combine both along 
####time to get a single curve describing the whole situation of both. They have the same epidemiological parameters and
####epidemic onset but different initial states.)
##################model constrution
data("SIR_reporting")
SIR_region_deter <- SIR_reporting
SIR_region_deter$theta.names <- c("R0","D_inf","RR")
SIR_region_deter$state.names <- c("S1","I1","R1","Inc1","S2","I2","R2","Inc2")
SIR_region_deter$simulate <- function(theta,init.state,times) {
  
  SIR_region1 <- function(time, state, parameters) {# ��T==0����������ͬʱ��ʼ��������
    
    ## parameters
    beta <- parameters[["R0"]] / parameters[["D_inf"]]
    nu <- 1 / parameters[["D_inf"]]
    
    ## states
    S1 <- state[["S1"]]
    I1 <- state[["I1"]]
    R1 <- state[["R1"]]
    Inc1 <- state[["Inc1"]]
    
    N1 <- S1 + I1 + R1
    
    dS1 <- -beta * S1 * I1/N1
    dI1 <- beta * S1 * I1/N1 - nu * I1
    dR1 <- nu * I1
    dInc1 <- beta*S1*I1/N1
    
    ## states
    S2 <- state[["S2"]]
    I2 <- state[["I2"]]
    R2 <- state[["R2"]]
    Inc2 <- state[["Inc2"]]
    
    N2 <- S2 + I2 + R2
    
    dS2 <- -beta * S2 * I2/N2
    dI2 <- beta * S2 * I2/N2 - nu * I2
    dR2 <- nu * I2
    dInc2 <- beta*S2*I2/N2
    
    return(list(c(dS1, dI1, dR1, dInc1, dS2, dI2, dR2, dInc2)))
  }
  
  SIR_region2 <- function(time, state, parameters) {# ��T!=0�����������������п�ʼʱ�����T��
    
    ## parameters
    beta <- parameters[["R0"]] / parameters[["D_inf"]]
    nu <- 1 / parameters[["D_inf"]]
    
    ## states
    S <- state[["S"]]
    I <- state[["I"]]
    R <- state[["R"]]
    inc <- state[["inc"]]
    
    N <- S + I + R
    
    dS <- -beta * S * I/N
    dI <- beta * S * I/N - nu * I
    dR <- nu * I
    dinc <- beta*S*I/N
    
    return(list(c(dS, dI, dR, dinc)))
  }
  
  ######################use the value of T to decide which kind of calculation to be implemented
  T <- theta[["T"]] #T���������������п�ʼʱ��ļ������
  if(T==0) { # ��T==0����������ͬʱ��ʼ��������
    ###����4�У���74�е���77�У��Ǵ����
    ## put incidence at 0 in init.state
    #init.state["Inc1"] <- 0
    #init.state["Inc2"] <- 0
    #traj <- as.data.frame(ode(init.state, times, SIR_region1, theta, method = "ode45"))
    
    init.state3 <- c('S1' = init.state[["S1"]], 'I1' = init.state[["I1"]], 'R1' = init.state[["R1"]],'Inc1' = init.state[["Inc1"]],
                     'S2' = init.state[["S2"]], 'I2' = init.state[["I2"]], 'R2' = init.state[["R2"]],'Inc2' = init.state[["Inc2"]])
    # put incidence at 0 in init.state
    init.state3["Inc1"] <- 0
    init.state3["Inc2"] <- 0
    
    traj <- as.data.frame(ode(init.state3, times, SIR_region1, theta, method = "ode45"))
    
    # compute incidence of each time interval
    traj$Inc1 <- c(0, diff(traj$Inc1))
    traj$Inc2 <- c(0, diff(traj$Inc2))
    
    traj$Inc <- traj$Inc1 + traj$Inc2
    
    return(traj)
  }
  else {
    if(T>=1){ # ��T!=0�����������������п�ʼʱ�����T��
      init.state1 <- c('S' = init.state[["S1"]], 'I' = init.state[["I1"]], 'R' = init.state[["R1"]],'inc' = init.state[["Inc1"]])
      # put incidence at 0 in init.state
      init.state1["inc"] <- 0
      
      traj1 <- as.data.frame(ode(init.state1, times, SIR_region2, theta, method = "ode45"))
      # compute incidence of each time interval
      traj1$inc <- c(0, diff(traj1$inc))
      
      
      #when the time after T
      if(T < max(times)) {
        Tafter = times[T:max(times)]
      }
      else { print("Wrong, T is larger than max(times).")}
      
      init.state2 <- c('S' = init.state[["S2"]], 'I' = init.state[["I2"]], 'R' = init.state[["R2"]],'inc' = init.state[["Inc2"]])
      # put incidence at 0 in init.state
      init.state2["inc"] <- 0
      
      traj2 <- as.data.frame(ode(init.state2, Tafter, SIR_region2, theta, method = "ode45"))
      # compute incidence of each time interval
      traj2$inc <- c(0, diff(traj2$inc))
      
      #combine the Inc of traj1 and traj2
      for(i in 1:(T-1)) {
        traj1$Inc[i] <- traj1$inc[i] # Inc is the sum of inc of subregion 1 and 2
        traj1$inc2[i] <- 0 # inc2 is the inc of subregion 2
        traj1$inc1[i] <- traj1$inc[i] # inc1 is the inc of subregion 1
      }
      for(j in T:max(times)) {
        traj1$Inc[j] <- traj2$inc[j-(T-1)] + traj1$inc[j] # Inc is the sum of inc of subregion 1 and 2
        traj1$inc2[j] <- traj2$inc[j-(T-1)] # inc2 is the inc of subregion 2
        traj1$inc1[j] <- traj1$inc[j] # inc1 is the inc of subregion 1
      }
      
      return(traj1)
    }
  }
  
}

SIR_region_deter$rPointObs <- function(model.point, theta){
  
  ## the prevalence is observed through a Poisson process
  #obs.point <- rpois(n = 1, lambda = model.point[["Inc"]] * theta[["RR"]])
  #the prevalence is observed through a normal distribution process
  obs.point <- rnorm(n = 1, mean = model.point[["Inc"]] * theta[["RR"]], sd = 2)
  
  return(c(obs = obs.point))
}

SIR_region_deter$dPointObs <- function(data.point, model.point, theta, log = FALSE){
  
  ## the prevalence is observed through a Poisson process with a reporting rate
  #return(dpois(x = data.point[["obs"]], lambda = model.point[["Inc"]] * theta[["RR"]], log = log))
  #the prevalence is observed through a normal distribution process
  return(dnorm(x = data.point[["obs"]], mean = model.point[["Inc"]] * theta[["RR"]], sd = 0.1, log = log))
}

epi_time = c(1:100)

# simulation when time interval between onsets equals to nonzero.
#init.state <- c('S1' = 1000, 'I1' = 5, 'R1' = 1,'Inc1' = 0,'S2' = 2978, 'I2' = 2, 'R2' = 3,'Inc2' = 0)
lam = 2
init.state <- c('S2' = 27700/lam, 'I2' = 5, 'R2' = 1,'Inc2' = 0,'S1' = 27700, 'I1' = 2, 'R1' = 3,'Inc1' = 0)
#init.state <- c('S1' = 27700, 'I1' = 2, 'R1' = 3,'Inc1' = 0,'S2' = 27700/lam, 'I2' = 5, 'R2' = 1,'Inc2' = 0)
init.theta <- c(R0 = 3.2, D_inf = 3.6,D_lat =0, RR = 0.75, T=15) ##��T=0����ͼʱҪ�õ�310�е���315�еĴ���
#�ı�time interval T��basic reproductino number R0
#init.theta <- c(R0 = 1.9, D_inf = 3.6, RR = 0.75, T=30)
#simulation <- SIR_region_deter$simulate(theta = init.theta, init.state = init.state, times = epi_time)
library(plyr)
obs.traj <- rTrajObs(fitmodel = SIR_region_deter, theta = init.theta, init.state = init.state, times = epi_time)

sum(obs.traj$Inc)#total number of infected people of two regions
sum(obs.traj$Inc)/(27700+27700/lam)#attack rate
##################################################
### plot designed by Bo Xu
par(mar=c(5, 4, 2, 4))
plot(x = epi_time, y = obs.traj$Inc, type = "l", lty=1 ,col = "red", lwd=2, xlab = "Day", ylab = "Number of cases")#Inc is the sum of the incidence of subregion 1 and 2 calculated by model
lines(x=epi_time, y=obs.traj$inc1, type = "l",lty=2, col="dark gray", lwd=2 )#Inc1 is incidence of subregion 1 calculated by model
lines(x=c(init.theta[["T"]]:max(epi_time)), y=obs.traj$inc2[init.theta[["T"]]:max(epi_time)], type = "l", lty=3, col="dark gray", lwd=2 )#Inc2 is incidence of subregion 2 calculated by model
lines(x = epi_time, y = obs.traj$obs, type = "l", lty=1,col = "black",lwd=2)#obs is from Inc considering reporting rate and poisson observation process
#abline(v = init.theta[["T"]], type = "l", lty=2, col="blue",lwd=1)
legend("topright", legend = c("Incidence","Incidence of subregion 1","Incidence of subregion 2","Incidence observation"), bty = "n", lty = c(1,2,3,1), lwd = c(2,2,2,2),col = c("red","dark gray","dark gray","black") )
###############################################################

##################################################################################################################
##�ж�˫��
##########################����3
#GT = init.theta[["D_lat"]] + init.theta[["D_inf"]]#generation time
GT = init.theta[["D_inf"]]#generation time
final_size = sum(obs.traj$Inc)
i = which.max(obs.traj$Inc) #���ֵλ�ڴ������ҵĵڼ�������һ����time��
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
    if((Si>0)&&(Sj>0)){ # &&�Ǳ������߼������㣬&���������߼�������
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
  if( (i+2) <= length(obs.traj$Inc) ){#���������㣬�����TWr����
    #combine
    TW = rbind(TWl,TWr)
  }else{
    TW = TWl
  }
}else{
  TW = TWr
}

#ind1 = which.max(TW[,8]) #��TW[,8]�����ֵ��ֹһ��������ֵ�ظ�������which.max(TW[,8])ֻ�᷵��һ��λ��ֵ����which(TW[,8] == max(TW[,8]))��ȫ����������λ��ֵ
###����9�е��޸����ڶ�ȡTristan da cunha��ʵ�ʷ������ݺ�������
ind1 = which(TW[,8] == max(TW[,8])) #TW[,8]��P2; ind1����¼�Ķ���TW[,8]��P2�����ֵ���ڵ�λ��
if(length(ind1) == 1){ #���TW[,8]��P2�����ֵֻ������һ��λ�����Ҹ����ֵ����1,��ֱ��ѡȡ��
  if(TW[ind1,8] != 1) {ind = ind1}
  if(TW[ind1,8] == 1) {
    TW[ind1,11]
    #������Ϊ�㣬����֮���ͷ����
    TW[ind1,8] = 0
    ind66 = which(TW[,8] == max(TW[,8])) #TW[,8]��P2; ind66����¼�Ķ��Ǹ���֮��TW[,8]��P2�����ֵ���ڵ�λ��
    if(length(ind66) > 1){#���TW[,8]��P2�в�ֹһ����ȵ����ֵ
      ind77 = which(TW[ind66,11] == max(TW[ind66,11])) #ind66����¼�Ķ��Ǹ���֮��TW[,8]��P2�����ֵ���ڵ�λ��,��ind66����TW[,11]��W2�����ֵ��ind66�е�λ��
      if(length(ind77) == 1){ #����ind66�У�TW[,11]��W2�����ֵֻ��һ��λ�ó��֣���ѡȡ��
        ind00 = ind66[ind77]
      }
      if(length(ind77) > 1){ #����ind66�У�TW[,11]��W2�����ֵ�ڶ��λ�ó��֣�������Щλ�����ҵ�TW[,7]��|i-j|�����Ǹ�λ��
        ind88 = which.max(TW[ind66[ind77],7])
        ind00 = ind66[ind77][ind88]
      }
    }
    if(length(ind66) == 1){
      ind00 = ind66
    }
    #�Ƚ�ind1��ind00���Զ�Ӧ��TW[,11]��W2�Ĵ�С
    if(TW[ind1,11] >= TW[ind00,11]){
      ind = ind1
      TW[ind1,8] = 1 #��֮ǰ��Ϊ0��TW[ind,8]���û�ԭֵ1
    }
    if(TW[ind1,11] < TW[ind00,11]){
      ind = ind00
    }
    #��TW[,8]��P2�Ķ����ȵ����ֵ����1��ǰ����,��֮ǰ��Ϊ0��TW[ind1,8]���û�ԭֵ1
    TW[ind1,8] = 1
  }
}

if(length(ind1) > 1){ #���TW[,8]��P2�����ֵ�����ڲ�ֹһ��λ����
  #####################################################
  #��TW[,8]��P2�Ķ����ȵ����ֵ����1
  if(TW[ind1,8][1] == 1){ #P2����1����vij����0
    ind2 = which(TW[ind1,11] == max(TW[ind1,11])) #ind1����¼�Ķ���TW[,8]��P2�����ֵ���ڵ�λ��,��ind1����TW[,11]��W2�����ֵ��ind1�е�λ��
    if(length(ind2) == 1){ #����ind1�У�TW[,11]��W2�����ֵֻ��һ��λ�ó��֣���ѡȡ��
      ind = ind1[ind2]
    }
    if(length(ind2) > 1){ #����ind1�У�TW[,11]��W2�����ֵ�ڶ��λ�ó��֣�������Щλ�����ҵ�TW[,7]��|i-j|�����Ǹ�λ��
      ind3 = which.max(TW[ind1[ind2],7])
      ind = ind1[ind2][ind3]
    }
    
    #����ȫ����Ϊ�㣬����֮���ͷ����
    TW[ind1,8] = 0
    ind11 = which(TW[,8] == max(TW[,8])) #TW[,8]��P2; ind11����¼�Ķ��Ǹ���֮��TW[,8]��P2�����ֵ���ڵ�λ��
    if(length(ind11) > 1){#���TW[,8]��P2�в�ֹһ����ȵ����ֵ
      ind22 = which(TW[ind11,11] == max(TW[ind11,11])) #ind11����¼�Ķ��Ǹ���֮��TW[,8]��P2�����ֵ���ڵ�λ��,��ind11����TW[,11]��W2�����ֵ��ind11�е�λ��
      if(length(ind22) == 1){ #����ind11�У�TW[,11]��W2�����ֵֻ��һ��λ�ó��֣���ѡȡ��
        ind0 = ind11[ind22]
      }
      if(length(ind22) > 1){ #����ind11�У�TW[,11]��W2�����ֵ�ڶ��λ�ó��֣�������Щλ�����ҵ�TW[,7]��|i-j|�����Ǹ�λ��
        ind33 = which.max(TW[ind11[ind22],7])
        ind0 = ind11[ind22][ind33]
      }
    }
    if(length(ind11) == 1){
      ind0 = ind11
    }
    
    #�Ƚ�ind��ind0���Զ�Ӧ��TW[,11]��W2�Ĵ�С
    if(TW[ind,11] >= TW[ind0,11]){
      ind = ind
      TW[ind,8] = 1 #��֮ǰ��Ϊ0��TW[ind,8]���û�ԭֵ1
    }
    if(TW[ind,11] < TW[ind0,11]){
      ind = ind0
    }
    #��TW[,8]��P2�Ķ����ȵ����ֵ����1��ǰ����,��֮ǰ��Ϊ0��TW[ind1,8]���û�ԭֵ1
    TW[ind1,8] = 1
  }
  
  #####################################################
  #��TW[,8]��P2�Ķ����ȵ����ֵ������1
  if(TW[ind1,8][1] != 1){
    ind2 = which(TW[ind1,11] == max(TW[ind1,11])) #ind1����¼�Ķ���TW[,8]��P2�����ֵ���ڵ�λ��,��ind1����TW[,11]��W2�����ֵ��ind1�е�λ��
    if(length(ind2) == 1){ #����ind1�У�TW[,11]��W2�����ֵֻ��һ��λ�ó��֣���ѡȡ��
      ind = ind1[ind2]
    }
    if(length(ind2) > 1){ #����ind1�У�TW[,11]��W2�����ֵ�ڶ��λ�ó��֣�������Щλ�����ҵ�TW[,7]��|i-j|�����Ǹ�λ��
      ind3 = which.max(TW[ind1[ind2],7])
      ind = ind1[ind2][ind3]
    }
  }
  #ind2 = which(TW[ind1,11] == max(TW[ind1,11])) #TW[,11]��W2
  ##ind2 = which.max(TW[ind1,11]) #��TW[ind1,11]�����ֵ��ֹһ������ֻȡ����һ���������������
  ##ind = ind1[ind2]
}

ind
P2 = TW$P2[ind] #(Sj-vij)/Sj, Sj is the value of the second peak
W2 = TW$W2[ind] #(Sj-vij)/(Si-vij)

thre1 = 0.5
thre2 = 10
thre3 = 3 #�ڵ�����16��ģ�͵�17�����ߣ�roadNetworkģ����country��city�������ߣ�֮��ȷ����
#TW$TW[ind] #2-peak metric
#TW$`|j-i|`[ind] #time gap between two peaks
if(P2 >= thre1) {
  if((TW$`|j-i|`[ind] >= thre3*GT) && (W2 >= (1/thre2))){peakNum = 2}
  else {
    if((TW$`|j-i|`[ind] < thre3*GT) && (W2 >= (1/thre2))) {peakNum = 1.8} #time gap between two peaks is not long enough to generate a second generation
    else {peakNum = 1} #ֻҪW2 < (1/thre2)����ô����Ϊ�ǵ���
  }
}

if((P2 < thre1)&(P2>0)){
  if((TW$`|j-i|`[ind] >= thre3*GT) && (W2 >= (1/thre2))) {peakNum = 1.5} #the valley is not deep enough to be looked as an interwave period, but the time gap is enough
  else {
    if((TW$`|j-i|`[ind] < thre3*GT) && (W2 >= (1/thre2))) {peakNum = 1.2}
    else {peakNum = 1} #ֻҪW2 < (1/thre2)����ô����Ϊ�ǵ���
  }
}
if(P2==0) {peakNum = 1}#ֻҪP2==0����ôҲ��Ϊ�ǵ���
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
        else {pN = 1} #ֻҪW2 < (1/thre2)����ô����Ϊ�ǵ���
      }
    }
    if((P2 < thre1)&(P2>0)){
      if((TW$`|j-i|`[k] >= thre3*GT) && (W2 >= (1/thre2))) {pN = 1.5} #the valley is not deep enough to be looked as an interwave period, but the time gap is enough
      else {
        if((TW$`|j-i|`[k] < thre3*GT) && (W2 >= (1/thre2))) {pN = 1.2}
        else {pN = 1} #ֻҪW2 < (1/thre2)����ô����Ϊ�ǵ���
      }
    }
    if(P2==0) {pN = 1}#ֻҪP2==0����ôҲ��Ϊ�ǵ���
    
    if(pN > peakNum) {
      tk = cbind(tk,k)
      tpN = cbind(tpN,pN)
    }
  }
  
  if((length(tk) > 1) && (length(tpN) > 1)){
    #tk = tk[,-1]
    #tpN = tpN[,-1]
    peakNum = max(tpN)
    ind44=which(tpN == max(tpN)) #��tpN���ֵ���ֵ�����λ�ö�����ind44��
    if(length(ind44)==1) { #��tpN���ֵֻ��һ��λ�ó���
      ind = tk[ind44]
    }
    if(length(ind44) > 1) { #��tpN���ֵ�ڶ��λ�ó���
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

###############################################################
###############################################################
# simulation when time interval between onsets equals to zero.
init.state <- c('S1' = 2978, 'I1' = 2, 'R1' = 3,'Inc1' = 0,'S2' = 1000, 'I2' = 60, 'R2' = 1,'Inc2' = 0)
#init.state <- c('S2' = 1000, 'I2' = 60, 'R2' = 1,'Inc2' = 0,'S1' = 2978, 'I1' = 2, 'R1' = 3,'Inc1' = 0)
init.theta <- c(R0 = 1.9, D_inf = 3.6, RR = 0.75, T=0)
#simulation <- SIR_region_deter$simulate(theta = init.theta, init.state = init.state, times = epi_time)
obs.traj <- rTrajObs(fitmodel = SIR_region_deter, theta = init.theta, init.state = init.state, times = epi_time)
plot(x = epi_time, y = obs.traj$Inc, type = "l", lty=2 ,col = "black", lwd=3, xlab = "time(days)", ylab = "Incidence from model")#Inc is the sum of the incidence of subregion 1 and 2 calculated by model
lines(x=epi_time, y=obs.traj$Inc1, type = "l",lty=1, col="blue", lwd=1 )#Inc1 is incidence of subregion 1 calculated by model
lines(x=epi_time, y=obs.traj$Inc2, type = "l", lty=1, col="green", lwd=1 )#Inc2 is incidence of subregion 2 calculated by model
lines(x = epi_time, y = obs.traj$obs, type = "l", lty=1,col = "red",lwd=1)#obs is from Inc considering reporting rate and poisson observation process
#abline(v = init.theta[["T"]], type = "l", lty=2, col="gray",lwd=0.5)
legend("topright", legend = c("incidence of the whole region","incidence of subregion 1","incidence of subregion 2","observation of incidence"),lty = c(2,1,1,1), lwd = c(3,1,1,1),col = c("black","blue","green","red") )
