#########################################################################################################################################
#########################################################################################################################################
##beta1(amplitude of fluctuating part of beta)和beta0(average value of beta)
step1 = 0.05 #步长
step2 = 0.1
begin1 = 0 #begin 必须是零，否则会出错！！！
begin2 = 0.1
pbeta1 = seq(begin1, 1, by=step1) #pbeta1是为beta1设定的取值范围
pbeta0 = seq(begin2, 1.8, by=step2) #pbeta0是为beta0设定的取值范围
PNmat = matrix(99, nrow = length(pbeta1), ncol = length(pbeta0)) #存储number of waves
ATmat = matrix(0, nrow = length(pbeta1), ncol = length(pbeta0)) #存储attack rate
epi_time = c(151:350)
N=27700
library(plyr)

ptm = proc.time() #计时开始
for(b1 in pbeta1){
  for(b0 in pbeta0){
    
    init.state <- c('S' = 0.794*N, 'E'=0*N, 'I' = 0.006*N, 'R' = 0.2*N,'Inc' = 0*N)
    init.theta <- c(beta0 = b0, beta1 = b1, D_lat = 1, D_inf = 3, RR = 0.75)
    #obs.traj <- PTR$simulate(init.theta = init.theta, init.state = init.state, times = epi_time)
    obs.traj <- rTrajObs(fitmodel = PTR, theta = init.theta, init.state = init.state, times = epi_time)
    
    Nreal = init.state[["S"]]
    #attack rate
    AT = sum(obs.traj$Inc)/Nreal
    #ATmat[-4+(r01/step1),-4+(r02/step2)] = AT
    A1 = ((b1-begin1)/step1)+1
    B1 = ((b0-begin2)/step2)+1
    #print(ATmat[A1,B1])
    ATmat[A1,B1] <- AT
    if(b0 < b1) {
      ATmat[A1,B1] <- -1
    }
    ####################################################################################################################
    ##判断two wave
    ##########################方案2
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
    #print("ind=")
    #ind
    #print("peakNum=")
    #peakNum
    
    PNmat[((b1-begin1)/step1)+1,((b0-begin2)/step2)+1] <- peakNum
    if(b0 < b1) {
      PNmat[((b1-begin1)/step1)+1,((b0-begin2)/step2)+1] <- -1
    }
  }
}

proc.time() - ptm #计时结束

#plot method 1
source("http://www.phaget4.org/R/myImagePlot.R")
#setwd("E:/科研/paper/bi-modal/figures")
#myImagePlot(PNmat, xLabels=pbeta0,yLabels=pbeta1, title=c("Number of waves"), zlim=c(1,2)) 
myImagePlot(PNmat[,8:ncol(ATmat)], xLabels=pbeta0[8:length(pbeta0)],yLabels=pbeta1, title=c("Number of waves"), zlim=c(1,2))
#myImagePlot(t(PNmat), xLabels=pC11,yLabels=pA22, title=c("Number of waves"), zlim=c(1,2)) 
#myImagePlot(ATmat, xLabels=pbeta0,yLabels=pbeta1, title=c("Attack rate"), zlim=c(0,max(ATmat))) 
myImagePlot(ATmat[,8:ncol(ATmat)], xLabels=pbeta0[8:length(pbeta0)],yLabels=pbeta1, title=c("Attack rate"), zlim=c(min(abs(ATmat[,8:ncol(ATmat)])),max(ATmat))) 

length(which(PNmat==2))
length(which(PNmat==1.8))
length(which(PNmat==1.5))
length(which(PNmat==1.2))
length(which(PNmat==1))

#确定双峰对应的参数组合
ind2 = which(PNmat==2) 
Num_two_param = length(ind2) #能得到双峰的参数组合的数目
two_param = matrix(data = 0, nrow = Num_two_param, ncol = 2) #存放 能得到双峰的参数组合 的矩阵
colnames(two_param)=c('beta1','beta0')
for(i in c(1:Num_two_param)) {
  if((ind2[i] %% nrow(PNmat)) == 0) {
    two_param[i,1] = max(pbeta1)
    two_param[i,2] = begin2 + step2 * (ind2[i] %/% nrow(PNmat) - 1)
  }
  else {
    two_param[i,1] = begin1 + step1 * ((ind2[i] %% nrow(PNmat)) - 1) # "%%"表示取余。计算得到beta1
    two_param[i,2] = begin2 + step2 * (ind2[i] %/% nrow(PNmat)) # "%/%"表示整除。计算得到beta0
  }
}

#(beta1,beta0)的取值对(0.35,0.3),(0.7,0.6),(0.8,0.7)是不符合"beta1小于等于beta0"的，要去除
two_param = two_param[-c(1,10,13),]
cor(two_param[,1], two_param[,2], method = "pearson")
cor.test(two_param[,1], two_param[,2], method = "pearson")
lm.sol = lm(two_param[,2] ~ two_param[,1])
summary(lm.sol) #Multiple R-squared:0.8339, p-value:< 2.2e-16

lowAttackRate_Param = read.csv("E:/科研/paper/bi-modal/figures/5Periodic-Transmission-Rate/lowAttackRate_Param_b1&b0.csv")
cor.test(lowAttackRate_Param[,1],lowAttackRate_Param[,2],method = "pearson")
lm.sol = lm(lowAttackRate_Param$beta0 ~ lowAttackRate_Param$beta1)
#lm.sol = lm(lowAttackRate_Param$beta1 ~ lowAttackRate_Param$beta0)
summary(lm.sol) #Multiple R-squared:0.994, p-value:1.06e-14


#保存到本地
write.csv(PNmat,file = "E:/科研/paper/bi-modal/figures/5Periodic-Transmission-Rate/PNmat_b1&b0.csv") #存储number of waves
write.csv(ATmat, file = "E:/科研/paper/bi-modal/figures/5Periodic-Transmission-Rate/ATmat_b1&b0.csv") #存储attack rate
write.csv(two_param, file = "E:/科研/paper/bi-modal/figures/5Periodic-Transmission-Rate/twoWaveParam_b1&b0.csv") #存储 能得到双峰的参数组合
