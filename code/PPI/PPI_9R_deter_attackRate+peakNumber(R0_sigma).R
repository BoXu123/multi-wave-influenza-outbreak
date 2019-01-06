#R0和alpha
step1 = 0.25 #步长
step2 = 0.025
begin1 = 0 #begin 必须是零，否则会出错！！！
begin2 = 0
pR0 = seq(begin1, 15, by=step1)
psigma = seq(begin2, 1, by=step2)
PNmat = matrix(99, nrow = length(pR0), ncol = length(psigma)) #存储number of waves
ATmat = matrix(99, nrow = length(pR0), ncol = length(psigma)) #存储attack rate
epi_time = c(1:100)

ptm = proc.time() #计时开始
for(r0 in pR0){
  for(sig in psigma){
    init.state <- c('S' = 27700, 'E' = 0, 'I' = 1, 'R1' = 0, 'R2' = 0,'R3' = 0, 'R4' = 0,'R5' = 0, 'R6' = 0,'R7' = 0, 'R8' = 0,'R9' = 0, 'L'= 6,'RL' = 0,'Inc' = 0)
    init.theta <- c(R0 = r0, D_lat = 1.8, D_inf = 1.0, D_imm = 5.7, sigma = sig, rho = 0.68)
    obs.traj <- rTrajObs(fitmodel = PPI_9R_deter, theta = init.theta, init.state = init.state, times = epi_time)
    
    Nreal = init.state[["S"]] + sum(obs.traj$RL)
    AT = sum(obs.traj$Inc)/Nreal
    A1 = ((r0-begin1)/step1)+1
    B1 = ((sig-begin2)/step2)+1
    #print(ATmat[A1,B1])
    ATmat[A1,B1] <- AT
    
    ##################################################################################################################
    ##判定双峰
    ##########################方案2
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
    thre2 = 6 #此处设为10的话，出的图与设6是一样的 #W2的最大值为0.194970448，thre2要大于1/0.194970448=5.13才能使peakNum为2
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
    
    #Number of peak
    #PN = cbind(PN,peakNum)
    #PNmat[-4+(r01/step1),-4+(r02/step2)] = peakNum
    PNmat[((r0-begin1)/step1)+1,((sig-begin2)/step2)+1] <- peakNum
    
  }
}

proc.time() - ptm #计时结束

#plot method 1
source("http://www.phaget4.org/R/myImagePlot.R")
#setwd("E:/科研/paper/bi-modal/figures")
myImagePlot(PNmat, xLabels=psigma, yLabels=pR0, title=c("Number of waves"), zlim=c(1,2)) 
#myImagePlot(t(PNmat), xLabels=pC11,yLabels=pA22, title=c("Number of waves"), zlim=c(1,2)) 
myImagePlot(ATmat, xLabels=psigma, yLabels=pR0, title=c("Attack rate"), zlim=c(min(ATmat),max(ATmat))) 

length(which(PNmat==2))
length(which(PNmat==1.8))
length(which(PNmat==1.5))
length(which(PNmat==1.2))
length(which(PNmat==1))

#确定双峰对应的参数组合
ind2 = which(PNmat==2) 
#ind2 = which(PNmat==1.8) #这里没有峰数为2的，只有峰数为1.8的
Num_two_param = length(ind2) #能得到双峰的参数组合的数目
two_param = matrix(data = 0, nrow = Num_two_param, ncol = 2) #存放 能得到双峰的参数组合 的矩阵
colnames(two_param)=c('R0','sigma')
for(i in c(1:Num_two_param)) {
  if((ind2[i] %% nrow(PNmat)) == 0) {
    two_param[i,1] = max(pR0)
    two_param[i,2] = begin2 + step2 * (ind2[i] %/% nrow(PNmat) - 1)
  }
  else {
    two_param[i,1] = begin1 + step1 * ((ind2[i] %% nrow(PNmat)) - 1) # "%%"表示取余。计算得到R0
    two_param[i,2] = begin2 + step2 * (ind2[i] %/% nrow(PNmat)) # "%/%"表示整除。计算得到sigma
  }
}

#保存到本地
write.csv(PNmat,file = "E:/科研/paper/bi-modal/figures/13Partially Protective Immunity (PPI)/PNmat_R0&sigma.csv") #存储number of waves
write.csv(ATmat, file = "E:/科研/paper/bi-modal/figures/13Partially Protective Immunity (PPI)/ATmat_R0&sigma.csv") #存储attack rate
write.csv(two_param, file = "E:/科研/paper/bi-modal/figures/13Partially Protective Immunity (PPI)/twoWaveParam_R0&sigma.csv") #存储 能得到双峰的参数组合

#########################################################################################################################################
#########################################################################################################################################
######plot epidemic curves (Begin)
sig = 0.35 #固定sig为0.35，R0取值为3.5、5、7.5
r0 = 3.5

epi_time = c(1:100)
init.state <- c('S' = 27700, 'E' = 0, 'I' = 1, 'R1' = 0, 'R2' = 0,'R3' = 0, 'R4' = 0,'R5' = 0, 'R6' = 0,'R7' = 0, 'R8' = 0,'R9' = 0, 'L'= 6,'RL' = 0,'Inc' = 0)
init.theta <- c(R0 = r0, D_lat = 1.8, D_inf = 1.0, D_imm = 5.7, sigma = sig, rho = 0.68)
obs.traj <- rTrajObs(fitmodel = PPI_9R_deter, theta = init.theta, init.state = init.state, times = epi_time)

r=0
g=0
b=255
par(mar=c(5, 4, 2, 4))
lines(x = epi_time, y = obs.traj$Inc, type = "l", lty=1 ,col = rgb(r/255,g/255,b/255), lwd=2)

plot(x = epi_time, y = obs.traj$Inc, type = "l", lty=1 ,col = rgb(r/255,g/255,b/255), lwd=2, xlab = "Day", ylab = "Number of cases")
legend("topleft", legend = c("σ = 0.35"),bty = "n")
legend("topright", legend = c("R0 = 3.5","R0 = 5","R0 = 7.5"),bty = "n", lty = c(1,1,1), lwd = c(2,2,2),col = c("#0000FF","#7F7F80","#FFFF00") )
legend(x=60, y=4000, legend = c("PN=1","PN=1.5","PN=2"),bty = "n")
##PN为1的对应RGB颜色为(0,0,255),PN为1.2的对应RGB颜色为(51,51,204),PN为1.5的对应RGB颜色为(127,127,128)，PN为1.8的对应RGB颜色为(204,204,51),PN为2的对应RGB颜色为(255,255,0)
##PN为1的对应颜色为"#0000FF"，PN为1.2的对应颜色为"#3333CC"，PN为1.5的对应颜色为"#7F7F80"，PN为1.8的对应颜色为"#CCCC33"，PN为2的对应颜色为"#FFFF00"
######plot epidemic curves (End)