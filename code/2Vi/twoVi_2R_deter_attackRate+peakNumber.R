#########################################################################################################################################
#########################################################################################################################################
##R01和R02
step1 = 0.2 #步长
step2 = 0.2
begin1 = 0 #begin 必须是零，否则会出错！！！
begin2 = 0
pR01 = seq(begin1, 5, by=step1) #pR01是为R01设定的取值范围
pR02 = seq(begin2, 5, by=step2) #pR02是为R02设定的取值范围
PNmat = matrix(99, nrow = length(pR01), ncol = length(pR02)) #存储number of waves
ATmat = matrix(0, nrow = length(pR01), ncol = length(pR02)) #存储attack rate
epi_time = c(1:100)

ptm = proc.time() #计时开始
for(r01 in pR01){
  for(r02 in pR02){
    init.state <- c('S' = 27700, 'E1' = 0, 'I1' = 20, 'T11' = 0, 'T12' = 0, 'L1'= 3, 'E2'=0, 'I2'=10, 'T21' = 0, 'T22' = 0, 'L12' = 0, 'Inc1' = 0, 'Inc2' = 0, 'Inc' = 0)
    init.theta <- c(R01 = r01, R02=r02, D_lat = 1.5, D_inf = 2, D_imm = 0.3, rho = 0.65)
    obs.traj <- rTrajObs(fitmodel = twoVi_2R_deter, theta = init.theta, init.state = init.state, times = epi_time)
    #Nreal = init.state[["S"]]+init.state[["E1"]]+init.state[["I1"]]+init.state[["T11"]]+init.state[["T12"]]+init.state[["L1"]]+init.state[["E2"]]+init.state[["I2"]]+init.state[["T21"]]+init.state[["T22"]]+init.state[["L12"]]
    Nreal = init.state[["S"]]+sum(obs.traj$IncL1)
    #attack rate
    #AT = cbind(AT,sum(obs.traj$Inc)/Nreal)
    AT = sum(obs.traj$Inc)/Nreal
    #ATmat[-4+(r01/step1),-4+(r02/step2)] = AT
    A1 = ((r01-begin1)/step1)+1
    B1 = ((r02-begin2)/step2)+1
    #print(ATmat[A1,B1])
    ATmat[A1,B1] <- AT
    #print(ATmat[A1,B1])
    
    ##################################################################################################################
    ##判定双峰
    ##########################方案3
    #GT = init.theta[["D_lat"]] + init.theta[["D_inf"]]#generation time
    #GT = (init.theta[["D_lat"]] + init.theta[["D_inf"]] + init.theta[["D_imm"]]) * 2 #generation time
    GT = init.theta[["D_lat"]] + 0.5*init.theta[["D_inf"]] #generation time
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
    
    #Number of peak
    #PN = cbind(PN,peakNum)
    #PNmat[-4+(r01/step1),-4+(r02/step2)] = peakNum
    PNmat[((r01-begin1)/step1)+1,((r02-begin2)/step2)+1] <- peakNum
  }
  
}
proc.time() - ptm #计时结束

#PN = PN[,-1]
#plot(x=pC11, y=PN, type = "p", lty=1 ,col = "dark green", lwd=2, xlab = "Transmission rate among children", ylab = "Attack rate")

#plot method 1
source("http://www.phaget4.org/R/myImagePlot.R")
#setwd("E:/科研/paper/bi-modal/figures")
myImagePlot(PNmat, xLabels=pR01,yLabels=pR02, title=c("Number of waves"), zlim=c(1,2)) 
#myImagePlot(t(PNmat), xLabels=pC11,yLabels=pA22, title=c("Number of waves"), zlim=c(1,2)) 
myImagePlot(ATmat, xLabels=pR01,yLabels=pR02, title=c("Attack rate"), zlim=c(min(ATmat),max(ATmat))) 

length(which(PNmat==2))
length(which(PNmat==1.8))
length(which(PNmat==1.5))
length(which(PNmat==1.2))
length(which(PNmat==1))
sumpix = length(which(PNmat==1))+length(which(PNmat==1.2))+length(which(PNmat==1.5))+length(which(PNmat==1.8))
par(mar=c(5, 4, 2, 4))
plot(x=c(1,1.2,1.5,1.8,2), y=c(length(which(PNmat==1))/sumpix*100,length(which(PNmat==1.2))/sumpix*100,length(which(PNmat==1.5))/sumpix*100,length(which(PNmat==1.8))/sumpix*100,length(which(PNmat==2))/sumpix*100),
     type = "h",lty=1 ,col = "red", lwd=2, xlab = "peakNum", ylab = "Frequency")

#确定双峰对应的参数组合
ind2 = which(PNmat==2) 
Num_two_param = length(ind2) #能得到双峰的参数组合的数目
two_param = matrix(data = 0, nrow = Num_two_param, ncol = 2) #存放 能得到双峰的参数组合 的矩阵
colnames(two_param)=c('R01','R02')
for(i in c(1:Num_two_param)) {
  if((ind2[i] %% nrow(PNmat)) == 0) {
    two_param[i,1] = max(pR01)
    two_param[i,2] = begin2 + step2 * (ind2[i] %/% nrow(PNmat) - 1)
    }
  else {
    two_param[i,1] = begin1 + step1 * ((ind2[i] %% nrow(PNmat)) - 1) # "%%"表示取余。计算得到R01
    two_param[i,2] = begin2 + step2 * (ind2[i] %/% nrow(PNmat)) # "%/%"表示整除。计算得到R02
  }
}

#保存到本地
write.csv(PNmat,file = "E:/科研/paper/bi-modal/figures/8Two-viruses/PNmat_R01&R02.csv") #存储number of waves
write.csv(ATmat, file = "E:/科研/paper/bi-modal/figures/8Two-viruses/ATmat_R01&R02.csv") #存储attack rate
write.csv(two_param, file = "E:/科研/paper/bi-modal/figures/8Two-viruses/twoWaveParam_R01&R02.csv") #存储 能得到双峰的参数组合

#################################################################################
#################################################################################
####################plot epidemic curve (Begin)
epi_time = c(1:100)
r01 = 2 #(r01,r02)取(4,4)、(2,2)
r02 = 2

init.state <- c('S' = 27700, 'E1' = 0, 'I1' = 20, 'T11' = 0, 'T12' = 0, 'L1'= 3, 'E2'=0, 'I2'=10, 'T21' = 0, 'T22' = 0, 'L12' = 0, 'Inc1' = 0, 'Inc2' = 0, 'Inc' = 0)
init.theta <- c(R01 = r01, R02=r02, D_lat = 1.5, D_inf = 2, D_imm = 0.3, rho = 0.65)
obs.traj <- rTrajObs(fitmodel = twoVi_2R_deter, theta = init.theta, init.state = init.state, times = epi_time)

r=0
g=0
b=255
par(mar=c(5, 4, 2, 4))
lines(x = epi_time, y = obs.traj$Inc, type = "l", lty=2 ,col = rgb(r/255,g/255,b/255), lwd=2)

plot(x = epi_time, y = obs.traj$Inc, type = "l", lty=1 ,col = rgb(r/255,g/255,b/255), lwd=2, xlab = "Day", ylab = "Number of cases")
#legend("topleft", legend = c("R01=1.6"),bty = "n")
legend("topright", legend = c("R01 = 4, R02 = 4","R01 = 2, R02 = 2"),bty = "n", lty = c(1,1), lwd = c(2,2),col = c("#FFFF00","#0000FF") )
legend(x=65, y=3000, legend = c("PN=2","PN=1"),bty = "n")

#"#0000FF","#FFFF00","#7F7F80","#FFFF00"
####################plot epidemic curve (End)