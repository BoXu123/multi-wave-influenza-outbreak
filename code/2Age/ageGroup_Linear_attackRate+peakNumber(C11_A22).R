#AT=0 #Attack Rate
#PN = 0 #number of peak
#########################################################################################################################################
#########################################################################################################################################
##C11和A22 CA12固定为5
step = 5 #步长
pC11 = seq(0, 300, by=step)#C11表示孩子们一年中有效接触次数,C对应children，pC11表示设定的C11的取值范围
pA22 = seq(0, 150, by=step)#A22表示成人一年中的有效接触次数，A对应Adults,pA22表示设定的A22的取值范围
PNmat = matrix(0, nrow = length(pC11), ncol = length(pA22)) #存储number of waves
ATmat = matrix(0, nrow = length(pC11), ncol = length(pA22)) #存储attack rate
CA12 = 5
#A22 = 53
epi_time = c(1:200)
N = 27700

ptm = proc.time() #计时开始
for(C11 in pC11){
  for(A22 in pA22){

    init.state <- c('N'=N,'S1a' = 0.3, 'E1a'=0, 'I1a' = 0.0012, 'R1a' = 0,'S2a' = 0.22, 'E2a'=0, 'I2a' = 0.001, 'R2a' = 0.02,'Inc1' = 0,'Inc2' = 0,'Inc' = 0)
    init.theta <- c(beta11 = C11/365, beta12 = CA12/365, beta21 = CA12/365, beta22 = A22/365, l = 0, D_lat = 365/150, D_inf = 365/10, RR=0.75, T=27, Tint=18, 
                    a1=-0.0023, b1=0.3023, a2=-0.01439, b2=0.665)
    #Nreal = N*(init.state[["S1a"]] +init.state[["E1a"]] + init.state[["I1a"]] +init.state[["R1a"]] + init.state[["S2a"]] +init.state[["E2a"]] + init.state[["I2a"]] +init.state[["R2a"]])
    Nreal = N*(init.state[["S1a"]] + init.state[["S2a"]])
    obs.traj <- rTrajObs(fitmodel = age.Linear, theta = init.theta, init.state = init.state, times = epi_time)
    #attack rate
    #AT = cbind(AT,sum(obs.traj$Inc)/Nreal)
    AT = sum(obs.traj$Inc)/Nreal
    ATmat[1+C11/step,1+A22/step] = AT
    
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
    
    
    #Number of peak
    #PN = cbind(PN,peakNum)
    PNmat[1+C11/step,1+A22/step] = peakNum
  }
  
}
proc.time() - ptm #计时结束

#PN = PN[,-1]
#plot(x=pC11, y=PN, type = "p", lty=1 ,col = "dark green", lwd=2, xlab = "Transmission rate among children", ylab = "Attack rate")

#plot method 1
source("http://www.phaget4.org/R/myImagePlot.R")
#setwd("E:/科研/paper/bi-modal/figures")
myImagePlot(PNmat, xLabels=pA22,yLabels=pC11, title=c("Number of waves"), zlim=c(1,2)) 
#myImagePlot(t(PNmat), xLabels=pC11,yLabels=pA22, title=c("Number of waves"), zlim=c(1,2)) 
myImagePlot(ATmat, xLabels=pA22,yLabels=pC11, title=c("Attack rate"), zlim=c(min(ATmat),max(ATmat))) 

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
two_param = matrix(data = 0, nrow = Num_two_param, ncol = 2) #存放 能得到双峰的参数组合
colnames(two_param)=c('C11_Child','A22_Adult')
for(i in c(1:Num_two_param)) {
  if((ind2[i] %% nrow(PNmat)) == 0) {
    two_param[i,1] = max(pC11)
    two_param[i,2] = step * (ind2[i] %/% nrow(PNmat) - 1)
  }
  else {
    two_param[i,1] = step * ((ind2[i] %% nrow(PNmat)) - 1) # "%%"表示取余。计算得到C11
    two_param[i,2] = step * (ind2[i] %/% nrow(PNmat)) # "%/%"表示整除。计算得到CA12
  }
}

#保存到本地
write.csv(PNmat,file = "E:/科研/paper/bi-modal/figures/7Two-Age-Group/PNmat_C11&A22.csv") #存储number of waves
write.csv(ATmat, file = "E:/科研/paper/bi-modal/figures/7Two-Age-Group/ATmat_C11&A22.csv") #存储attack rate
write.csv(two_param, file = "E:/科研/paper/bi-modal/figures/7Two-Age-Group/twoWaveParam_C11&A22.csv") #存储 能得到双峰的参数组合
#参数的相关系数
test1=read.csv("E:/科研/paper/bi-modal/figures/7Two-Age-Group/twoWaveParam_C11&A22.csv")
cor.test(test1$C11_Child, test1$A22_Adult,method = "pearson")


######################################################################################################
######################################################################################################
######plot epidemic curve (Begin)
CA12 = 5

C11 = 190 #固定C11为135，取A22为15、55、80、105
A22 = 15
epi_time = c(1:200)
N = 27700
init.state <- c('N'=N,'S1a' = 0.3, 'E1a'=0, 'I1a' = 0.0012, 'R1a' = 0,'S2a' = 0.22, 'E2a'=0, 'I2a' = 0.001, 'R2a' = 0.02,'Inc1' = 0,'Inc2' = 0,'Inc' = 0)
init.theta <- c(beta11 = C11/365, beta12 = CA12/365, beta21 = CA12/365, beta22 = A22/365, l = 0, D_lat = 365/150, D_inf = 365/10, RR=0.75, T=27, Tint=18, 
                a1=-0.0023, b1=0.3023, a2=-0.01439, b2=0.665)
obs.traj <- rTrajObs(fitmodel = age.Linear, theta = init.theta, init.state = init.state, times = epi_time)

r=0
g=0
b=255
par(mar=c(5, 4, 2, 4))
lines(x = epi_time[-c(length(epi_time), length(epi_time)-1)], y = obs.traj$Inc, type = "l", lty=1 ,col = rgb(r/255,g/255,b/255), lwd=2)

plot(x = epi_time[-c(length(epi_time), length(epi_time)-1)], y = obs.traj$Inc, type = "l", lty=1 ,col = rgb(r/255,g/255,b/255), lwd=2, xlab = "Day", ylab = "Number of cases")
legend("topleft", legend = c("C11=190","CA12 = 5"),bty = "n")
legend("topright", legend = c("A22=15","A22=55","A22=80","A22=105"),bty = "n", lty = c(1,1,1), lwd = c(2,2,2),col = c("#0000FF","#FFFF00","#7F7F80","#3333CC") )
legend(x=130, y=105, legend = c("PN=1","PN=2","PN=1.5","PN=1.2"),bty = "n")
######plot epidemic curve (End)


#########################################################################################################################################
#########################################################################################################################################
#plot method 2
install.packages("corrplot")
require(corrplot)
corrplot(1/PNmat, method = "color")
corrplot(ATmat, method = "color")

#plot method 3
library(reshape2)
library(ggplot2)
#PNmat
longData<-melt(PNmat)
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Transmission rate among adults", y="Transmission rate among children", title="Number of waves") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
#ATmat
longData<-melt(PNmat)
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Transmission rate among adults", y="Transmission rate among children", title="Number of waves") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
