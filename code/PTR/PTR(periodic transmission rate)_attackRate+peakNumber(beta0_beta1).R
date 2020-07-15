#########################################################################################################################################
#########################################################################################################################################
##beta1(amplitude of fluctuating part of beta)��beta0(average value of beta)
step1 = 0.05 #����
step2 = 0.1
begin1 = 0 #begin �������㣬��������������
begin2 = 0.1
pbeta1 = seq(begin1, 1, by=step1) #pbeta1��Ϊbeta1�趨��ȡֵ��Χ
pbeta0 = seq(begin2, 1.8, by=step2) #pbeta0��Ϊbeta0�趨��ȡֵ��Χ
PNmat = matrix(99, nrow = length(pbeta1), ncol = length(pbeta0)) #�洢number of waves
ATmat = matrix(0, nrow = length(pbeta1), ncol = length(pbeta0)) #�洢attack rate
epi_time = c(151:350)
N=27700
library(plyr)

ptm = proc.time() #��ʱ��ʼ
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
    ##�ж�two wave
    ##########################����2
    GT = init.theta[["D_lat"]] + 0.5*init.theta[["D_inf"]]#generation time
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

proc.time() - ptm #��ʱ����

#plot method 1
source("http://www.phaget4.org/R/myImagePlot.R")
#setwd("E:/����/paper/bi-modal/figures")
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

#ȷ��˫���Ӧ�Ĳ������
ind2 = which(PNmat==2) 
Num_two_param = length(ind2) #�ܵõ�˫��Ĳ�����ϵ���Ŀ
two_param = matrix(data = 0, nrow = Num_two_param, ncol = 2) #��� �ܵõ�˫��Ĳ������ �ľ���
colnames(two_param)=c('beta1','beta0')
for(i in c(1:Num_two_param)) {
  if((ind2[i] %% nrow(PNmat)) == 0) {
    two_param[i,1] = max(pbeta1)
    two_param[i,2] = begin2 + step2 * (ind2[i] %/% nrow(PNmat) - 1)
  }
  else {
    two_param[i,1] = begin1 + step1 * ((ind2[i] %% nrow(PNmat)) - 1) # "%%"��ʾȡ�ࡣ����õ�beta1
    two_param[i,2] = begin2 + step2 * (ind2[i] %/% nrow(PNmat)) # "%/%"��ʾ����������õ�beta0
  }
}

#(beta1,beta0)��ȡֵ��(0.35,0.3),(0.7,0.6),(0.8,0.7)�ǲ�����"beta1С�ڵ���beta0"�ģ�Ҫȥ��
two_param = two_param[-c(1,10,13),]
cor(two_param[,1], two_param[,2], method = "pearson")
cor.test(two_param[,1], two_param[,2], method = "pearson")
lm.sol = lm(two_param[,2] ~ two_param[,1])
summary(lm.sol) #Multiple R-squared:0.8339, p-value:< 2.2e-16

lowAttackRate_Param = read.csv("E:/����/paper/bi-modal/figures/5Periodic-Transmission-Rate/lowAttackRate_Param_b1&b0.csv")
cor.test(lowAttackRate_Param[,1],lowAttackRate_Param[,2],method = "pearson")
lm.sol = lm(lowAttackRate_Param$beta0 ~ lowAttackRate_Param$beta1)
#lm.sol = lm(lowAttackRate_Param$beta1 ~ lowAttackRate_Param$beta0)
summary(lm.sol) #Multiple R-squared:0.994, p-value:1.06e-14


#���浽����
write.csv(PNmat,file = "E:/����/paper/bi-modal/figures/5Periodic-Transmission-Rate/PNmat_b1&b0.csv") #�洢number of waves
write.csv(ATmat, file = "E:/����/paper/bi-modal/figures/5Periodic-Transmission-Rate/ATmat_b1&b0.csv") #�洢attack rate
write.csv(two_param, file = "E:/����/paper/bi-modal/figures/5Periodic-Transmission-Rate/twoWaveParam_b1&b0.csv") #�洢 �ܵõ�˫��Ĳ������