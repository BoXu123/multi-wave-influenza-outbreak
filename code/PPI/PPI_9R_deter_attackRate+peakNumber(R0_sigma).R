#R0��alpha
step1 = 0.25 #����
step2 = 0.025
begin1 = 0 #begin �������㣬��������������
begin2 = 0
pR0 = seq(begin1, 15, by=step1)
psigma = seq(begin2, 1, by=step2)
PNmat = matrix(99, nrow = length(pR0), ncol = length(psigma)) #�洢number of waves
ATmat = matrix(99, nrow = length(pR0), ncol = length(psigma)) #�洢attack rate
epi_time = c(1:100)

ptm = proc.time() #��ʱ��ʼ
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
    ##�ж�˫��
    ##########################����2
    #GT = init.theta[["D_lat"]] + init.theta[["D_inf"]]#generation time
    #GT = init.theta[["D_lat"]] + init.theta[["D_inf"]] + init.theta[["D_imm"]]#generation time
    GT = init.theta[["D_lat"]] + 0.5*init.theta[["D_inf"]]
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
      #combine
      TW = rbind(TWl,TWr)
    }else{
      TW = TWr
    }
    
    #ind1 = which.max(TW[,8]) #��TW[,8]�����ֵ��ֹһ��������ֵ�ظ�������which.max(TW[,8])ֻ�᷵��һ��λ��ֵ��
    #��which(TW[,8] == max(TW[,8]))��ȫ����������λ��ֵ
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
    thre2 = 6 #�˴���Ϊ10�Ļ�������ͼ����6��һ���� #W2�����ֵΪ0.194970448��thre2Ҫ����1/0.194970448=5.13����ʹpeakNumΪ2
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
    
    #Number of peak
    #PN = cbind(PN,peakNum)
    #PNmat[-4+(r01/step1),-4+(r02/step2)] = peakNum
    PNmat[((r0-begin1)/step1)+1,((sig-begin2)/step2)+1] <- peakNum
    
  }
}

proc.time() - ptm #��ʱ����

#plot method 1
source("http://www.phaget4.org/R/myImagePlot.R")
#setwd("E:/����/paper/bi-modal/figures")
myImagePlot(PNmat, xLabels=psigma, yLabels=pR0, title=c("Number of waves"), zlim=c(1,2)) 
#myImagePlot(t(PNmat), xLabels=pC11,yLabels=pA22, title=c("Number of waves"), zlim=c(1,2)) 
myImagePlot(ATmat, xLabels=psigma, yLabels=pR0, title=c("Attack rate"), zlim=c(min(ATmat),max(ATmat))) 

length(which(PNmat==2))
length(which(PNmat==1.8))
length(which(PNmat==1.5))
length(which(PNmat==1.2))
length(which(PNmat==1))

#ȷ��˫���Ӧ�Ĳ������
ind2 = which(PNmat==2) 
#ind2 = which(PNmat==1.8) #����û�з���Ϊ2�ģ�ֻ�з���Ϊ1.8��
Num_two_param = length(ind2) #�ܵõ�˫��Ĳ�����ϵ���Ŀ
two_param = matrix(data = 0, nrow = Num_two_param, ncol = 2) #��� �ܵõ�˫��Ĳ������ �ľ���
colnames(two_param)=c('R0','sigma')
for(i in c(1:Num_two_param)) {
  if((ind2[i] %% nrow(PNmat)) == 0) {
    two_param[i,1] = max(pR0)
    two_param[i,2] = begin2 + step2 * (ind2[i] %/% nrow(PNmat) - 1)
  }
  else {
    two_param[i,1] = begin1 + step1 * ((ind2[i] %% nrow(PNmat)) - 1) # "%%"��ʾȡ�ࡣ����õ�R0
    two_param[i,2] = begin2 + step2 * (ind2[i] %/% nrow(PNmat)) # "%/%"��ʾ����������õ�sigma
  }
}

#���浽����
write.csv(PNmat,file = "E:/����/paper/bi-modal/figures/13Partially Protective Immunity (PPI)/PNmat_R0&sigma.csv") #�洢number of waves
write.csv(ATmat, file = "E:/����/paper/bi-modal/figures/13Partially Protective Immunity (PPI)/ATmat_R0&sigma.csv") #�洢attack rate
write.csv(two_param, file = "E:/����/paper/bi-modal/figures/13Partially Protective Immunity (PPI)/twoWaveParam_R0&sigma.csv") #�洢 �ܵõ�˫��Ĳ������

#########################################################################################################################################
#########################################################################################################################################
######plot epidemic curves (Begin)
sig = 0.35 #�̶�sigΪ0.35��R0ȡֵΪ3.5��5��7.5
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
legend("topleft", legend = c("�� = 0.35"),bty = "n")
legend("topright", legend = c("R0 = 3.5","R0 = 5","R0 = 7.5"),bty = "n", lty = c(1,1,1), lwd = c(2,2,2),col = c("#0000FF","#7F7F80","#FFFF00") )
legend(x=60, y=4000, legend = c("PN=1","PN=1.5","PN=2"),bty = "n")
##PNΪ1�Ķ�ӦRGB��ɫΪ(0,0,255),PNΪ1.2�Ķ�ӦRGB��ɫΪ(51,51,204),PNΪ1.5�Ķ�ӦRGB��ɫΪ(127,127,128)��PNΪ1.8�Ķ�ӦRGB��ɫΪ(204,204,51),PNΪ2�Ķ�ӦRGB��ɫΪ(255,255,0)
##PNΪ1�Ķ�Ӧ��ɫΪ"#0000FF"��PNΪ1.2�Ķ�Ӧ��ɫΪ"#3333CC"��PNΪ1.5�Ķ�Ӧ��ɫΪ"#7F7F80"��PNΪ1.8�Ķ�Ӧ��ɫΪ"#CCCC33"��PNΪ2�Ķ�Ӧ��ɫΪ"#FFFF00"
######plot epidemic curves (End)