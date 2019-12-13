library(readr)
setwd("D:/MEGA/Publication 1/PILOT4 (2)/btp_final_2/btp_final/MA5/")
UF_PILOT4_3_new_25_PQN_min_noRSD0_3_GLOG <- read_csv("D:/MEGA/Publication 1/PILOT4 (2)/btp_final_2/btp_final/MA5/bucket_table_led_new.csv")
UF_PILOT4_3_withexcl=UF_PILOT4_3_new_25_PQN_min_noRSD0_3_GLOG[order(UF_PILOT4_3_new_25_PQN_min_noRSD0_3_GLOG$GROUP),]

jack=UF_PILOT4_3_withexcl
require("plyr")
count(jack, c("GROUP"))
jack[jack==0]=NA
pearl=jack[,-which(colSums(is.na(jack[which(jack$GROUP=="PHEO"),]))>=2 & colSums(is.na(jack[which(jack$GROUP=="PA"),]))>=2)]
#write.csv(file="bucket_table_led_new_20.csv", pearl)

#scale to MA peak (5.99)
captain=pearl
captain[is.na(captain)]=0
captain=captain[,order(colnames(captain))]
captain=cbind(captain$X1, captain$GROUP, captain[,1:(as.numeric(ncol(captain))-2)])
#captain$`5.99`=NULL
#captain=cbind(captain, rowSums(captain[3:(as.numeric(ncol(captain)))]))
#sparrow=cbind(captain[,1:2], captain[,3:(as.numeric(ncol(captain)))]/captain[,(as.numeric(ncol(captain)))])
#sparrow=sparrow[,-ncol(sparrow)]
sparrow=cbind(captain[,1:2], captain[,3:(as.numeric(ncol(captain)))]/captain$`5.99`)
sparrow$`5.99`<- NULL
#write.csv(file="bucket_table_led_new_20_MA.csv", sparrow)

#PQN
black=sparrow[,3:ncol(sparrow)]
dauntless=apply(black, 2, median)
check=t(t(black)/dauntless)
check2=apply(check, 1, median)
checkmate=black/check2
beckett=cbind(sparrow[,1:2], checkmate)
cutler=beckett
#write.csv(file="bucket_table_led_new_20_MA_PQN.csv", beckett)



library("LMGene")
library(Biobase)
library(tools)
library(readr)

ENP_vanilla_forGLOG=cutler
PLASMA_QC_PQN=ENP_vanilla_forGLOG[which(ENP_vanilla_forGLOG$`captain$GROUP`=="PA"),]
QCled_PQN=as.matrix(t(PLASMA_QC_PQN))
QCled_PQN=QCled_PQN[-2,]
colnames(QCled_PQN)=QCled_PQN[1,]
QCled_PQN=QCled_PQN[-1,]
QCr=apply(QCled_PQN, 1,as.numeric)
QCr=t(QCr)
QCmonster=as.factor(c(rep(1:2, each=4), 1))
QCdose=as.numeric(c(rep(1, times=9)))
QCled_list=list("monster"=QCmonster, "dose"=QCdose)
QCled.eS=neweS(QCr, QCled_list)
tranpar <- tranest(QCled.eS)
tranpar

PLASMA_PQN=ENP_vanilla_forGLOG
led_PQN=as.matrix(t(PLASMA_PQN))
led_PQN=led_PQN[-2,]
colnames(led_PQN)=led_PQN[1,]
led_PQN=led_PQN[-1,]
r=apply(led_PQN, 1,as.numeric)
r=t(r)
monster=as.factor(c(1:18))
dose=as.numeric(c(rep(1, times=18)))
led_list=list("monster"=monster, "dose"=dose)
led.eS=neweS(r, led_list)
trled.eS <- transeS(led.eS, tranpar$lambda, tranpar$alpha)
kostakis=exprs(trled.eS)
colnames(kostakis)=as.factor(colnames(kostakis))
colnames(kostakis)=colnames(led_PQN)
final=t(kostakis)
final=cbind(cutler[,1:2], final[,1:as.numeric(ncol(final))])
final=final[,-1]
final=cbind(UF_PILOT4_3_withexcl[,1], final)
final=final[order(final[,2]),]
#write.csv(final, file="bucket_table_led_new_20_MA_PQN_GLOG.csv")

#MixOmics models
PLASMA_GLOG_forMIXOMIX = final
ex=matrix(as.numeric(unlist(PLASMA_GLOG_forMIXOMIX[,-(1:2)])), nrow=nrow(PLASMA_GLOG_forMIXOMIX))
colnames(ex)=colnames(PLASMA_GLOG_forMIXOMIX[,-(1:2)])
row.names(ex)=PLASMA_GLOG_forMIXOMIX$X1
why=as.factor(PLASMA_GLOG_forMIXOMIX$`captain$GROUP`)
#ex=RFmarkerDetector::paretoscale(ex, exclude=F)
library(mixOmics)
innames=as.numeric(colnames(ex))
innames=round(innames, 2)
innames=as.character(innames)
colnames(ex)=innames
pca.ENP = pca(ex, ncomp = 2, center = TRUE, scale = FALSE)
plotIndiv(pca.ENP, group = why, ind.names = FALSE,legend = T, title= 'PCA LED AMIX MA PQN GLOG')
plotLoadings(pca.ENP, comp = 1, ndisplay = 30, title= 'PCA LED AMIX MA PQN GLOG')
ENP.plsda <- plsda(ex, why, ncomp = 10, scale = F)
perf.plsda <- perf(ENP.plsda, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 25, dist = "mahalanobis.dist")
ENP.plsda <- plsda(ex, why, ncomp = min(which(perf.plsda$error.rate$BER[,1]==min(perf.plsda$error.rate$BER[,1]))), scale = F)
mitch=vip(ENP.plsda)
write.csv(mitch, file="LED_AMIX_MA_PQN_GLOG_VIPs.csv")

#CV2
x <- ex
X <- x
Y <- why

hay=function(x) {
  {
    samp12=x
    test <- samp12
    train <- setdiff(1:nrow(X), test)
    plsda.train <- plsda(X[train, ], Y[train], ncomp = 10, scale = F)
    perf.plsda.train <- perf(plsda.train, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 25, dist = "mahalanobis.dist")
    plsda.train <- plsda(X[train, ], Y[train], ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), scale = F)
    test.predict <- predict(plsda.train, t(as.matrix(X[test,])), dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1])))]                         #number of components
    well=Y[test]==Prediction
  }
  return(well)
}

NC=lapply((1:18), hay)          #number of repeat
NC=unlist(NC)
