#remove peaks not present in at least 80% of samples in each group
setwd("D:/MEGA/Publication 1/PILOT4 (2)/LED_UR_16082018/")
library(readr)
UF_PILOT4_3_withexcl <- read_csv("D:/MEGA/Publication 1/PILOT4 (2)/LED_UR_16082018/LED_UR_PILOT4_3_excl_final.csv")
jack=UF_PILOT4_3_withexcl
require("plyr")
count(jack, c("GROUP"))
pearl=jack[,-which(colSums(is.na(jack[which(jack$GROUP=="PHEO"),]))>=2 & colSums(is.na(jack[which(jack$GROUP=="PA"),]))>=2)]
#write.csv(file="LED_PILOT4_3_excl_new_20.csv", pearl)

#scale to MA peak (5.99)
captain=pearl
captain[is.na(captain)]=0
captain=captain[,order(colnames(captain))]
captain=cbind(captain$SAMPLE, captain$GROUP, captain[,1:(as.numeric(ncol(captain))-2)])
#captain$`5.99` <- NULL
#captain=cbind(captain, rowSums(captain[3:(as.numeric(ncol(captain)))]))
#sparrow=cbind(captain[,1:2], captain[,3:(as.numeric(ncol(captain)))]/captain[,(as.numeric(ncol(captain)))])
#sparrow=sparrow[,-ncol(sparrow)]
sparrow=cbind(captain[,1:2], captain[,3:(as.numeric(ncol(captain)))]/captain$`5.99`)
sparrow$`5.995043486`<- NULL
#write.csv(file="LED_PILOT4_3_excl_new_20_MA.csv", sparrow)

#PQN
black=sparrow[,3:ncol(sparrow)]
black[black==0]=NA
dauntless=apply(black, 2, median, na.rm=T)
check=t(t(black)/dauntless)
check2=apply(check, 1, median, na.rm=T)
checkmate=black/check2
checkmate[is.na(checkmate)]=0
beckett=cbind(sparrow[,1:2], checkmate)
#write.csv(file="LED_PILOT4_3_excl_new_20_MA_PQN0.csv", beckett)

#kNN
library(impute)
ncheck=checkmate
ncheck[ncheck==0]=NA
ncheckO.imputed <- impute.knn(as.matrix(ncheck), k=5, rowmax = 1, colmax = 1, maxp=nrow(ncheck))
cutler=cbind(beckett[,1:2], ncheckO.imputed$data)
#write.csv(file="LED_PILOT4_3_excl_new_20_MA_PQN0_knn.csv", cutler)

#glog
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
#write.csv(final, file="LED_PILOT4_3_excl_new_20_MA_PQN0_knn_GLOG.csv")

#MixOmics models
PLASMA_GLOG_forMIXOMIX = final
ex=matrix(as.numeric(unlist(PLASMA_GLOG_forMIXOMIX[,-(1:2)])), nrow=nrow(PLASMA_GLOG_forMIXOMIX))
colnames(ex)=colnames(PLASMA_GLOG_forMIXOMIX[,-(1:2)])
row.names(ex)=PLASMA_GLOG_forMIXOMIX$SAMPLE
why=as.factor(PLASMA_GLOG_forMIXOMIX$`captain$GROUP`)
#ex=RFmarkerDetector::paretoscale(ex, exclude=F)
library(mixOmics)
innames=as.numeric(colnames(ex))
innames=round(innames, 3)
innames=as.character(innames)
colnames(ex)=innames
pca.ENP = pca(ex, ncomp = 2, center = TRUE, scale = FALSE)
plotIndiv(pca.ENP, group = why, ind.names = FALSE,legend = T, title= 'PCA LED SPEAQ MA PQN GLOG')
plotLoadings(pca.ENP, comp = 1, ndisplay = 30, title= 'PCA LED SPEAQ MA PQN GLOG')
ENP.plsda <- plsda(ex, why, ncomp = 10, scale = F)
perf.plsda <- perf(ENP.plsda, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 25, dist = "mahalanobis.dist")
ENP.plsda <- plsda(ex, why, ncomp = min(which(perf.plsda$error.rate$BER[,1]==min(perf.plsda$error.rate$BER[,1]))), scale = F)

#check for correlation with sequence (stability-degradation)
seq_check=match(row.names(ex), UF_PILOT4_3_withexcl$SAMPLE)
mitch=predict(ENP.plsda, ex)
predictions_uncorrected=as.matrix(mitch$class$mahalanobis.dist[,2])
predictions_uncorrected[predictions_uncorrected=="PA"]=1
predictions_uncorrected[predictions_uncorrected=="PHEO"]=2
predictions_uncorrected=as.matrix(predictions_uncorrected)
predictions_uncorrected=as.numeric(predictions_uncorrected)
predictions_uncorrected=as.matrix(predictions_uncorrected)
cor.test(seq_check, predictions_uncorrected)
seq_check[seq_check<9]=1
seq_check[seq_check==9]=1
seq_check[seq_check>9]=2
cor.test(seq_check, predictions_uncorrected)
#no correlations, same holds for rest of models, including AMIX!
mitch=vip(ENP.plsda)
write.csv(mitch, file="LED_PILOT4_3_excl_new_20_MA_PQN0_knn_GLOG_VIPS.csv")

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

save.image("D:/MEGA/Publication 1/PILOT4 (2)/LED_UR_16082018/LED_SPEAQ_P4.RData")