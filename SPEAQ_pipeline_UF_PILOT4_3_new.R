library("batman", lib.loc="~/R/win-library/3.4")

brukerdata<-readBruker("C:/Users/z072108/Downloads/PILOT4/UF_16082018")
new_bd=brukerdata[-c(rep(1:31586), rep(61825:62539), rep(63126:68002), rep(71189:72684), rep(73595:74765), rep(75416:75870), rep(88161:89461), rep(90763:91737), rep(96615:131072)),]
spectra=as.matrix(t(new_bd))
colnames(spectra)=spectra[1,]
spectra=spectra[-1,]
garrosh=as.data.frame(row.names(spectra))
duplicated(garrosh)
ppm=as.numeric(colnames(spectra))
indx <- grepl('QC', row.names(spectra))
grom=as.numeric(indx)
grom=as.factor(grom)
ENP=list("spectra"=spectra, "ppm"=ppm, "QCs"=grom)

library("speaq", lib.loc="~/R/win-library/3.4")

Y.peaks=getWaveletPeaks(Y.spec = ENP$spectra, X.ppm = ENP$ppm, baselineThresh = 0, SNR.Th = 3, nCPU = -1, include_nearbyPeaks = TRUE)
Y.grouped=PeakGrouper(Y.peaks = Y.peaks, min.samp.grp = 1, grouping.window.width = 200)
Y.grouped_3=Y.grouped[-which(Y.grouped$peakSNR<3),]
Y.filled=PeakFilling(Y.grouped = Y.grouped_3, Y.spec=ENP$spectra, max.index.shift = 8, nCPU = -1)
Y21.filled=Y.filled
ENP_features_new=BuildFeatureMatrix(Y21.filled)
ENP_features_newppm=BuildFeatureMatrix(Y21.filled, var="peakPPM")
nzmean <- function(x) {
  zvals <- x==0
  if(all(zvals)) 0 else mean(x[!zvals])
}
PEAKSHIFT=ENP_features_newppm
PEAKSHIFT2=apply(PEAKSHIFT, 2, nzmean)
PEAKSHIFT=rbind(PEAKSHIFT, (as.numeric(nrow(PEAKSHIFT)))+1)
PEAKSHIFT[nrow(PEAKSHIFT),]=PEAKSHIFT2
colnames(ENP_features_new)=PEAKSHIFT[nrow(PEAKSHIFT),]
row.names(ENP_features_new)=row.names(spectra)
jack=ENP_features_new[,-c(5:9, 20, 21, 25, 82, 83, 128, 129, 136, 138)]   #remove peaks from histidine and citrate AND MA 3.54, METHANOL
will=jack
will[will==0]<-NA
will=will[,-which(colSums(is.na(will))>4)]


setwd("C:/Users/z072108/Downloads/PILOT4/UF_16082018")
write.csv(will, file="UF_PILOT4_3_new.csv")

#total intensity
UF_PILOT4_3_new_tint[UF_PILOT4_3_new_tint==0]<-NA
smith=UF_PILOT4_3_new_tint[,-which(colSums(is.na(UF_PILOT4_3_new_tint))>0)]

write.csv(smith, file="UF_PILOT4_3_for_quartiles.csv")
