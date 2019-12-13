#load batman from r 3.4

brukerdata<-readBruker("C:/Users/z072108/Downloads/PILOT4/CPMG_16082018")
new_bd=brukerdata[-c(rep(1:31599), rep(61838:62552), rep(63138:68014), rep(71202:72696), rep(73608:74777), rep(75429:75883), rep(88174:89474), rep(90775:91750), rep(96628:131072)),]
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

#load speaq from r 3.6

Y.peaks=getWaveletPeaks(Y.spec = ENP$spectra, X.ppm = ENP$ppm, baselineThresh = 0, SNR.Th = 3, nCPU = -1, include_nearbyPeaks = TRUE)
Y.grouped=PeakGrouper(Y.peaks = Y.peaks, min.samp.grp = 1, grouping.window.width = 200)
Y.grouped_3=Y.grouped[-which(Y.grouped$peakSNR<3),]
Y.filled=PeakFilling(Y.grouped = Y.grouped_3, Y.spec=ENP$spectra, max.index.shift = 8, nCPU = -1)
ENP_features_new=BuildFeatureMatrix(Y.filled)
ENP_features_newppm=BuildFeatureMatrix(Y.filled, var="peakPPM")
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
#jack=ENP_features_new[,-c(75, 74, 72, 71, 47, 48, 8)]   #remove peaks from histidine, citrate, maleic acid, 3.54 dummy

setwd("C:/Users/z072108/Downloads/PILOT4/CPMG_16082018")
write.csv(ENP_features_new, file="CPMG_PILOT4_3.csv")
