install.packages(c("MALDIquant", "MALDIquantForeign",
                   "sda", "crossval"))
install.packages("sda")
install.packages('crossval')
setwd("/home/vdarbot/Bureau/metabolo/")
library(MALDIquant)
library(MALDIquantForeign)
library(sda)
library(crossval)
directory_data <- "Converted"
data(fiedler2009subset)
fiedler2009subset
## import
files <- importMzMl(directory_data)
sample_names <- gsub(list.files(directory_data),pattern="\\.mzml$", replacement="")
metaData(fiedler2009subset)
## check raw data
## any empty spectra? (empty spectra are ignored in subsequent baseline
## correction/peak detection; you could find/remove them by calling
## findEmptyMassObjects/removeEmptyMassObjects)
## see ?isEmpty, ?findEmptyMassObjects, ?removeEmptyMassObjects
any(sapply(files, isEmpty))
# FALSE

## do length of spectra differ? (if they differ you have to adjust the
## corresponding halfWindowSize in subsequent baseline correction and peak
## detection.)
table(sapply(files, length))
# 42388
#    16

## import

## preprocessing
## sqrt transform (for variance stabilization)
spectra_sqrt <- transformIntensity(files, method="sqrt")
plot(spectra_sqrt[[1]])
spectra_smooth <- smoothIntensity(spectra_sqrt, method="SavitzkyGolay", halfWindowSize=20)
plot(spectra_smooth[[1]])
spectra_baseline <- removeBaseline(spectra_smooth, method="SNIP", iterations=100)
plot(spectra_baseline[[1]])
## calibrate (normalize) intensities (different calibration methods available)
## see ?calibrateIntensity
spectra_calibr <- calibrateIntensity(spectra_baseline, method="TIC")
plot(spectra_calibr[[1]])
## spectra alignment
## (the spectra alignment is peak based, maybe you need to adjust
## halfWindowSize, SNR, tolerance, warpingMethod)
## see ?alignSpectra
spectra_align <- alignSpectra(spectra_calibr,
                        halfWindowSize=20, SNR=2,
                        tolerance=0.002, warpingMethod="lowess")
plot(spectra_align[[1]])

metaData(spectra_align[[1]])
sample_names
strain_1 <- grepl(pattern="200064", x=sample_names)
strain_2 <- grepl(pattern="220020", x=sample_names)
BTS <- grepl(pattern="BTS", sample_names)



classes <- factor(ifelse(strain_1, "200064", ifelse(strain_2, "220020", "BTS")),
                  levels=c("200064", "220020", "BTS"))
## 2. average technical replicates
## see ?averageMassSpectra
avgSpectra <- averageMassSpectra(spectra_align, labels= classes, method="mean")
avgSpectra
metaData()
## run peak detection
## (maybe you need to adjust halfWindowSize [decreasing it for high resolution
## spectra] and SNR [a higher value increase the True-Positive-Rate but decrease
## sensitivity])
## see ?detectPeaks, ?estimateNoise
peaks <- detectPeaks(avgSpectra, method="MAD", halfWindowSize=20, SNR=2)
plot(peaks[[1]])

## bin peaks
## (After alignment peak positions (mass) are similar but not identical. Binning
## is needed to make similar peak mass values identical.)
## see ?binPeaks
peaks <- binPeaks(peaks, tolerance=0.002)



## 2. export expression/training matrix
## (and fill missing peaks by interpolated values)
## see ?intensityMatrix
featureMatrix <- intensityMatrix(peaks, avgSpectra)

rownames(featureMatrix) <-  unique(classes)
colnames(featureMatrix) <-
  round(as.double(colnames(featureMatrix)),2)
classes <- factor(classes)
ddar <- sda.ranking(Xtrain=featureMatrix, L=sample_names, fdr=FALSE, diagonal=TRUE, verbose=FALSE)
ddar
classes

install.packages("pvclust")
library(pvclust)
pv <- pvclust(t(featureMatrix),
              method.hclust="ward.D2",
              method.dist="euclidean")
plot(pv, print.num=FALSE)

# Next, we repeat the above clustering on the data set containing only the best two
# top-ranking peaks:
  top <- ddar[1:2, "idx"]
distanceMatrixTop <- dist(featureMatrix[, top],
                          method="euclidean")
hClustTop <- hclust(distanceMatrixTop, method="complete")
plot(hClustTop, hang=-1)
?(binPeaks)
##
predfun.dda <- function(Xtrain, Ytrain, Xtest, Ytest, negative) {
  dda.fit <- sda(Xtrain, Ytrain, diagonal=TRUE, verbose=FALSE)
  ynew <- predict(dda.fit, Xtest, verbose=FALSE)$class
  return(confusionMatrix(Ytest, ynew, negative=negative))
}
# set seed to get reproducible results
set.seed(1234)
cv.out <- crossval(predfun.dda,
                   X=featureMatrix[, top],
                   Y=classes,
                   K=10, B=20,
                   negative="control",
                   verbose=FALSE)
diagnosticErrors(cv.out$stat)

data(singh2002)

# training data
Xtrain = singh2002$x
Ytrain = singh2002$y

######################################### 
# feature ranking (diagonal covariance) #
#########################################

# ranking using t-scores (DDA)
ranking.DDA = sda.ranking(Xtrain, Ytrain, diagonal=TRUE)
