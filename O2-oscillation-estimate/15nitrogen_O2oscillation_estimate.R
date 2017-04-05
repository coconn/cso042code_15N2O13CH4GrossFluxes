# 15nitrogen_O2oscillation_estimate.R
# 
# for the 15N2O lab experiment, what anoxic-oxic oscillation treatment should we use?
# determine the most common length of fluctuation between oxygenated and not
#
# O2 - redox - GHG project
# CS O'Connell, UCB, Silver Lab
#
# see also: Combine-array-sensor-data.R, which points to where the O2 csv file lives


########################################################################
# BRING IN DATA

library(ggplot2)
library(plyr)
library(dplyr)
library(lubridate)

# summarySE using plyr
source("~/Documents/GITHUB/RPersonalFunctionsChristine/summarySE.r")

# where does the O2 csv live
pathdata = "~/Documents/GITHUB/cso044code_HotSpotsHotMoments/HotSpotsHotMomentsAnalysis/HotSpotsHotMoments-Data-Raw/Sensors/SurfaceDataArchive/"
# where to save figures
pathsavefigs = "~/Documents/GITHUB/cso044code_HotSpotsHotMoments/HotSpotsHotMoments15N/"

# bring in O2 data
O2hourly <- as.data.frame(read.csv(paste(pathdata,"O2hourly",".csv",sep=""), stringsAsFactors=FALSE))


########################################################################
# FORMAT AND CLEAN DATA

# make factors where needed
O2hourly$SensorID <- as.factor(O2hourly$SensorID)
O2hourly$TopoLocation <- as.factor(O2hourly$TopoLocation)

# get rid of any O2 values above 0.23
O2hourly$O2pct[O2hourly$O2pct>22.5] <- NA

# deal with date formatting
O2hourly$TIMESTAMP2 <- ymd_hms(O2hourly$TIMESTAMP2)

# subset working for dates?
datesub <- O2hourly[O2hourly$TIMESTAMP2 > ymd(20160601),]
str(datesub)
# yup, working; all good


########################################################################
# EXPLORATORY GRAPHS
# to take a look

ggplot(datesub[datesub$SensorID=="Sensor06",],aes(x=TIMESTAMP2,y=O2pct)) + 
  geom_point() +
  labs(x="Date",y="O2 concentration")

png(file = paste(pathsavefigs, "O2sensors_oscillations.png", sep=""),width=10,height=10,units="in",res=400)
ggplot(datesub,aes(x=TIMESTAMP2,y=O2pct)) + 
  geom_point() +
  labs(x="Date",y="O2 concentration") + 
  facet_wrap(~SensorID, nrow = 5, ncol = 7)
dev.off()


########################################################################
# WEIRD O2 SPIKE IN OCT 2016

datesub2 <- O2hourly[O2hourly$TIMESTAMP2 > ymd(20161010) & O2hourly$TIMESTAMP2 < ymd(20161030),]

ggplot(datesub2[datesub2$SensorID=="Sensor06",],aes(x=TIMESTAMP2,y=O2pct)) + 
  geom_point() +
  labs(x="Date",y="O2 concentration")

png(file = paste(pathsavefigs, "O2sensors_oscillations_zoomed.png", sep=""),width=10,height=10,units="in",res=400)
ggplot(datesub2,aes(x=TIMESTAMP2,y=O2pct)) + 
  geom_point() +
  labs(x="Date",y="O2 concentration") + 
  facet_wrap(~SensorID, nrow = 5, ncol = 7)
dev.off()

png(file = paste(pathsavefigs, "O2sensors_Sensor06_goodex.png", sep=""),width=4,height=3,units="in",res=400)
ggplot(datesub2[datesub2$SensorID=="Sensor06",],aes(x=TIMESTAMP2,y=O2pct)) + 
  geom_point() +
  labs(x="Date",y="O2 concentration")
dev.off()


########################################################################
# ATTEMPT PERIODOGRAM

# see https://onlinecourses.science.psu.edu/stat510/node/71

# sort by time
out <- datesub2[order(datesub2$TIMESTAMP2),]
# only one sensor
out2 <- out[out$SensorID=="Sensor06",]
# there are weirdly hella duplicates (does the CR1000 not record minutes?)
out3 <- out2[!duplicated(out2$TIMESTAMP2),]

# only time and O2
out4 <- out3[,c(2,11)]

# time difference
out4$TIMESTAMP2[1] - out4$TIMESTAMP2[dim(out4)[1]]
timepds <- difftime(out4$TIMESTAMP2[dim(out4)[1]], out4$TIMESTAMP2[1], units="hours")
as.numeric(timepds) # Time difference of 458 hours

# time series context
# O2 percent, measured every 1 hour for 458 hours (with missing data)

xB <- out4[,2]
# NAs are a problem
xB[is.na(xB)] <- 0
# back to suggested code
FF <- abs(fft(xB)/sqrt(458))^2
howmany <- 458/2+1
howmany2 <- 458/2
P <- (4/458)*FF[1:howmany] # Only need the first (n/2)+1 values of the FFT result.
f <- (0:howmany2)/458 # this creates harmonic frequencies from 0 to .5 in steps of 1/128.
plot(f, P, type="l") # This plots the periodogram; type = “l” creates a line plot.  Note: l is lowercase L, not number 1.
# this looks incredibly fucked up

# what about the other example?
# sunspots=scan("sunspots.dat")
# plot(sunspots,type="b")
# x = diff(sunspots)
# I = abs(fft(x)/sqrt(458))^2
# P = (4/458)*I[1:230]
# freq = (0:229)/458
# plot(freq,P,type="l")

# try with my data
O2 <- out4[,2]
O2[is.na(O2)] <- 9 # note that this is a really fucked up way to solve this problem... I need a better solution since this only works for this one sensor
# instead, take the average of the points around the NAs or something
plot(O2,type="b")
x = diff(O2)
#x = O2 to see if detrending helps or not
plot(x,type="b")
I = abs(fft(x)/sqrt(458))^2
P = (4/458)*I[1:230]
freq = (0:229)/458
plot(freq,P,type="l")
# miraculously, this does not look fucked up

# where is the peak in the prob plot
indexmax <- which(P==max(P))
freqmax <- freq[indexmax]
periodtmp <- 1/freqmax
periodtmp/24

# so for this short time period, for sensor 6, the period is ~90 hours, or 3.8 days


# reference code from https://onlinecourses.science.psu.edu/stat510/node/71
# 
# # Here’s the code.  The value 128 in a few of the lines is the sample size.  The fast Fourier Transform (fft) includes a frequency value of 0 and goes all the way to 1.  The second line, creates the basis for what our textbook authors call the unscaled periodogram.  We only need to go to 0.5 for the periodogram so in the third line below we pick off the first (n/2)+1 = 65 elements of FF (the object created in the second line).  Also in the third line, the multiplication by 4/128 (=4/n) creates the scaled periodogram as described on pages 68-70 of the book.  In a sense, this scaling is optional because we would see the same pattern even if we did not use this constant.
# x = scan("cortex.dat")
# FF = abs(fft(x)/sqrt(128))^2
# P = (4/128)*FF[1:65] # Only need the first (n/2)+1 values of the FFT result.
# f = (0:64)/128 # this creates harmonic frequencies from 0 to .5 in steps of 1/128.
# plot(f, P, type="l") # This plots the periodogram; type = “l” creates a line plot.  Note: l is lowercase L, not number 1.
# 
# sunspots=scan("sunspots.dat")
# plot(sunspots,type="b")
# x = diff(sunspots)
# I = abs(fft(x)/sqrt(458))^2
# P = (4/458)*I[1:230]
# freq = (0:229)/458
# plot(freq,P,type="l")
# 


