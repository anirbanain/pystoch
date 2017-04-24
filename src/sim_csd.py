import numpy as np
import math

GPStime_start = 1126626073
GPStime_end = 1126712237
segDuration = 26
nSegment = np.int(math.floor((GPStime_end-GPStime_start)/segDuration))

fHigh = 1726
fLow = 20
deltaF = 0.25
nFreqBin = np.int(math.floor((fHigh-fLow)/deltaF) + 1)

#csd = np.random.randn(np.int(nFreqBin),np.int(nSegment))
csd = np.random.randn(nSegment,nFreqBin)


#f = np.arange(fLow, fHigh+deltaF, deltaF)
# mul = f[None,:] * tau[:,:,None]
