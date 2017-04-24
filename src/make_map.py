import numpy as np
import math
import cmath

GPStime_start = 1126626073
GPStime_end = 1126712237
segDuration = 26
nSegment = np.int(math.floor((GPStime_end-GPStime_start)/segDuration))

fHigh = 1726
fLow = 20
deltaF = 0.25
nFreqBin = np.int(math.floor((fHigh-fLow)/deltaF) + 1)

f = np.arange(fLow, fHigh+deltaF, deltaF)

#for t_gps in np.arange(GPStime_start, GPStime_end, segDuration):
#    for f in np.range(fLow, fHigh, deltaF):
#        cmath.exp(2.0 * np.pi * complex(0,1) * f * t_delay)

#cmath.exp(2.0 * np.pi * complex(0,1) * f * t_delay)

phase = (2.0 * np.pi * complex(0,1)) * np.array(t_delay)[:,:,None] * f[None,:]

csd_mat = np.reshape(numpy.matlib.repmat(csd,1,npix),(nSegment,npix,nFreqBin))
f_factor = np.sum(cmath.exp(phase) * csd_mat,axis=2)

map = np.sum(f_factor * combined_antenna_response,axis=0)
