import math
import cmath
import time
import healpy as hp
import numpy as np
import numpy.matlib
import pycbc.detector
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

start = time.time()

# Parameters for the High Resolution Map
nside = 1
npix = hp.nside2npix(nside)
pix = np.arange(npix)

(theta, phi) = \
hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))

dec = (np.pi/2) - theta
ra = phi

# The two LIGO detectors
H1 = pycbc.detector.Detector('H1')
L1 = pycbc.detector.Detector('L1')

GPStime_start = 1126626073
#GPStime_end = 1126626173
GPStime_end = 1126712237
segDuration = 26
nSegment = np.int(math.floor((GPStime_end-GPStime_start)/segDuration)+1)

fLow = 20
#fHigh = 22
fHigh = 1726
deltaF = 0.25
nFreqBin = np.int(math.floor((fHigh-fLow)/deltaF) + 1)

combined_antenna_response = []
t_delay = []

for t_gps in np.arange(GPStime_start, GPStime_end, segDuration):

    (H1_plus, H1_cross) = np.vectorize(H1.antenna_pattern)(ra, dec, 0, t_gps)
    (L1_plus, L1_cross) = np.vectorize(L1.antenna_pattern)(ra, dec, 0, t_gps)

    combined_antenna_t = (H1_plus * L1_plus) + (H1_cross * L1_cross)
    combined_antenna_response.append(combined_antenna_t)

    t_delay_t = np.vectorize(H1.time_delay_from_detector)(L1,ra,dec,t_gps)
    t_delay.append(t_delay_t)

end = time.time()
print 'Initial Loop done in', (end-start)
start = time.time()
csd = np.random.randn(nSegment,nFreqBin)
end = time.time()
print 'CSD generated in', (end-start)

f = np.arange(fLow, fHigh+deltaF, deltaF)
start = time.time()
phase = (2.0 * np.pi * complex(0,1)) * np.array(t_delay)[:,:,None] * f[None,:]
end = time.time()
print 'Phase Calculated in ', (end-start)
start = time.time()
csd_mat = np.reshape(numpy.matlib.repmat(csd,1,npix),(nSegment,npix,nFreqBin))
end = time.time()
print 'CSD coverted in matrix in ', (end-start)
start = time.time()
f_factor = np.sum((np.vectorize(cmath.exp)(phase) * csd_mat),axis=2)
end = time.time()
print 'intermediate map calculated in ', (end-start)
map = np.sum(f_factor * combined_antenna_response,axis=0)
end = time.time()
print (end-start)
