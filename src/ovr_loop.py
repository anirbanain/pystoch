import math
import cmath
import time
import healpy as hp
import numpy as np
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
GPStime_end = 1126712237
segDuration = 26
nSegment = math.floor((GPStime_end-GPStime_start)/segDuration)

fHigh = 1726
fLow = 20
deltaF = 0.25
nFreqBin = math.floor((fHigh-fLow)/deltaF) + 1

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
print (end-start)
