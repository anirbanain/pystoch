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

param_input = open("../input/param_all.txt","r")

ifo1 = "H1"
ifo2 = "L1"

nside = 1

GPStime_start = 1126626073
#GPStime_end = 1126626173
GPStime_end = 1126712237
segDuration = 26

fLow = 20
#fHigh = 21
fHigh = 1726
deltaF = 0.25

# Parameters for the High Resolution Map

npix = hp.nside2npix(nside)
pix = np.arange(npix)

(theta, phi) = \
hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))

dec = (np.pi/2) - theta
ra = phi

# The two LIGO detectors
ifo1 = pycbc.detector.Detector(ifo1)
ifo2 = pycbc.detector.Detector(ifo2)


nSegment = np.int(math.floor((GPStime_end-GPStime_start)/segDuration)+1)
nFreqBin = np.int(math.floor((fHigh-fLow)/deltaF) + 1)

combined_antenna_response = []
t_delay = []

for t_gps in np.arange(GPStime_start, GPStime_end, segDuration):

    (ifo1_plus, ifo1_cross) = np.vectorize(ifo1.antenna_pattern)(ra, dec, 0, t_gps)
    (ifo2_plus, ifo2_cross) = np.vectorize(ifo2.antenna_pattern)(ra, dec, 0, t_gps)

    combined_antenna_t = (ifo1_plus * ifo2_plus) + (ifo1_cross * ifo2_cross)
    combined_antenna_response.append(combined_antenna_t)

    t_delay_t = np.vectorize(ifo1.time_delay_from_detector)(ifo2,ra,dec,t_gps)
    t_delay.append(t_delay_t)

end = time.time()
print 'time taken for the overlap reduction function seed matrices',end-start,'sec'
start = time.time()


csd = np.random.randn(nSegment,nFreqBin)
t_delay = np.array(t_delay)
combined_antenna_response = np.array(combined_antenna_response)

# for f in np.arange(fLow, fHigh+deltaF, deltaF):
#     for t in np.arange(GPStime_start, GPStime_end, segDuration):
#         print t
#         print f

map_final_mat = []

ii = 0
f = 2*np.pi*(np.arange(fLow, fHigh+deltaF, deltaF))
for t in np.arange(GPStime_start, GPStime_end, segDuration):

    #exp_term = np.vectorize(cmath.exp)\
    #(2*np.pi*complex(0,1)*(t_delay[ii,:][:,None] * f[None,:]))

    #thisphase = 2*np.pi*(t_delay[ii,:][:,None] * f[None,:])
    np.shape(t_delay)
    np.shape(f)
    thisphase = np.outer(t_delay, f)
    sine_mat = np.vectorize(math.sin)(thisphase)
    cos_mat = np.vectorize(math.cos)(thisphase)
    exp_term =  np.vectorize(complex)(cos_mat,sine_mat)
    #(2*np.pi*complex(0,1)*(t_delay[ii,:][:,None] * f[None,:]))

    f_mat = exp_term * (numpy.matlib.repmat(csd[ii,:],npix,1))
    map_t =  np.sum(f_mat, axis=1)

    ii =ii+1
    if (ii%100) == 0:
        #print (100*(ii/nSegment)),'% Done'
        print ii, 'segments done in', time.time()

    map_final_mat.append(map_t)

print ii, 'segments done'
map_final = np.sum(map_final_mat,axis=0)

end = time.time()
print 'total processing and post-processing time',end-start,'sec for nside =',nside

#phase = (2.0 * np.pi * complex(0,1)) * np.array(t_delay)[:,:,None] * f[None,:]
#csd_mat = np.reshape(numpy.matlib.repmat(csd,1,npix),(nSegment,npix,nFreqBin))
#f_factor = np.sum((np.vectorize(cmath.exp)(phase) * csd_mat),axis=2)
#map = np.sum(f_factor * combined_antenna_response,axis=0)
