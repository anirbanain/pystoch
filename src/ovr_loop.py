import sys
import time
import healpy as hp
import numpy as np
import numpy.matlib
import pycbc.detector
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ConfigParser

start = time.time()

config = ConfigParser.ConfigParser()
config.read("./input/parameter_all.ini")

nside = config.getint('parameters','nside')
ifo1 = config.get('parameters','ifo1')
ifo2 = config.get('parameters','ifo2')
GPStime_start = config.getint('parameters','GPStime_start')
GPStime_end = config.getint('parameters','GPStime_end')
segDuration = config.getint('parameters','segDuration')
fLow = config.getint('parameters','fLow')
fHigh = config.getint('parameters','fHigh')
deltaF = config.getfloat('parameters','deltaF')

npix = hp.nside2npix(nside)
pix = np.arange(npix)

(theta, phi) = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))

dec = (np.pi/2) - theta
ra = phi

# The two LIGO detectors
ifo1 = pycbc.detector.Detector(ifo1)
ifo2 = pycbc.detector.Detector(ifo2)


nSegment = np.int(np.floor((GPStime_end-GPStime_start)/segDuration)+1)
nFreqBin = np.int(np.floor((fHigh-fLow)/deltaF) + 1)

combined_antenna_response = []
t_delay = []

for t_gps in np.arange(GPStime_start, GPStime_end, segDuration):

    (ifo1_plus, ifo1_cross) = np.vectorize(ifo1.antenna_pattern)(ra, dec, 0, t_gps)
    (ifo2_plus, ifo2_cross) = np.vectorize(ifo2.antenna_pattern)(ra, dec, 0, t_gps)

    combined_antenna_t = (ifo1_plus * ifo2_plus) + (ifo1_cross * ifo2_cross)
    combined_antenna_response.append(combined_antenna_t)

    t_delay_t = np.vectorize(ifo1.time_delay_from_detector)(ifo2,ra,dec,t_gps)
    t_delay.append(t_delay_t)

del theta,phi,ra,dec,ifo1_plus,ifo2_plus,ifo1_cross,ifo2_cross,t_delay_t,combined_antenna_t

end = time.time()
print 'time taken for the overlap reduction function seed matrices',end-start,'sec'


csd = np.random.randn(nSegment,nFreqBin)
t_delay = np.array(t_delay)
combined_antenna_response = np.array(combined_antenna_response)

map_final_mat = []

f = np.arange(fLow, fHigh+deltaF, deltaF)
f2 = 2*np.pi*complex(0,1)*f

start = time.time()
for ii, t in enumerate(np.arange(GPStime_start, GPStime_end, segDuration)):

    tf_mat = t_delay[ii,:][:,None] * f2[None,:]
    exp_term = np.exp(tf_mat)

    f_mat = exp_term * (numpy.matlib.repmat(csd[ii,:],npix,1))
    map_t =  np.sum(f_mat, axis=1)

    if ((ii%100) == 0 and ii != 0):
        print ii, 'segments done. Time per segment =', (time.time()-start)/ii, 'sec'
        #sys.exit()

    map_final_mat.append(map_t)
#    sys.exit()

print 'All segments (',ii+1,') done.'

map_final = np.sum(map_final_mat * combined_antenna_response,axis=0)

del exp_term,f_mat,map_t

end = time.time()
print 'total processing and post-processing time',end-start,'sec for nside =',nside
