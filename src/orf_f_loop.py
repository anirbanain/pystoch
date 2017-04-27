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


nSegment = np.int(np.floor((GPStime_end-GPStime_start)/segDuration))
nFreqBin = np.int(np.floor((fHigh-fLow)/deltaF) + 1)

combined_antenna_response = []
t_delay = []
print 'Calculation of gamma starts'
for t_gps in np.arange(GPStime_start, GPStime_end, segDuration):

    (ifo1_plus, ifo1_cross) = np.vectorize(ifo1.antenna_pattern)(ra, dec, 0, t_gps)
    (ifo2_plus, ifo2_cross) = np.vectorize(ifo2.antenna_pattern)(ra, dec, 0, t_gps)

    combined_antenna_t = (ifo1_plus * ifo2_plus) + (ifo1_cross * ifo2_cross)
    combined_antenna_response.append(combined_antenna_t)

    t_delay_t = np.vectorize(ifo1.time_delay_from_detector)(ifo2,ra,dec,t_gps)
    t_delay.append(t_delay_t)

hp.mollview(t_delay_t, title = "Time Delay Map  High Res")
plt.savefig('./output/Time_Delay_Map.png',dpi = 300)
plt.close('all')
hp.mollview(combined_antenna_t, title = "combined_antenna_t")
plt.savefig('./output/Combined_Antenna_Response.png',dpi = 300)
plt.close('all')
# hp.mollview(np.absolute(plus), title = "Total Plus High Res")
# plt.savefig('./output/Total Plus HR.png',dpi = 200)
# plt.close('all')
# hp.mollview(np.absolute(cross), title = "Total Cross High Res")
# plt.savefig('./output/Total Cross HR.png',dpi = 200)

del theta,phi,ra,dec,ifo1_plus,ifo2_plus,ifo1_cross,ifo2_cross,t_delay_t,combined_antenna_t

end = time.time()
print 'time taken for the overlap reduction function seed matrices',end-start,'sec'


csd = np.random.randn(nSegment,nFreqBin)
t_delay = np.array(t_delay)
combined_antenna_response = np.array(combined_antenna_response)

map_final_mat = []

t = np.arange(GPStime_start, GPStime_end, segDuration)

#f = np.arange(fLow, fHigh+deltaF, deltaF)
#f2 = 2*np.pi*complex(0,1)*f
print 'Total freq bins =', nFreqBin

exp_mat = np.exp(2*np.pi*complex(0,1)*fLow*t_delay)                             #1.65s
exp_df = np.exp(2*np.pi*complex(0,1)*deltaF*t_delay)                            #1.51s
start = time.time()
for ii, f in enumerate(np.arange(fLow, fHigh+deltaF, deltaF)):                  #2.63s
#    phase = 2*np.pi*complex(0,1)*f*t_delay
#    exp_mat = exp_running
    csd_mat = np.transpose(numpy.matlib.repmat(np.transpose(csd[:,ii]),npix,1)) #380ms
    mat_fta = combined_antenna_response * exp_mat * csd_mat                     #1.59s
    exp_mat = exp_mat * exp_df                                                  #415ms
    map_f = np.sum(mat_fta, axis=0)

    if ((ii%50) == 0 and ii != 0):
        print 'upto', f,'Hz done. Time per freq bin =', (time.time()-start)/ii, 'sec'

        hp.mollview(np.absolute(map_f),title='',cbar=False)
        plt.savefig('./output/map_upto_'+str(f).zfill(3)+'_Hz.png',dpi = 300)
        plt.close('all')
        #sys.exit()

    map_final_mat.append(map_f)
    #sys.exit()

print 'All freq bins (',ii+1,') done.'

map_final = np.sum(map_final_mat ,axis=0)

hp.mollview(np.absolute(map_final), title = " ")
plt.savefig('./output/Map.png',dpi = 300)

end = time.time()
print 'total processing and post-processing time',end-start,'sec for nside =',nside
