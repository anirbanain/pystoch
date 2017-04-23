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
nside_hres = 128
npix_hres = hp.nside2npix(nside_hres)
pix_hres = np.arange(npix_hres)

(theta_hres, phi_hres) = \
hp.pix2ang(nside_hres, np.arange(hp.nside2npix(nside_hres)))

dec_hres = (np.pi/2) - theta_hres
ra_hres = phi_hres

# The two LIGO detectors
H1 = pycbc.detector.Detector('H1')
L1 = pycbc.detector.Detector('L1')

# Segment time and frequency, will be looped later
t_gps = 1126626073
f = 1726

# Plus and Cross Antenna Pattern for the detectors
(H1_plus, H1_cross) = np.vectorize(H1.antenna_pattern)(ra_hres, dec_hres, 0, t_gps)
(L1_plus, L1_cross) = np.vectorize(L1.antenna_pattern)(ra_hres, dec_hres, 0, t_gps)

plus = H1_plus * L1_plus
cross = H1_cross * L1_cross

#overlap reduction function
gamma_Ifta = 4 * np.pi * (plus + cross)* \
np.vectorize(cmath.exp)(2.0 * np.pi * complex(0,1) * f * \
(np.vectorize(H1.time_delay_from_detector)(L1,ra_hres,dec_hres,t_gps)))

#Free Some RAM
H1_plus = None
H1_cross = None
L1_plus = None
L1_cross = None
plus = None
cross = None

#t_delay = np.vectorize(H1.time_delay_from_detector)(L1,ra_hres,dec_hres,t_gps)
#phase = 2.0 * np.pi * complex(0,1) * f * t_delay
#exp_factor = np.vectorize(cmath.exp)(phase)
#gamma_Ifta = 4 * np.pi * (plus + cross) * exp_factor

end = time.time()

print end - start

#hp.mollview(t_delay, title = "Time Delay Map  High Res")
#plt.savefig('../output/Time_Delay_Map.png',dpi = 200)
#plt.close('all')
#hp.mollview(exp_factor.real, title = "exp_factor High Res")
#plt.savefig('../output/exp_factor_HR.png',dpi = 200)
#plt.close('all')
#hp.mollview(np.absolute(plus), title = "Total Plus High Res")
#plt.savefig('../output/Total Plus HR.png',dpi = 200)
#plt.close('all')
#hp.mollview(np.absolute(cross), title = "Total Cross High Res")
#plt.savefig('../output/Total Cross HR.png',dpi = 200)

# Parameters for the Low Resolution Maps
nside2 = 16
npix2 = hp.nside2npix(nside2)
pix2 = np.arange(npix2)

(theta, phi) = hp.pix2ang(nside2, np.arange(hp.nside2npix(nside2)))

#ang_rot  = np.pi/4

earth_rotation = (2 * np.pi) / 86164.1

start = time.time()

ii = 0
for t in np.arange(1126626073, 1126712237, 96):
    pix_rot = np.array(hp.ang2pix(nside2,theta,phi - (earth_rotation * 96 * ii)))
    gamma_Ifta_rot = list(gamma_Ifta[pix_rot])

    hp.mollview(np.absolute(gamma_Ifta_rot),title=" ",cbar=False)
    plt.savefig('../output/gamma_Ifta_rot_'+str(ii).zfill(3)+'.jpg')

    plt.close('all')

    ii = ii + 1

end = time.time()
print (end-start)
