import math
import healpy as hp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pycbc.detector

#plt.ioff()

nside_hres = 256
npix_hres = hp.nside2npix(nside_hres)
pix_hres = np.arange(npix_hres)

(theta_hres, phi_hres) = hp.pix2ang(nside_hres, np.arange(hp.nside2npix(nside_hres)))

dec_hres = (np.pi/2) - theta_hres
ra_hres = phi_hres

H1 = pycbc.detector.Detector('H1')

(H1_plus_hres, H1_cross_hres) = np.vectorize(H1.antenna_pattern)(ra_hres, dec_hres, 0, 1126626073)

H1_plus_hres = np.array(H1_plus_hres)

#pycharm edit
#fig = plt.figure()
hp.mollview(np.absolute(H1_plus_hres), title = "H1 Plus High Res")
plt.savefig('./H1_Plus_HR.png')

H1_cross_hres = np.array(H1_cross_hres)

hp.mollview(np.absolute(H1_cross_hres), title = "H1 Cross High Res")
plt.savefig('./H1_Cross_HR.png')


#nside = 4
#npix = hp.nside2npix(nside)
#pix = np.arange(npix)

#(theta, phi) = hp.pix2ang(nside, np.arrange(hp.nside2npix(nside)))

#dec = (np.pi/2) - theta
#ra = phi

#ang_rot = np.pi/8

#pix_rot = np.array(hp.ang2pix(nside_hres, theta, phi - ang_rot))

#map_rot = list(H1_plus_hres[pix_rot])

#hp.mollview(np.absolute(H1_plus_hres), title = "H1 Plus")
#plt.savefig('./H1_Plus_Rot.png')

plt.close('all')
