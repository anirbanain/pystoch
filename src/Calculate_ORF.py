import math
import cmath
import healpy as hp
import numpy as np
import pycbc.detector
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Parameters for the High Resolution Map
nside_hres = 256
npix_hres = hp.nside2npix(nside_hres)
pix_hres = np.arange(npix_hres)

(theta_hres, phi_hres) = \
hp.pix2ang(nside_hres, np.arange(hp.nside2npix(nside_hres)))

dec_hres = (np.pi/2) - theta_hres
ra_hres = phi_hres

# The two LIGO detectors
H1 = pycbc.detector.Detector('H1')
L1 = pycbc.detector.Detector('L1')

# Plus and Cross Antenna Pattern for the detectors
(H1_plus, H1_cross) = np.vectorize(H1.antenna_pattern)(ra, dec, 0, 1126626073);
(L1_plus, L1_cross) = np.vectorize(L1.antenna_pattern)(ra, dec, 0, 1126626073);

plus = H1_plus * L1_plus
cross = H1_cross * L1_cross
