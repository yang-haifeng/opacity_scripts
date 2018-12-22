#!/usr/bin/env python

from bhmie import *
import numpy as np
from scipy.interpolate import interp1d

def compute_opac_mie(refidx, agarr, wgtarr, lam, theta=None, errtol=0.01, verbose=False):
    """
    Calculate the opacity for a given refractive index and a grain size distribution.
    refidx: complex refractive index m = n + i k
    agarr : array of grain sizes
    wgtarr: array of weights for averaging. n(a)*m(a)*da
    """
    if theta is None:
        angles = np.array([0.,90.,180.])
    else:
        angles = theta
    nang = len(angles)
    nagr = len(agarr)
    kabs  = 0.0
    kscat = 0.0
    gscat = 0.0
    siggeo = np.pi*agarr**2
    mgrain = 4*np.pi/3.0 * dens * agarr**3

    for i in range(nagr):
        x = 2*np.pi*agarr[i]/lam
        S1, S2, Qext, Qabs, Qsca, Qback, gsca = bhmie(x, refidx, angles)

        kabs  += wgt[i] * Qabs*siggeo[i]
        kscat += wgt[i] * Qabs*siggeo[i]
        gscat += wgt[i] * gsca
