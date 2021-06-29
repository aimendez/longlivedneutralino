######################################################################
'''
Aeff_decay Integral in [arXiv:1910.12839v2]
@Authors:
    - Pablo Candia (pcandia@uc.cl)
    - Giovanna Cottin ( @uai.cl)
    - Andres Mendez (aimendez@uc.cl)
    - Victor Munoz (vmmunoz2@uc.cl)
'''
######################################################################

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import RectBivariateSpline
from scipy import integrate
from math import pi as PI
import os

# Definitions in Appendix A of [1910.12839]

def delta_l1(H, R, theta, phi, r):
    return np.minimum(H/(np.abs(np.cos(theta))),(R*np.sqrt(1. - r**2*np.sin(phi)**2/R**2) + r*np.cos(phi))/(np.abs(np.sin(theta))))

def delta_l2(H, R, theta, phi, x):
    return np.minimum((H-x)/(np.abs(np.cos(theta))),2*R*np.cos(phi)/np.abs(np.sin(theta)))

def integrand_top(H, R, theta, phi, r, ctau):
    return r*(1. - np.exp(-delta_l1(H, R, theta, phi, r)/ctau))

def integrand_side(H, R, theta, phi, x, ctau):
    return R*(1. - np.exp(-delta_l2(H, R, theta, phi, x)/ctau))

def integral_top(H, R, theta, ctau):
    return integrate.nquad(lambda r, phi: integrand_top(H, R, theta, phi, r, ctau),[[0,R],[0,2*PI]])

def integral_side(H, R, theta, ctau):
    return integrate.nquad(lambda x, phi: integrand_side(H, R, theta, phi, x, ctau),[[0,H],[-PI/2.,PI/2.]])


def GetAeffDecay(detector):
	if detector == 'SK':
		H = 0.04
		R = 0.02
	elif detector == 'ANTARES':
		H = 0.48
		R = 0.09
	
	#check if file already exist:
	if os.path.exists(f"./datafiles/Aeff/{detector}_Aeff.txt") :
		return None

	ctau_list = np.logspace(-8,17,800)
	cos_theta_list = np.linspace(-0.9,0.9,10)
	A_eff_list = []
	for cos_theta in cos_theta_list:
	    theta = np.arccos(cos_theta)
	    sin_theta = np.sin(theta)
	    for ctau in ctau_list:
	        A_top = np.abs(cos_theta)*integral_top(H, R, theta, ctau)[0]
	        A_side = np.abs(sin_theta)*integral_side(H, R, theta, ctau)[0]
	        A_eff_list.append([A_top, A_side, cos_theta, ctau, (A_top + A_side)*1e10]) # Save total area in cm^2
	
	return A_eff_list

if __name__ == '__main__':
	detector = 'SK'
	A_eff_list = GetAeffDecay(detector)
	if A_eff_list is not None:
		np.savetxt(f"./datafiles/Aeff/{detector}_Aeff.txt", A_eff_list)