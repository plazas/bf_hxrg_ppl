import numpy as np
import math

# Formula from Chris hirata
def calculate_B_gaussian(sigma,aNN, aD):
    const = 1./(2*np.sqrt(2)*sigma)
    num = 2*np.sqrt(2)*np.exp(-1./(8*sigma**2))
    den = np.sqrt(np.pi)*sigma*math.erf(const)
    factor = (aNN + 0.5*(1 + (math.erf(3*const))/(math.erf(const)))*aD)

    return -(num/den)*factor


# Approximate an Airy with a gaussian of sigma ~ 0.42\Lambda*f-number

# Fig. 5, Choi and Hirata
ppm = 1.e-6
aNN = 0.29*ppm
aD = 0.038*ppm
f_number = 8
sigma_PSF_PPL = 0.42*1000*f_number


B = calculate_B_gaussian(sigma_PSF_PPL, aNN, aD)

print ("B", B)


B0 = calculate_B_gaussian(sigma_PSF_PPL, aNN, 0.0)

print ("B0 (no diagonal terms)", B0)

print (-4*aNN)
