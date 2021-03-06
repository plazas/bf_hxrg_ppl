import numpy as np
import astropy.io.fits as pf
import timeit

start = timeit.timeit()

print ("Reading files")
spotsFile = "/projector/aplazas/stacked/2021APR14/2021APR14-median-stacked-spots-PPL-2021-02-18.fits"
#spotsFile = "/projector/aplazas/stacked/2021MAR29/spots_median_stacked.fits"
hdulist = pf.open(spotsFile)
medianSpots = 2**16 - 1 - hdulist[0].data
medianSpotsHeader = hdulist[0].header

#flatsFile = "/projector/aplazas/stacked/2021MAR29/flats_median_stacked.fits"
flatsFile = "/projector/aplazas/stacked/2021APR14/2021APR14-median-stacked-flats-PPL-2021-02-18.fits"

hdulist = pf.open(flatsFile)
medianFlats = 2**16 - 1 - hdulist[0].data
medianFlatsHeader = hdulist[0].header

#darksFile = "/projector/aplazas/stacked/2021MAR29/darks_median_stacked.fits"
darksFile = "/projector/aplazas/stacked/2021APR14/2021APR14-median-stacked-darks-PPL-2021-02-18.fits"
hdulist = pf.open(darksFile)
medianDarks = 2**16 - 1 - hdulist[0].data
medianDarksHeader = hdulist[0].header

end = timeit.timeit()
print(end - start)



start = timeit.timeit()
print ("Reading coeff. matrices")
#c2Matrix = np.genfromtxt("/projector/aplazas/NL_MATRICES/2021MAR29/NL_C2_Coeffs.dat")
#c3Matrix = np.genfromtxt("/projector/aplazas/NL_MATRICES/2021MAR29/NL_C3_Coeffs.dat")
#c2Matrix = np.genfromtxt("/projector/aplazas/NL_MATRICES/2021APR02/NL_C2_Coeffs_2021APR02.dat")
#c3Matrix = np.genfromtxt("/projector/aplazas/NL_MATRICES/2021APR02/NL_C3_Coeffs_2021APR02.dat")
c2Matrix = np.genfromtxt("/projector/aplazas/NL_MATRICES/2021APR14/NL_C2_Coeffs_2021APR14.dat")
c3Matrix = np.genfromtxt("/projector/aplazas/NL_MATRICES/2021APR14/NL_C3_Coeffs_2021APR14.dat")

end = timeit.timeit()
print(end - start)


def correctNL (medianImageCube):
    correctionCoeffs=[c2Matrix, c3Matrix]
    correction = np.zeros_like(medianImageCube)
    for order, coeff in enumerate(correctionCoeffs):
        correction+= coeff * np.power(medianImageCube, order+2)
    correctedImageCube = medianImageCube + correction
    return correctedImageCube

start = timeit.timeit()
print ("Performing correction")
correctedSpots = correctNL(medianSpots)
correctedDarks = correctNL(medianDarks)
correctedFlats = correctNL(medianFlats)
end = timeit.timeit()
print(end - start)


print ("Flats")
print (medianFlats[-1][500,500], correctedFlats[-1][500,500])
print (medianFlats[-1][2000,2000], correctedFlats[-1][2000,2000])

print ("Darks")
print (medianDarks[-1][500,500], correctedDarks[-1][500,500])
print (medianDarks[-1][2000,2000], correctedDarks[-1][2000,2000])

print ("Spots")
print (medianSpots[-1][500,500], correctedSpots[-1][500,500])
print (medianSpots[-1][2000,2000], correctedSpots[-1][2000,2000])

print ("Writting NL-corrected stacked files to disk")

correctedSpots = 2**16 - 1 - correctedSpots
correctedFlats = 2**16 - 1 - correctedFlats
correctedDarks = 2**16 - 1 - correctedDarks


outFileFlats = "/projector/aplazas/stacked/2021APR14/2021APR14-median-stacked-flats-PPL-2021-02-18-NL-corrected.fits"
outFileSpots =  "/projector/aplazas/stacked/2021APR14/2021APR14-median-stacked-spots-PPL-2021-02-18-NL-corrected.fits"
outFileDarks =  "/projector/aplazas/stacked/2021APR14/2021APR14-median-stacked-darks-PPL-2021-02-18-NL-corrected.fits"

pf.writeto (outFileSpots,
            correctedSpots, header=medianSpotsHeader, overwrite=True)

pf.writeto (outFileDarks,
            correctedDarks, header=medianDarksHeader, overwrite=True)

pf.writeto (outFileFlats,
            correctedFlats, header=medianFlatsHeader, overwrite=True)




