import numpy as np
import astropy.io.fits as pf
import timeit

start = timeit.timeit()

print ("Reading files")
spotsFile = "/projector/aplazas/stacked/2021MAR29/spots_median_stacked.fits"
hdulist = pf.open(spotsFile)
medianSpots = hdulist[0].data
medianSpotsHeader = hdulist[0].header

flatsFile = "/projector/aplazas/stacked/2021MAR29/flats_median_stacked.fits"
hdulist = pf.open(flatsFile)
medianFlats = hdulist[0].data
medianFlatsHeader = hdulist[0].header

darksFile = "/projector/aplazas/stacked/2021MAR29/darks_median_stacked.fits"
hdulist = pf.open(darksFile)
medianDarks = hdulist[0].data
medianDarksHeader = hdulist[0].header

end = timeit.timeit()
print(end - start)



start = timeit.timeit()
print ("Reading coeff. matrices")
c2Matrix = np.genfromtxt("/projector/aplazas/NL_MATRICES/2021MAR29/NL_C2_Coeffs.dat")
c3Matrix = np.genfromtxt("/projector/aplazas/NL_MATRICES/2021MAR29/NL_C3_Coeffs.dat")
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

print ("Writting corrected stacked files to disk")

pf.writeto ("/projector/aplazas/stacked/2021MAR29/spots_median_stacked_NL_corrected.fits",
            correctedSpots, header=medianSpotsHeader)

pf.writeto ("/projector/aplazas/stacked/2021MAR29/darks_median_stacked_NL_corrected.fits",
            correctedDarks, header=medianDarksHeader)

pf.writeto ("/projector/aplazas/stacked/2021MAR29/flats_median_stacked_NL_corrected.fits",
            correctedFlats, header=medianFlatsHeader)




