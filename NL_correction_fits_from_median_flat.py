from multiprocessing import Pool
from multiprocessing import cpu_count
from lsst.cp.pipe.utils import (funcPolynomial, irlsFit)
import numpy as np
import astropy.io.fits as pf
from matplotlib.backends.backend_pdf import PdfPages
import pylab as plt
import matplotlib.font_manager as fm
import os
import sigma_clip

def plotCoeffsMatrix (coeffs, pdfPages, title=''):
    """Make a histogram and a 2D plot of NL-correction coefficients
    
    Parameters
    ----------
    coeffs : `numpya.array`
        NL correction coefficients
    """
    fig = plt.figure()
    sigmaCut = 3
    meanCoeffs, scatterCoeffs, mask = sigma_clip.sigma_clip(coeffs.flatten(), niter=10, 
                                                   nsig=sigmaCut, get_indices=True)
    print ("Mean and std: ", meanCoeffs, scatterCoeffs)
    print ("Median: ", np.median(coeffs))
    prop = fm.FontProperties(size=7)
    loc_label = 'upper right'
    
    ax=fig.add_subplot(111)
    plt.imshow(coeffs, cmap='viridis', interpolation='nearest', origin='lower',
               vmin=meanCoeffs - sigmaCut*meanCoeffs,
               vmax=meanCoeffs + sigmaCut*meanCoeffs)
    plt.colorbar()
    ax.set_title (r"correction coeff", size=11)
    fig.suptitle(title)
    pdfPages.savefig()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bins, patches_out = ax.hist(coeffs.flatten()[mask], 50, facecolor='red',
            alpha=0.75, label=f"Mean: {meanCoeffs:.4e} \n Scatter: {scatterCoeffs:.4e}")
    ax.set_title('Histogram after %g-sigma clipping' %sigmaCut, size=10)
    ax.legend(loc=loc_label , fancybox=True, ncol=1, numpoints=1, prop = prop)
    ax.tick_params(axis='both', which='major', labelsize=11.5)
    fig.suptitle(title)
    #plt.tight_layout()
    pdfPages.savefig()


def fit_pixel_ramp(medianFlatsCube, expTimes, counter=0, detSize=(4096, 4096)):
    index_x, index_y = np.unravel_index(counter, detSize)
    polyFit, polyFitErr, chiSq, weights = irlsFit([0.0, 100.0, -1e-6, 1e-10], expTimes,
                                                  medianFlatsCube[:, index_x, index_y],
                                                  funcPolynomial)
    # Correction coefficients (Eq. 37 of Jenna's paper, 3rd BFE paper from OSU group)
    k1 = polyFit[1]
    correctionCoeffs = [-coeff/(k1**order) for order, coeff in enumerate(polyFit)]
    return correctionCoeffs[2:] #(polyFit, polyFitErr, chiSq, weights)

def nl_function (index):
    """Function to be passed to nultiprocess.
    Must have single index as input.
    """
    index_x, index_y = np.unravel_index(index, (detSizeX, detSizeY))
    corrCoeffs = fit_pixel_ramp (medianFlatsCube, expTimes, counter=index, detSize=(detSizeX, detSizeY))
    return index_x, index_y, corrCoeffs[0], corrCoeffs[1] #polyFit, polyFitErr, chiSq, weights)

c2MatrixFile = "./NL_C2_Coeffs.dat"
c3MatrixFile = "./NL_C3_Coeffs.dat"
if not (os.path.exists(c2MatrixFile) and os.path.exists(c3MatrixFile)):
    #If coeff. matrices don't exist in disk, run the fits.
    # Read in median (of 48 exposures) flat and dark
    hdulist = pf.open("/project/plazas/PPL/H4RG/output/2021MAR26/PPL-data-2021-02-18/stacked/flats_median_stacked.fits")
    medianFlatsCube = 2**16 - 1 - hdulist[0].data
    headerFlatsCube = hdulist[0].header

    hdulist = pf.open("/project/plazas/PPL/H4RG/output/2021MAR26/PPL-data-2021-02-18/stacked/darks_median_stacked.fits")
    medianDarksCube = 2**16 - 1 - hdulist[0].data
    headerDarksCube = hdulist[0].header

    framtime = headerFlatsCube['FRAMTIME']
    nframes = headerFlatsCube['NFRAMES']
    expTimes = np.array([i*framtime for i in range(nframes)])/1000.
    detSizeX, detSizeY = headerFlatsCube['DETSIZEX'], headerFlatsCube['DETSIZEY']

    processes = cpu_count()
    use = 20
    print("I have ", processes, "cores here. Using: %g" % use)
    pool = Pool(processes=use)
    npix_total = detSizeX*detSizeY
    results = pool.map(nl_function, list(range(npix_total)))
    pool.close()
    pool.join()

    results = np.array(results)

    c2Matrix = results[:,2].reshape((detSizeX,detSizeY))
    c3Matrix = results[:,3].reshape((detSizeX,detSizeY))

    np.savetxt(c2MatrixFile, c2Matrix)
    np.savetxt(c3MatrixFile, c3Matrix)
else:
    # Don't do 16 million fits if the matrices are already in disk.
    print ("Reading files from disk")
    c2Matrix = np.genfromtxt (c2MatrixFile)
    c3Matrix = np.genfromtxt (c3MatrixFile)

pp = PdfPages("outNLCoeffs.pdf")
plotCoeffsMatrix (c2Matrix, pp, title="C2 coefficient (NL corr. polynomial, order=3)")
plotCoeffsMatrix (c3Matrix, pp, title="C3 coefficient (NL corr. polynomial, order=3)")
pp.close()

