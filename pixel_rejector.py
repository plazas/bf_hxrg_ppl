#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.stats import norm
from scipy import ndimage


def estimateSigma(data=None, mu=0., maxOrder=None, minOrder=None, n_iter=1000):
    factor = 0.
    for i in range(n_iter):
        sample = np.random.randn(data.size) + mu
        factor = factor + np.std(sample[minOrder:maxOrder])
    factor = factor*1./n_iter
    return factor * np.std(data[minOrder:maxOrder])


def getBounds(data=None, center=None, thresh=60, median=False):
    from scipy.special import erf, erfc
    from scipy.optimize import minimize
    # First, get mode.
    if center == None:
        lo, hi = np.percentile(data, [10, 90])
        bins = np.linspace(lo, hi, 151)
        h, _ = np.histogram(data, bins)
        ff = np.zeros(21)
        ff[5:16] = 1./11.
        yh = np.convolve(ff, h, mode='same')
        maxind = np.argmax(yh)
        mode = np.mean([bins[maxind], bins[maxind+1]])
    elif median == False:
        mode = center
    else:
        print("using the median to center the distribution")
        mode = np.median(data)
    print("assuming mode is:", mode)
    # Now assume that the data are gaussian around this point.
    dist = data - mode
    outer_bound_upper = np.percentile(np.abs(dist[dist > 0]), thresh)
    outer_bound_lower = np.percentile(np.abs(dist[dist < 0]), thresh)
    print("upper, lower bounds:", outer_bound_upper, outer_bound_lower)

    window_bound = np.min([outer_bound_upper, outer_bound_lower])
    these = ((data - mode) > -window_bound) & ((data-mode) < window_bound)
    initial_sigma = (outer_bound_upper + outer_bound_lower)/2.

    def gaussian_logL(sigma):
        det = np.sqrt(2*np.pi*sigma**2)
        norm = 2./erf(window_bound/np.abs(sigma)/np.sqrt(2))
        #logL_each = np.log(np.exp(-(data[these] - mode)**2/sigma**2/2.)/det * norm)
        logL_each = -(data[these]-mode)**2/sigma**2 / \
            2. - np.log(det) + np.log(norm)
        return -np.mean(logL_each)
    res = minimize(gaussian_logL, initial_sigma)
    sigma_est = np.abs(res.x[0])
    ncore = 1./erf(window_bound/np.abs(sigma_est)/np.sqrt(2)) * np.sum(these)

    return sigma_est, ncore, mode, mode + outer_bound_upper, mode - outer_bound_lower


def evaluate_logL(data, sigma=None, mu=None):
    det = np.sqrt(2*np.pi*sigma**2)
    logL_each = -(data - mu)**2/sigma**2/2. - np.log(det)
    return logL_each


def evaluate_pdf(data, sigma=None, mu=None):
    det = np.sqrt(2*np.pi*sigma**2)
    logL_each = -(data - mu)**2/sigma**2/2. - np.log(det)
    L = np.exp(logL_each)
    return L


def get_test_data(mu=None, sigma=None):
    data1 = sigma * np.random.randn(197000) + mu
    data2 = .3*np.random.randn(3000) - 3.
    data = np.concatenate((data1, data2))
    return data


def get_real_data(dfile='./2017-01-13/master-dark.fits'):
    data = fits.getdata(dfile).astype('float')
    #data = data[data != 0]
    return data


def determine_quality_threshold(objval, bins, qhist, nkeep, forgive=None):
    # return array of indices into original object array of those that pass quality cuts.
    bin_inds = np.digitize(objval, bins)-1

    objQual = np.zeros_like(objval)
    for ind in np.unique(bin_inds[(bin_inds > -1) & (bin_inds < qhist.size)]):
        objQual[bin_inds == ind] = qhist[ind]
    if (np.min(bin_inds) == -1) | (np.max(bin_inds) >= qhist.size):
        objQual[(bin_inds == -1) | (bin_inds >= qhist.size)] = np.max(qhist)
    keep_frac = nkeep*100./objval.size
    thresh = np.percentile(objQual, keep_frac)
    if forgive is not None:
        thresh = thresh*(1+forgive)
    bad_inds = objQual > thresh
    return thresh, bad_inds


def find_anomalies(data, center=None, forgive=.25):
    sigma_est, ncore, mode, upper_bd, lower_bd = getBounds(
        data=data, thresh=80, median=False, center=center)
    print("best sigma guess for sigma is:", sigma_est)
    print("number of core elements is:", ncore)
    dist = norm(loc=mode, scale=sigma_est)
    low_bd, hi_bd = dist.interval(1 - 1./data.size)
    bins = np.linspace(low_bd, hi_bd, 150)
    h, _ = np.histogram(data, bins=bins, density=True)
    bin_centers = (bins[0:-1] + bins[1:])/2.
    L = evaluate_pdf(bins, sigma=sigma_est, mu=mode)
    Lmean = (L[0:-1] + L[1:])/2. * data.size
    excess = np.abs(h*data.size - Lmean)
    excess_err = np.sqrt(Lmean)
    exthresh, mask = determine_quality_threshold(
        data, bins, excess/excess_err, ncore, forgive=forgive)
    mask_open = ndimage.binary_opening(mask)
    mask_adj = ndimage.binary_closing(mask_open)

    return mask, mask_adj


def main(argv):
    sigma_true = 1.0
    mu_true = 0.
    forgive = 0.25
    data = get_test_data(mu=mu_true, sigma=sigma_true)
    #data = get_real_data()
    sigma_est, ncore, mode, upper_bd, lower_bd = getBounds(
        data=data, thresh=80, median=False)
    logL = evaluate_logL(data, sigma=sigma_est, mu=mode)
    # print "bound is:",bd
    print("best sigma guess, sigma actual is:", sigma_est, sigma_true)
    print("number of core elements is:", ncore)

    # Now that we have the model fit
    # and the log-likelihoods, figure
    # out which points don't belong.
    # Which points have a low likelihood,
    # given all of the others?

    from scipy.stats import norm
    dist = norm(loc=mode, scale=sigma_est)
    # thresh = 1./npts
    low_bd, hi_bd = dist.interval(1 - 1./data.size)
    print(low_bd, hi_bd)
    bins = np.linspace(low_bd, hi_bd, 150)
    h, _ = np.histogram(data, bins=bins, density=True)
    bin_centers = (bins[0:-1] + bins[1:])/2.
    L = evaluate_pdf(bins, sigma=sigma_est, mu=mode)
    Lmean = (L[0:-1] + L[1:])/2. * data.size
    plt.plot(bin_centers, h*data.size)
    plt.plot(bin_centers, Lmean)
    plt.axvline(upper_bd, color='red', linestyle='--')
    plt.axvline(lower_bd, color='red', linestyle='--')
    plt.yscale('log')
    plt.show()

    # To establish the thresh, let's figure out how many objects we'll
    # actually want to throw out.
    ngood = ncore*1./data.size
    # Let's assume that the Poisson error is 1./sqrt(n),
    # for n in the histogram bins.
    excess = np.abs(h*data.size - Lmean)
    excess_err = np.sqrt(Lmean)
    exthresh, bad_inds = determine_quality_threshold(
        data, bins, excess/excess_err, ncore)
    exthresh = exthresh*(1+forgive)

    mask = np.zeros_like(data).astype('int')
    mask[bad_inds] = 1
    # Let's try simplifying the mask.
    from scipy import ndimage
    mask_open = ndimage.binary_opening(bad_inds)
    mask_close = ndimage.binary_closing(mask_open)
    mask2 = np.zeros_like(data).astype('int')
    mask2[mask_close] = 1
    #mask[bad_inds & (data < mode)] = 2
    # fits.writeto('dark_mask.fits',mask,clobber=True)
    # fits.writeto('dark_mask2.fits',mask2,clobber=True)
    plt.plot(bin_centers, excess / excess_err)
    plt.axhline(exthresh, color='red', linestyle='--')
    plt.show()


if __name__ == "__main__":
    import pdb
    import traceback
    import sys
    try:
        main(sys.argv)
    except:
        thingtype, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
