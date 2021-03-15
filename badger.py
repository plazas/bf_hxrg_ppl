#!/usr/bin/env python
import os
import numpy as np
from astropy.io import fits
import sextractor_engine as seng
import pixel_rejector as pr
import glob


def quad_correct(pixel_series=None, bias=None, flat=None, ivar_series=None, time_series=None, nl_coeff=None):
    # Compute the polynomial coefficients for this pixel.
    # This is only a second-order nonlinearity correction,
    #   so nl_coeff should be a scalar.
    if ivar_series == None:
        ivar_series = pixel_series*0.+1.
    if flat == None:
        flat = 1.
    t = np.array([x*flat for x in time_series])
    a = 2 * nl_coeff*nl_coeff*np.sum(t*t*t*t * ivar_series)
    b = 3 * nl_coeff * np.sum(t*t*t*ivar_series)
    c = np.sum((1 - 2 * nl_coeff * (pixel_series - bias))
               * t * t * ivar_series)
    d = -np.sum((pixel_series - bias) * t * ivar_series)
    roots = np.roots([a, b, c, d])
    # There should be multiple roots, unless we're very lucky.
    # We need to choose the root that is: a) Real, and b) maximizes the likelihood.
    # It's theoretically possible that there is only one! But in this case, it will be repeated.
    real_roots, nreal = np.unique(
        roots[np.imag(roots) == 0], return_counts=True)
    if np.sum(nreal) == 1:
        # if there is only only one real root, we've found our answer.
        fest = roots[np.imag(roots) == 0]
        logL = -np.sum(ivar_series * ((pixel_series - (bias + fest * t)))**2)
    if len(nreal) > 1:
        # if there's more than one unique real root:
        logL = np.zeros(len(nreal))
        for i, thisroot in enumerate(real_roots):
            logL[i] = -np.sum(ivar_series * ((pixel_series - (bias +
                                                              thisroot*t +
                                                              nl_coeff*(thisroot*t)**2)))**2)
        # which real root maximizes L?
        fest = real_roots[logL == np.max(logL)]
    return np.real(fest), np.max(logL)


def alpha_estimator(flux=None, pixel_series=None, bias=None, flat=None, ivar_series=None, time_series=None, diagnostics=False):
    # Given knowledge of the flux, estimate the nonlinearity.
    if ivar_series == None:
        ivar_series = pixel_series*0.+1.
    if flat == None:
        flat = 1.
    if bias == None:
        bias = pixel_series[0]
    if time_series == None:
        time_series = np.array(len(pixel_series))
    t = time_series * flat
    a = (flux*t)**2
    b = (pixel_series - bias - flux*t)
    d = np.sum(a**2 * ivar_series)
    est = np.sum(a * ivar_series * b) * 1. / d
    if diagnostics == False:
        return est
    else:
        return t, a, b, d


def crudeLinearFlux(pixel=None, timeseries=None, polynomial_order=2):
    # Fit a crude polynomial of some order to the pixel time series.
    coeff = np.polyfit(timeseries, pixel, polynomial_order)
    linear = coeff[-2]
    return linear


def makeNonlinearityTable(ramp=None, ivar=None, timeseries=None):
    # "ramp" and "ivar" are 3d numpy arrays arranged as [time, npix1, npix2]
    # "timeseries" is a list or array of times which will be used to get the flux.
    #   If this is not given, we just assume the ramp is correctly ordered; then the flux units
    #   will be data units per sample, rather than data units per time.
    # --------------------------------------------------
    # We will assume that the first entry in "ramp" is the bias frame.
    if ivar == None:
        ivar = ramp*0. + 1.
    bias = ramp[0, :, :]
    nframes = ramp.shape[0]
    if timeseries == None:
        timeseries = np.arange(nframes)
    else:
        # the estimator assumes that ideally, data = f*time + bias.
        # We'll do the dark-current correction downstream, since really data = (f+dc)*time+bias
        # so time needs to start at zero.
        timeseries = timeseries - timeseries[0]
    t = timeseries
    # make a 3d array of times.
    tgrid = np.zeros(ramp.shape)
    # make a 3d array of biases.
    bgrid = np.zeros(ramp.shape)
    for i, t in enumerate(timeseries):
        tgrid[i, :, :] = timeseries[i]
        bgrid[i, :, :] = bias
    else:
        alpha = np.zeros((ramp.shape[1], ramp.shape[2]))
        for i in range(ramp.shape[1]):
            for j in range(ramp.shape[2]):
                # First figure out what we think the linear flux is in each pixel.
                flux = crudeLinearFlux(
                    pixel=ramp[:, i, j], timeseries=timeseries, polynomial_order=2)
                # alpha[i,j] = alpha_estimator(flux = flux,
                #                             pixel_series=ramp[:,i,j],
                #                             bias = bias[i,j],
                #                             ivar_series = ivar[:,i,j],
                #                             flat=1.,
                #                             time_series=timeseries)
                alpha[i, j] = np.polyfit(timeseries*flux, ramp[:, i, j], 2)[0]
    return alpha


def deriveNonlinearity(flatFiles=None, samples_per_ramp=None, ramps_per_file=None, times=None, time_axis=0, doplot=True):
    # The "flatFiles", "samples_per_ramp", "ramps_per_file", and "times" arguments
    #  here all need to be lists, and they all need to be the same length.
    nFiles = len(flatFiles)
    nl_coeff_est = []
    if doplot == True:
        import matplotlib.pyplot as plt
    index = 0
    for i, imFile in enumerate(flatFiles):
        image = fits.getdata(imFile).astype('float')
        flatRamps = makeRampFromStrip(image_orig=image,
                                      n_sample_per_ramp=samples_per_ramp[i],
                                      n_ramp_per_file=ramps_per_file[i],
                                      time_axis=time_axis[i])
        for j, fRamp in enumerate(flatRamps):
            this_nl_coeff = makeNonlinearityTable(
                ramp=fRamp, timeseries=times[i])
            nl_coeff_est.append(this_nl_coeff)
            if doplot == True:
                plt.imshow(this_nl_coeff, vmin=-3e-6,
                           vmax=3e-6, origin='lower')
                plt.colorbar()
                plt.savefig('nl_coeff-'+str(index)+'.pdf', format='pdf')
                plt.clf()
            index = index+1
    # stack the nonlinearity coefficients together
    nl_coeff_final = np.sum(np.stack(nl_coeff_est), axis=0) * \
        (1./(ramps_per_file[0] * nFiles))
    return nl_coeff_final


def reference_subtraction(image, mask, n_channel=32):
    maskdef = define_mask_bits()

    pix_per_channel = image.shape[0] / n_channel

    col = np.arange(mask.shape[0])
    reference_frame = np.zeros_like(image)
    #reference_frame = image * 0.
    channel = np.outer(np.ones(image.shape[0]), ((
        col - (col % pix_per_channel)) / pix_per_channel))

    for i in range(n_channel):
        reference_pixels = (
            (mask & maskdef['reference']) == maskdef['reference']) & (channel == i)
        reference_frame[channel == i] = np.median(image[reference_pixels])

    return image - reference_frame


def makeFrameFromRamp(ramp=None, flat=None, ivar=None, exptime=None, nl_coeff=None, getLogL=False, ramp_start=1,
                      dark=None, mask=None, roi=None, ref_subtract=True):
    # "ramp" and "ivar" are 3d numpy arrays arranged as [time, npix1,  npix2]
    # "timeseries" is a list or array of times which will be used to get the flux.
    #   If this is not given, we just assume the ramp is correctly ordered; then the flux units
    #   will be data units per sample, rather than data units per time.
    # --------------------------------------------------
    nframes = ramp.shape[0]
    timeseries = np.linspace(0, exptime, nframes)

    # If no mask provided, just use the default.
    if mask == None:
        mask = baseMask()
    # Do the appropriate reference image subtraction
    if ref_subtract is True:
        ramp_ref = ramp.copy()
        for i in range(nframes):
            ramp_ref[i, :, :] = reference_subtraction(ramp[i, :, :], mask=mask)
        ramp = ramp_ref
    # We will assume that the entry referenced by "ramp_start" is the bias frame, and ignore everything prior.
    # But we will require that there be at least two usable frames; ramp start is set to the provided value, or to len(ramp) -2
    if ramp.shape[0] - ramp_start < 2:
        ramp_start = ramp.shape[0] - 2
    if flat == None:
        flat = ramp[0, :, :] * 0. + 1.
    if ivar == None:
        ivar = ramp*0. + 1.
    ivar[0:ramp_start, :, :] = 0
    if dark == None:
        dark = ramp[0, :, :]*0.

    bias = ramp[ramp_start, :, :]
    timeseries = timeseries - timeseries[ramp_start]
    # the estimator assumes that ideally, data = f*time + bias.
    # We'll do the dark-current correction downstream, since really data = (f+dc)*time+bias
    # so time needs to start at zero.
    timeseries = timeseries - timeseries[0]
    t = timeseries
    # make a 3d array of times.
    tgrid = np.zeros(ramp.shape)
    bgrid = np.zeros(ramp.shape)
    fgrid = np.zeros(ramp.shape)
    dgrid = np.zeros(ramp.shape)
    logL = np.zeros((ramp.shape[1], ramp.shape[2]))
    for i, t in enumerate(timeseries):
        tgrid[i, :, :] = timeseries[i]
        fgrid[i, :, :] = flat
        bgrid[i, :, :] = bias
        dgrid[i, :, :] = timeseries[i] * dark
    if nl_coeff == None:
        f = np.sum((ramp - bias - dgrid)*tgrid * ivar, axis=0) / \
            (np.sum((flat*tgrid)**2 * ivar, axis=0))
        #logL = -np.sum( (ramp - (f *tgrid + bias ) )**2 * ivar, axis=0)
    else:
        f = np.zeros((ramp.shape[1], ramp.shape[2]))
        logL = np.zeros((ramp.shape[1], ramp.shape[2]))
        for i in range(ramp.shape[1]):
            for j in range(ramp.shape[2]):
                thisf, thisl = quad_correct(pixel_series=ramp[:, i, j],
                                            bias=bias[i, j],
                                            ivar_series=ivar[:, i, j],
                                            flat=flat[i, j],
                                            time_series=timeseries,
                                            nl_coeff=nl_coeff[i, j])
                f[i, j] = thisf
                logL[i, j] = np.real(thisl)
        #logL = -np.sum( (ramp - (bias + f*tgrid + nl_coeff * (f*tgrid) * (f*tgrid) ) )* ivar ,axis=0)
    # if the flat field response is zero, the inferred flux will be infinite. For now, just set it to zero.
    # Later, we'll add these pixels to a bad pixel mask.

    # if ref_subtract is True:
    #    f = reference_subtraction(f, mask = mask)
    maskdef = define_mask_bits()

    f[flat == 0] = 0.
    #f[~(mask & maskdef['good'] == maskdef['good'])] = 0.
    if getLogL == True:
        return f, logL
    else:
        return f


def organize_ramps(files=None,  expr='([0-9]){4}_([0-9]){4}'):
    # Determine the number of ramps.
    # Assign a ramp membership to each file.
    # Return this structure for processing.
    nFiles = len(files)
    fileStr = np.empty(nFiles, dtype=[('fileName', object), ('baseName', object), ('mode', object), (
        'ramp', np.int), ('sample', np.int), ('EXPTIME', np.float), ('FRAMTIME', np.float), ('LODFILE', object)])
    for i, iFile in enumerate(files):
        print("iFile: ", iFile)
        thisHdr = fits.getheader(iFile)
        fileStr['fileName'][i] = iFile
        fileStr['baseName'][i] = os.path.basename(thisHdr['BASE'])
        fileStr['ramp'][i] = thisHdr['IMAGEN']
        fileStr['EXPTIME'][i] = thisHdr['EXPTIME']
        fileStr['FRAMTIME'][i] = thisHdr['FRAMTIME']
        fileStr['LODFILE'][i] = thisHdr['LODFILE']
        if 'True' in thisHdr['FILMSTRP']:
        #if True:
            fileStr['mode'][i] = 'Filmstrip'
            fileStr['sample'][i] = thisHdr['RAMPLEN']
        else:
            fileStr['mode'][i] = 'Multiple'
            fileStr['sample'][i] = thisHdr['RAMPID']
    return fileStr


def getRampsFromFiles(files=None):
    # Get an index of the files.
    index = organize_ramps(files=files)
    rampNumbers = np.unique(index['ramp'])
    allRamps = []
    finalIndex = []
    for thisRamp in rampNumbers:
        these = index[(index['ramp'] == thisRamp)]
        these = these[np.argsort(these['sample'])]
        finalIndex.append(these[-1])
        # read in each sample associated with this ramp.
        thisRamp = []
        for this in these:
            if 'Multiple' in this['mode']:
                thisRamp.append(fits.getdata(this['fileName']).astype('float'))
            elif 'Filmstrip' in this['mode']:
                image = fits.getdata(this['fileName']).astype('float')
                thisRamp = makeRampFromStrip(image, n_sample_per_ramp=this['sample'])[
                    0]  # CDS only
        if 'Multiple' in these['mode']:
            thisRamp = np.stack(thisRamp)
        allRamps.append(thisRamp)
    finalIndex = np.stack(finalIndex)
    return allRamps, finalIndex


def makeRampFromStrip(image_orig=None, n_sample_per_ramp=None, n_ramp_per_file=1, time_axis=0):
    # Chop up a filmstrip (consisting of two or more concatenated samples):
    # How long is each thing?
    # Check that the long axis is actually divisible by n_sample.
    if time_axis == None:
        time_axis = 0
    if not (image_orig.shape[time_axis] % (n_sample_per_ramp * n_ramp_per_file)) == 0:
        raise Exception("the number of samples along the film strip is wrong.")
    ramp_length = image_orig.shape[time_axis] / n_ramp_per_file
    frame_length = image_orig.shape[time_axis] / \
        n_sample_per_ramp / n_ramp_per_file

    if time_axis == 0:
        image = image_orig.copy()
    if time_axis == 1:
        image = image_orig.T.copy()
    if time_axis > 1:
        raise Exception(
            "this script expects the image to be two-dimensional, but you have supplied an image with shape "+str(image_orig.shape))
    ramp_list = []
    for i in range(n_ramp_per_file):
        xstart = i * ramp_length
        xend = (i+1) * ramp_length
        thisImage = image[xstart:xend, :]
        thisFrames = []
        for j in range(n_sample_per_ramp):
            thisFrames.append(thisImage[j*frame_length:(j+1)*frame_length, :])
        frame = np.array(thisFrames)
        ramp_list.append(frame)
    return ramp_list


def updateHeader(header, inputFiles=None, masterDark=None, masterFlat=None, imageType=None, pixScale=18.0, maskdist=274.5):
    # The new header needs to have a few things:
    # 1. a keyword telling us what kind of image this is.
    # 2. If it's a calibration image, we need to
    #    store the paths to the calibration images this was made from.
    # 3. If it's a science image, we need to store the paths
    #    to the master calibration images.
    header.set('IMGTYPE', imageType, before='COMMENT')
    if masterDark is not None:
        header.set('MSTRDARK', masterDark, before='COMMENT',
                   comment='master dark used for this image.')
    if masterFlat is not None:
        header.set('MSTRFLAT', masterFlat, before='COMMENT',
                   comment='master flat used for this image.')
    if inputFiles is not None:
        header.set('HISTORY', '        ')
        header.set(
            'HISTORY', 'This image was made from the following image files:')
        if type(inputFiles) == type(''):
            # This image was made from a single input image.
            header.set('HISTORY', inputFiles)
        elif type(inputFiles) == type([]):
            for thisFile in inputFiles:
                header.set('HISTORY', thisFile)
    # Also, some fields that can be added to everything.
    header.set('PIXSCALE', pixScale, comment='pixel scale', before='COMMENT')
    header.set('MASKDIST', maskdist,
               comment='spot spacing, where relevant', before='COMMENT')
    return header


def makeCatalog(imageFile=None, flagFileName=None, outfileName=None,
                sconfig=None, spars=None, nnw=None, sfilter=None):
    if outfileName == None:
        outFileTemplate, ext = os.path.splitext(imageFile)
        outFileName = outFileTemplate+'-catalog.fits'

    s_engine = seng.SextractorEngine(IMAGE=imageFile, CATALOG_NAME=outFileName,
                                     c=sconfig,
                                     PARAMETERS_NAME=spars,
                                     FILTER_NAME=sfilter,
                                     STARNNW_NAME=nnw)
    s_engine.config['FLAG_IMAGE'] = flagFileName
    s_engine.run()
    print("Writing catalog from "+imageFile+" to  "+outFileName)
    catalog = fits.getdata(outFileName, ext=2)
    return catalog, outfileName


def combine(ramps=None, exptimes=None, weights=None, ref_subtract=False, dark=None, ramp_start=1):
    # Just a wrapper for makeFrameFromRamp.
    # exptimes and ivars needs to be either a list of the same length as ramps, or else just a single float/array.
    allFrames = []
    if len(exptimes) != len(ramps):
        exptimes = np.zeros(len(ramps)) + exptimes
    for i in range(len(ramps)):
        allFrames.append(makeFrameFromRamp(
            ramps[i], exptime=exptimes[i], ref_subtract=ref_subtract, dark=dark, ramp_start=ramp_start))
    finalFrame = np.average(np.stack(allFrames), axis=0, weights=weights)
    return finalFrame


def define_mask_bits():
    maskdef = np.empty(1, dtype=[('good', np.uint8), ('boundary', np.uint8),
                                 ('reference', np.uint8), ('hot', np.uint8),
                                 ('dead', np.uint8), ('noisy', np.uint8),
                                 ('anomalous', np.uint8), ('coating_defect', np.uint8)])
    maskdef['good'] = 2**0
    maskdef['boundary'] = 2**1
    maskdef['reference'] = 2**2
    maskdef['hot'] = 2**3
    maskdef['dead'] = 2**4
    maskdef['noisy'] = 2**5
    maskdef['anomalous'] = 2**6
    maskdef['coating_defect'] = 2**7
    return maskdef


def baseMask(roi=None):
    maskdef = define_mask_bits()
    # everything within 4 pixels of a boundary is a boundary pixel.
    mask = np.zeros((2048, 2048), dtype=np.uint8)
    mask[0:4, :] = maskdef['boundary']
    mask[-4:, :] = maskdef['boundary']  # + maskdef['reference']
    mask[-4:, :] = mask[-4:, :] + maskdef['reference']
    mask[:, 0:4] = maskdef['boundary']
    mask[:, -4:] = maskdef['boundary']
    #mask[mask == 0] = maskdef['good']
    # The last four pixels in a column are reference pixels.
    return mask


def makeCorrectedFrame(files_science=None, files_dark=None, files_flat=None,
                       do_science=True, do_flat=True, do_dark=True,
                       master_dark=None, master_flat=None,
                       nl_correction=None, nl_correction_outFile=None,
                       bad_pixels=None, bad_pixel_outFile=None,
                       fluence=False, ramp_start=2,
                       masterDarkOutFile=None,
                       masterFlatOutFile=None,
                       outFileTemplate=None):
    '''
    This is the driver for the reduction scripts. It takes lists of
    science, dark, and flat files. It makes a single master flat, and a
    single master dark, and uses these for reduction of all of the
    provided science frames.

    If the outFileTemplate keyword set, each reduced science ramp will be written as a separate
    file; the file names will be outFileTemplate+'-'+[original file_name]+'-'+[number]+'.fits'

    the nl_coeff can be one of a few things:
    1. None -- in this case, no nonlinearity correction is performed.
    2. a filename (a string) -- in this case, we read a file, and treat it as tabulated
    nonlinearity coefficients (npix_x, npix_y, 3)
    3. a 3d numpy array --  in this case, we treat this as tabulated nonlinearity 
    coefficients (npix_x, npix_y, 3).
    4. True -- in this case, derive the NL coefficients from the supplied flat fields. 
    (This only works if flats are supplied).

    the bad pixel mask can be one of a few things:
    1. None -- in this case, no mask is applied.
    2. a filename (a string) -- in this case, we read a file, and treat it as a mask (npix_x, npix_y)
    3. a 3d numpy array --  in this case, we treat this as the mask (npix_x, npix_y).
    4. True -- in this case, derive the pixel mask from whatever calibration data is available
    '''
    # Nonlinearity correction must be done prior to flat field estimation.
    #
    if nl_correction is not None:
        if type(nl_correction) == type(''):
            nl_coeff = fits.getdata(nl_correction).astype('float')
        elif nl_correction == True:
            nl_coeff = deriveNonlinearity(flatFiles=files_flat,
                                          samples_per_ramp=nsample_flat,
                                          ramps_per_file=nramp_flat,
                                          times=time_flat,
                                          time_axis=time_axis_science)
            if nl_correction_outFile is not None:
                fits.writeto(nl_coeff_file, nl_coeff, clobber=True)
        else:
            nl_coeff = nl_correction
    else:
        nl_coeff = None

    # combine the darks and flats.
    # this should scale files with different cadences appropriately, so the output is a self-consistent estimate
    # of the signal rate. Note that it's not totally consistent to try to combine flats with different lamp powers.
    if do_dark == True:
        if (files_dark is not None) and (master_dark is None):
            dark_ramps, darkStr = getRampsFromFiles(files=files_dark)
            masterDark = combine(
                ramps=dark_ramps, exptimes=darkStr['EXPTIME'], ramp_start=ramp_start)
        elif master_dark is not None:
            masterDark = fits.getdata(master_dark).astype('float')
    else:
        masterDark = None

    if do_flat == True:
        if (files_flat is not None) & (master_flat is None):
            flat_ramps, flatStr = getRampsFromFiles(files=files_flat)
            masterFlat = combine(
                ramps=flat_ramps, exptimes=flatStr['EXPTIME'], dark=masterDark, ramp_start=ramp_start)
        elif master_flat is not None:
            masterFlat = fits.getdata(master_flat).astype('float')
    else:
        masterFlat = None

    if bad_pixels == None:
        mask = baseMask()
        maskdef = define_mask_bits()

        # if masterDark != None:
        #    pmask_dark,_ = pr.find_anomalies(masterDark)
        #    mask[pmask_dark] = mask[pmask_dark] + maskdef['hot']
    # --------------------------------------------------
    # Write the master dark and master flat to disk.
    if (masterDarkOutFile is not None) and (do_dark == True):
        if type(files_dark) == type([]):
            darkheader = fits.getheader(files_dark[0])
        else:
            darkheader = fits.getheader(files_dark)
        darkheader = updateHeader(
            darkheader, inputFiles=files_dark, imageType='MASTER DARK')
        fits.writeto(masterDarkOutFile, masterDark,
                     header=darkheader, clobber=True)
    else:
        print("not writing dark to disk")
    if (masterFlatOutFile is not None) and (do_flat == True):
        if type(files_flat) == type([]):
            flatheader = fits.getheader(files_flat[0])
        else:
            flatheader = fits.getheader(files_flat)
        flatheader = updateHeader(
            flatheader, inputFiles=files_flat, imageType='MASTER FLAT')
        fits.writeto(masterFlatOutFile, masterFlat,
                     header=flatheader, clobber=True)
    else:
        print("not writing flat to file")
    # --------------------------------------------------
    if do_science == False:
        print("Skipping science frame reduction.")
        return None, None, None

    correctedFrames = []
    outFileNames = []
    outFlagFileNames = []
    scienceRamps, sciStr = getRampsFromFiles(files_science)
    for i, sciRamp in enumerate(scienceRamps):
        correctedFrame = -makeFrameFromRamp(ramp=sciRamp, flat=masterFlat,
                                            exptime=sciStr['EXPTIME'][i], dark=masterDark,
                                            nl_coeff=nl_coeff, ref_subtract=False, ramp_start=ramp_start)

        if fluence is True:
            correctedFrame = correctedFrame * sciStr[i]['EXPTIME']
        med = np.median(correctedFrame[np.isfinite(correctedFrame)])
        correctedFrame[~np.isfinite(correctedFrame)] = 0.
        correctedFrame[(mask != 0) & (mask != maskdef['good'])] = 0.
        correctedFrames.append(correctedFrame)
        if outFileTemplate is not None:
            baseName = os.path.split(sciStr[i]['baseName'])[-1].rstrip('-_')
            outFileName = outFileTemplate+'-'+baseName + \
                '-'+str(sciStr[i]['ramp']).zfill(4)+'.fits'
            outFlagFileName = outFileTemplate+'-'+baseName + \
                '_flags-'+str(sciStr[i]['ramp']).zfill(4)+'.fits'
            sciheader = fits.getheader(sciStr[i]['fileName'])
            sciheader = updateHeader(
                sciheader, inputFiles=sciStr[i]['fileName'], imageType='SCIENCE', masterDark=masterDarkOutFile, masterFlat=masterFlatOutFile)
            print("writing to: "+outFileName)
            fits.writeto(outFileName, correctedFrame,
                         header=sciheader, clobber=True)
            fits.writeto(outFlagFileName, mask, clobber=True)
            # fits.append(outFileName,mask)
            outFileNames.append(outFileName)
            outFlagFileNames.append(outFlagFileName)
    return correctedFrames, outFileNames, outFlagFileNames


def setup(config=None):
    '''
    This loads and parses the config file. It makes a dictionary which contains
    the keyword arguments needed to run makeCorrectedFrames

    output: a dictionary which can be used as input to
            makeCorrectedFrame(), resulting in corrected frames
            reduced in the manner described by the config file.
    '''
    # darks:
    sections = ['dark', 'flat', 'science']
    reducePars = {}
    for sec in sections:
        if config.has_section(sec):
            reducePars['do_'+sec] = (config.get(sec, 'do '+sec) == 'True')
            if sec in ['dark', 'flat']:
                if (config.get(sec, 'use master') == 'True') and config.has_option(sec, 'master file'):
                    reducePars['master_'+sec] = config.get(sec, 'master file')
            file_list = []
            filegroups = config.get(sec, 'files').split()
            for i, group in enumerate(filegroups):
                print("setting up parameters for "+sec+", file group "+group)
                fileNames = group.split(',')
                for fileName in fileNames:
                    file_list = file_list + glob.glob(fileName)
                reducePars['files_'+sec] = file_list
    reducePars['fluence'] = (config.get('science', 'fluence') == 'True')

    if config.has_option('output', 'master dark'):
        reducePars['masterDarkOutFile'] = config.get('output', 'master dark')
    if config.has_option('output', 'master flat'):
        reducePars['masterFlatOutFile'] = config.get('output', 'master flat')
    if config.has_option('output', 'output template'):
        reducePars['outFileTemplate'] = config.get('output', 'output template')
    return reducePars


def doplot(image, catalog, outfile='cat-overlay.pdf'):
    import matplotlib.pyplot as plt
    sigma_est, ncore, mode, hi, lo = pr.getBounds(data=image)
    vmin = mode - 3*sigma_est
    vmax = mode + 3*sigma_est
    fig, ax = plt.subplots(figsize=(16, 16))
    ax.imshow(image, vmin=vmin, vmax=vmax, origin='lower', cmap=plt.cm.Greys)
    keep = (catalog['IMAFLAGS_ISO'] == 0) & (catalog['FLAGS'] == 0)
    ax.plot(catalog[keep]['X_IMAGE'], catalog[keep]
            ['Y_IMAGE'], '+', ms=0.5, color='red')
    fig.savefig(outfile, format='pdf')
    pass


def main(argv):
    import argparse
    import configparser
    import socket
    parser = argparse.ArgumentParser(
        description='generate calibrated images of some data.')
    parser.add_argument('--config', default='default.config',
                        help="config file to use for processing.")
    args = parser.parse_args(argv[1:])
    print("using configuration from: "+args.config)
    config = configparser.ConfigParser()
    config.read(args.config)
    reductionArgs = setup(config=config)
    frames, frameFiles, flagFiles = makeCorrectedFrame(**reductionArgs)
    if config.has_section('sextractor'):
        if config.get('sextractor', 'do sextractor') == 'True':
            print("making catalogs.")
            for thisFrame, thisFlags in zip(frameFiles, flagFiles):
                sconfig = config.get('sextractor', 'sext config file')
                spars = config.get('sextractor', 'parameters file')
                nnw = config.get('sextractor', 'star nnw file')
                sfilter = config.get('sextractor', 'filter file')
                catalog, outfileName = makeCatalog(imageFile=thisFrame, flagFileName=thisFlags,
                                                   sconfig=sconfig, spars=spars,
                                                   nnw=nnw, sfilter=sfilter)
                image = fits.getdata(thisFrame)
                plotTemplate, ext = os.path.splitext(thisFrame)
                plotfile = plotTemplate+'-plot.pdf'
                doplot(image, catalog, outfile=plotfile)


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
