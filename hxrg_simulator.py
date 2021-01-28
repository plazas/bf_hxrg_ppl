#!/usr/bin/env python
'''
Simulation to generate HxRG images with which to test reduction scripts
Author: Chaz Shapiro (2017)
TODOS:  Make scriptable, add IPC?, add bad pixels?
'''

import numpy as np
import pyfits as pf
import os.path
import pdb
from scipy import ndimage
import galsim


def main(argv):

    #######################
    ### USER PARAMETERS ###

    # imgtypes = ['OBJECT','FLAT','BACKGROUND','DARK']  #Types of images to include; comment out unneeded types
    imgtypes = ['OBJECT', 'FLAT', 'BACKGROUND']

    basename = '/projector/aplazas/TESTJULY21_90_V9_OFFSET00_LOW_NOISE_NO_NL/img'
    #basename = '/projector/aplazas/PRUEBA/img'
    ramplen = 6			# number of samples up the ramp
    frametime = 3.0			# frame time in sec
    nexp = 90 			# number of ramps

    noise_seed = 123  # Changes the noise realization
    # Changes the device non-uniformity (don't change if simulating a single device)
    detector_seed = 456

    # Exposure parameters
    #exposure_file = 'test_pattern_2048x2048.fits'
    # exposure_file='./output/h2rg_sim_offset_00_random.fits'
    # exposure_file='./output/h2rg_sim_offset_00_05.fits'
    exposure_file = './output/h2rg_sim_offset_00.fits'
    # exposure_file='./output/h2rg_sim_offset_random.fits'
    # exposure_file="./output/h2rg_sim_offset_random_stamp_size_16.fits"
    exposure_file = './output/h2rg_sim_offset_norm_0_point_1.fits'
    flux = 9000 			# scale input image in e-/s/px
    flatflux = 9000 			# flat flux in e-/s/px
    # mean background in e-/s/px (affects object and flat)
    backflux = 1.0
    shotnoise_on = False		# Toggle shot noise

    # Detector parameters
    readnoise = 15.0  # e- rms per frame
    Idark = 5.0  # mean dark current in e-/s/px
    dark_scatter = 0.05  # 0.05	        #relative scatter (multiplies dark)
    n_channels = 32  # number of detector readout channels
    bias_offset = 0  # constant offset to bias frame
    # bias_file = './bias_frame_fixed.fits'  #FITS image of a bias (reset) frame
    bias_file = './bias_frame.fits'

    # Pixel response
    gain = -2.7  # e-/ADU  -- negative counts ADU down
    channel_gain_scatter = 0.05  # 0.05 	#relative scatter (multiplies gain)
    PRNU_scatter = 0.05  # pixel-wise QE variations
    # -1.e-6  		#mean 2nd order coefficient for quadratic nonlinearity model
    NL2_mean = -7.76e-7
    NL2_std = 0.34e-7  # 0.34e-7 #0.2e-6		#scatter in coefficient

    # IPC

    M = np.array([[0, -0.007, 0], [-0.009, 0.032, -0.009], [0, -0.007, 0]])
    I = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]])
    K = I-M

    # K=np.array([[1,0,0],[0,1,0],[0,0,1]])

    # BF, simple Power Law
    #galsim.cdmodel.PowerLawCD(self, n, r0, t0, rx, tx, r, t, alpha)
    """
 Initialize a power-law charge deflection model.

 The deflections from charges in the six pixels directly neighbouring a pixel border are 
 modelled independently by the parameters `r0`, `t0` (directly adjacent to borders between 
 two pixels in the same row=y / column=x) and `rx`, `tx` (pixels on the corner of pixel
 borders).

 Deflections due to charges further away are modelled as a power-law

    a = A * numpy.sin(theta) * (r_distance)**(-alpha)

 where `A` is a power-law amplitude (`r` for `a_l / a_b` and `t` `for a_b / a_t`), `theta` is
 the angle between the pixel border line and the line from border center to the other pixel
 center.

 Sign conventions are such that positive `r0`, `t0`, `rx`, `tx`, `r`, `t` correspond to
 physical deflection of equal charges (this is also how the `theta` above is defined).

 @param n      Maximum separation [pix] out to which charges contribute to deflection
 @param r0     a_l(0,-1)=a_r(0,+1) deflection coefficient along x direction
 @param t0     a_b(-1,0)=a_t(+1,0) deflection coefficient along y direction
 @param rx     a_l(-1,-1)=a_r(+1,+1) diagonal contribution to deflection along x direction
 @param tx     a_b(-1,-1)=a_t(+1,+1) diagonal contribution to deflection along y direction
 @param r      power-law amplitude for contribution to deflection along x from further away
 @param t      power-law amplitude for contribution to deflection along y from further away
 @param alpha  power-law exponent for deflection from further away
 """

    #cd = galsim.cdmodel.PowerLawCD(1, 5.e-7, 5.e-7, 1.5e-7, 1.5e-7, 2.5e-7, 2.5e-7, 1.3)
    # V1: 5, V2:4, V3:2.5, V4: 1.5, V5: 1.2, V6:1.0, V7: 1.1, V8: 1.1 + shot noise + centroid < 0.1, V9: 1.1 + no shot noise + centroid < 0.1
    cd = galsim.cdmodel.PowerLawCD(
        1, 1.1e-7, 1.1e-7, 1.0e-7, 1.0e-7, 0.0, 0.0, 0.0)

    ### END USER PARAMETERS ###
    ###########################

    ######################
    ### Meta stuff     ###

    outformat = '_'.join([basename, '%s', '%04i', '%04i']) + '.fits'
    #outformat = '_'.join(['%s','%04i','%04i']) + '.fits'

    # Store header for FITS file.  Note time units are in ms.
    header_dict = {
        'DETMODE': '32chSUTR', 'BASE': basename, 'RAMPLEN': ramplen, 'RAMP1ST': 0, 'RAMPLAST': ramplen-1, 'FRAMTIME': 1000*frametime, 'EXPTIME': (ramplen-1)*frametime*1000, 'FILMSTRP': 'False', 'CHANNELS': n_channels, 'NEXP': nexp, 'COMMENT': ''
    }
    header = pf.Header()
    for k in list(header_dict.keys()):
        header[k] = header_dict[k]

    def writeFITS(filename, data, **kwargs):
        '''Make the path if it doesn't exist.  Write the FITS file.
        '''
        path, fileonly = os.path.split(filename)

        # Using the 'make first, check later' method which avoids a race condition in parallel
        try:
            os.makedirs(path)
        except OSError:
            if not os.path.isdir(path):
                raise

        pf.writeto(filename, data, **kwargs)

    ######################
    ### Make Templates ###
    np.random.seed(seed=detector_seed)

    # Load exposure template and bias frame
    exposed_img = pf.getdata(exposure_file).astype('float64')
    bias = pf.getdata(bias_file).astype('float64') + bias_offset

    # bias=0.2*bias

    print(np.mean(exposed_img))

    print(np.mean(bias))
    # sys.exit(1)

    if exposed_img.shape != bias.shape:
        raise Exception('Bias and exposure dimensions must match')

    img_shape = bias.shape
    xsize, ysize = img_shape
    channel_width = img_shape[0]/n_channels

    x = np.arange(img_shape[0])
    y = np.arange(img_shape[1])
    [X, Y] = np.meshgrid(x, y)

    # reference pixels (do not accumulate charge)
    reference_pixels = np.zeros(img_shape)
    reference_pixels[4:-4, 4:-4] = 1.

    # Pixel response -- have the QE taper off in a corner
    # QE = np.exp(-(X/2./xsize)**2/2. -(Y/2./ysize)**2/2.)  #avoid integer division
    QE = 1

    # non-uniformity
    if PRNU_scatter > 0:
        QE *= np.random.normal(1, PRNU_scatter, size=img_shape)

    # save gains for later -- applies to noise, not just means
    channel_gains = gain*np.ones(n_channels)
    if channel_gain_scatter > 0:
        channel_gains *= np.random.normal(1.,
                                          channel_gain_scatter, size=n_channels)
    # Make copies to get a full row of gains: [g0,g0,g0...g1,g1,g1...gN,gN,gN]
    column_gains = np.tile(channel_gains, [channel_width, 1]).T.ravel()

    # Nonlinearity template
    NL_betas = NL2_mean*np.ones_like(bias)
    if NL2_std > 0:
        NL_betas += np.random.normal(0, NL2_std, size=img_shape)

    # Background pattern
    kx, ky = [2*np.pi/400, 2*np.pi/200.]
    A = .3
    background = backflux  # * (1. + A*np.sin(kx*X+ky*Y))

    # Dark current
    dark = Idark*np.ones_like(bias)
    if dark_scatter > 0:
        dark *= np.random.normal(1, dark_scatter, size=img_shape)

    ###################
    ### MAKE IMAGES ###
    # Careful when parallelizing!  Seed may get copied to each process
    np.random.seed(seed=noise_seed)

    for imgtype in imgtypes:

        header['TYPE'] = imgtype

        if imgtype == 'OBJECT':
            charge_template = (background+exposed_img*flux)*QE
        elif imgtype == 'FLAT':
            charge_template = (background+flatflux)*QE
        elif imgtype == 'BACKGROUND':
            charge_template = (background)*QE
        elif imgtype == 'DARK':
            charge_template = 0.
        else:
            raise Exception('Unknown image simulation type')
        #if imgtype=='OBJECT': pdb.set_trace()

        print("exposed_img: ", np.mean(exposed_img))
        print("flux: ", flux)

        # bias=0.2*(2**16-1-bias)  #temp
        # bias=2**16-1-bias
        # print  "bias: ", np.mean(bias)
        print("charge: ", np.mean(charge_template))
        charge_template += dark
        print(np.mean(charge_template), np.mean(dark))
        # charge_template += bias  #temp
        print(np.mean(charge_template), charge_template.shape)
        # sys.exit()
        charge_template *= frametime * reference_pixels

        ## SHOULD WE WRITE THE TEMPLATE TO A FILE TO AVOID PASSING IT TO THREADS? ##

        # variables to pass: header, charge_template, shotnoise_on, readnoise, column_gains, bias, outname

        if imgtype == 'BACKGROUND':
            nexp = 25

        for imagen in range(nexp):  # Loop over exposures

            header['IMAGEN'] = imagen

            # Create 0x, 1x, 2x... frametime worth of charge along axis 0 of array
            # We could move this step up one loop, but leave here to reduce data sent to parallel process
            cum_charge = np.array(
                [charge_template for i in range(header['RAMPLEN'])])
            cum_charge[0] = 0.

            #ytest, xtest=481-1 ,1953-1
            ytest, xtest = 161 - 1, 33 - 1
            # for a in cum_charge:
            #	print a.shape, a[ytest,xtest], np.min(a), np.max(a), np.mean(a)

            # Add shot noise to each cube slice
            if shotnoise_on:
                cum_charge = np.random.poisson(cum_charge).astype('float64')

            # Cumulatively sum charge AFTER shot noise so that later sums contain previous noise
            cum_charge = np.cumsum(cum_charge, axis=0)
            # for a in cum_charge:
            #        print a.shape, a[ytest, xtest]
            # sys.exit(1)

            #if imgtype=='OBJECT': pdb.set_trace()
            # print "Accumulated charge so far: "
            # for sample in cum_charge:
            #    print sample[ytest, xtest]
            # sys.exit()

            # TEMP: add bias here
            #cum_charge=cum_charge + bias/column_gains
            # print "Accumulated charge so far after bias addition: "
            # for sample in cum_charge:
            #    print sample[ytest, xtest]

            # Make sensed charge nonlinear -- IS THIS IN THE RIGHT PLACE??
            #cum_charge_no_nl=  cum_charge
            # print "Before applying NL: "
            # for a in cum_charge_no_nl:
            #    print "%.8f" %a[ytest, xtest]

            # Bias
            print("bias: ", (2**16-1-bias[ytest, xtest])*2.7)

            print("Before NL: ")
            for sample in cum_charge:
                print("%.7f" % sample[ytest, xtest])
            print(" ")

            cum_charge = cum_charge + NL_betas*cum_charge**2
            # print "after applying NL: "
            # for a in cum_charge:
            #    print a.shape, a[ytest, xtest]
            # sys.exit()

            # print " "

            def nl(obs, c2, c0):
                return (-1 + np.sqrt(1 - 4*c2*(c0-obs)))/(2*c2)
            c2 = NL2_mean
            # for sample in cum_charge:
            #    obs=sample[481, 1953]
            #    print obs, nl(obs, c2, 0.0)

            #if imgtype=='OBJECT': pdb.set_trace()

            ## IPC (AP)
            # print cum_charge.shape, K.shape
            temp_cum_charge = []
            for sample in cum_charge:
                print(np.mean(sample))

                #temp = sample
                # IPC
                temp = ndimage.filters.convolve(
                    sample, K, mode='constant', cval=0.0)
                # print sample[50,50], temp[50,50]

                # Apply BF
                temp = galsim.Image(temp)
                temp2 = cd.applyForward(temp)
                temp_cum_charge.append(temp2.array)

                # temp_cum_charge.append(temp)
            cum_charge = np.array(temp_cum_charge)

            # Add read noise in units of e-
            if readnoise > 0:
                cum_charge += np.random.normal(0,
                                               readnoise, size=cum_charge.shape)

            # gain (divide to convert e- to ADU)  Flips sign if gain<1
            cum_charge /= column_gains
            #cum_charge_no_nl/=  column_gains

            # Dump charge onto existing bias (reset) frame
            cum_charge += bias
            #cum_charge_no_nl = cum_charge_no_nl + bias*2.7

            # print "After bias: "
            # for a in cum_charge_no_nl:
            #    print "%.8f" %a[ytest, xtest]
            #    print "%.8f" %(bias[ytest,xtest]*2.7)
            #    print "%.8f" %(2.7*bias[ytest,xtest] + a[ytest, xtest])

            # for sample in cum_charge:
            #    obs=sample[481, 1953]
            #    c0=bias[481,1953]
            #    print obs, c0, nl(obs, c2, c0)
            # sys.exit()
            #if imgtype=='OBJECT': pdb.set_trace()

            print("Before quantization: ")
            for sample in cum_charge:
                print("%.7f" % sample[ytest, xtest])

            # Cutoff image at 0 ADU and convert to int
            # Really there is also quantization at the photon and dark current level, but ignoring that for now
            cum_charge = np.rint(cum_charge)
            cum_charge = np.clip(cum_charge, 0, 2**16-1).astype('uint16')

            print(" ")
            print("After quantization: ")
            for sample in cum_charge:
                print("%.7f, %.7f, %.7f" % (
                    sample[ytest, xtest], 2**16-1-sample[ytest, xtest], (2**16-1-sample[ytest, xtest])*2.7))

            #cum_charge_no_nl = np.clip (cum_charge_no_nl, 0, 2**16-1).astype('uint16')

            # a,b=[],[]
            # for sample, sample2 in zip(cum_charge, cum_charge_no_nl):
            #    print "sample2: %.7f " %sample2[ytest,xtest]

            #    obs=sample[ytest, xtest]
            #    c0=bias[ytest, xtest]*2.7
            #    print "c0: %.7f" %c0
            #    true=sample2[ytest, xtest] - c0
            #    print "sample2[ytest, xtest] - c0: %.7f" % (sample2[ytest, xtest] - c0)
            #    fraction=100*(nl(obs, c2, c0) - true)/true
            #    print "obs, c0, nl(obs, c2, c0), true, fraction: ", obs, c0, nl(obs, c2, c0), true, fraction
            #    a.append(true + c0)
            #    b.append(fraction)
            # print a, b
            # sys.exit()

                # Write all frames
            for i in range(header['RAMPLEN']):
                header['RAMPID'] = i
                outname = outformat % (header['TYPE'], imagen, i)
                writeFITS(outname, cum_charge[i], clobber=True, header=header)
                print('Saved ', outname)


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
