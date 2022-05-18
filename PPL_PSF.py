import galsim

pixel_scale = 1
stamp_size = 64   # power of 2, please
n_tiles = 64/stamp_size
n_res = 10  #Chaz's file is 1 micron

def write_psf_image (offset, output_file_name):
    # Setup the images:
    psf_image = galsim.ImageF(stamp_size * n_tiles, stamp_size * n_tiles)
    # Read PSF model
    psf_file = "/project/plazas/PPL/H4RG/PPL_PSF_chaz/chazPSF_lamda1_cd3_f11_pix1_noboxcar.fits" 
    psf_im = galsim.fits.read(psf_file)
    psf = galsim.InterpolatedImage(psf_im, scale=pixel_scale, flux=1.)
    #offset=(0.0, 0.0)
    psf.drawImage(psf_image, offset=offset,
                  scale=pixel_scale*n_res, use_true_center=False)
    galsim.fits.write(psf_image, output_file_name)



write_psf_image ((0.0,0.0), "PPS_PSF_no_offset.fits")
write_psf_image ((0.2,0.0), "PPS_PSF_02_offset_x.fits")
