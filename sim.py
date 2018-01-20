import sys
import os
import math
import numpy as np
import time
import galsim
import pyfits as pf

##simulations 2, all spots at center

def main(argv):

 # NL function
 f=lambda x,beta : x - beta*x*x


 stamp_size = 64   # power of 2, please
 n_tiles=2048/stamp_size  
 pixel_scale = 1.0
 n_res= 18 # Chaz's file is 1 micron 
 random_seed = 3339201


 if not os.path.isdir('output'):
    os.mkdir('output')

 file_name = os.path.join('output','h2rg_sim_offset_norm_0_point_1.fits')

 print "hola"

 nobj=n_tiles*n_tiles
 rng = galsim.BaseDeviate(random_seed+nobj)


 # Setup the images:
 psf_image = galsim.ImageF(stamp_size * n_tiles , stamp_size * n_tiles)

 #Read PSF model
 psf_file = os.path.join('data','chazPSF_lamda1_cd3_f11_pix1_noboxcar.fits')
 psf_im  = galsim.fits.read(psf_file)
 psf = galsim.InterpolatedImage(psf_im, scale = pixel_scale, flux = 1.)


 print "hola 2"
 ix_list = []
 iy_list = []
 for ix in range(n_tiles):
    for iy in range(n_tiles):
        ix_list.append(ix)
        iy_list.append(iy)
 # Build each postage stamp:
 f=open("sim_spots_position_random_norm_0_point_1.txt", 'w')
 for k in range(nobj):
    rng = galsim.BaseDeviate(random_seed+k)
    # Determine the bounds for this stamp and its center position.
    ix = ix_list[k]
    iy = iy_list[k]
    b = galsim.BoundsI(ix*stamp_size+1 , (ix+1)*stamp_size,
                           iy*stamp_size+1 , (iy+1)*stamp_size)

    

    sub_psf_image = psf_image[b]
    print " "
    print "Center in big image: ", b.center().x, b.center().y
    #line="%g %g \n" %(b.center().x, b.center().y)
    #f.write (line)
    #ud = galsim.UniformDeviate(rng)
    #Following lines to test effect of miscentering on B histogram
    random=np.random.uniform(low=-0.1, high=0.1)
    x = random
    y = np.sqrt (0.1**2 - x**2)
    offset = (x,y)
    #offset=(ud(), ud())
    #offset=(0.0, 0.0)
    #if k%2==0: 
    #    offset=(0.0,0.0)
    #else:
        #offset=(0.5,0.5)
    #    ud = galsim.UniformDeviate(rng)
    #    offset=(ud(), ud())
    print "ix, iy, dx, dy, norm: ", ix, iy, offset, np.sqrt(x**2 + y**2)
    final_x, final_y = b.center().x + offset[0], b.center().y + offset[1]
    print "Center in image + offset: ", final_x, final_y
    final_x_int, final_y_int = np.int (np.rint(final_x)), np.int (np.rint(final_y))
    print "Center in integer form: ", final_x_int, final_y_int
    line="%g %g \n" %(final_x_int, final_y_int)
    f.write (line)

    #continue
    psf.drawImage(sub_psf_image, offset=offset, scale=pixel_scale*n_res, use_true_center=False)
    #f.close()

 f.close()


 galsim.fits.write(psf_image, file_name)


if __name__ == "__main__":
    import pdb, traceback, sys
    try:
     main(sys.argv)
    except:
     thingtype,value,tb = sys.exc_info()
     traceback.print_exc()
     pdb.post_mortem(tb)
