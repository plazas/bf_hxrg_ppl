# bf_hxrg_ppl
Code to measure the BF effect in HXRG detectors using data from JPL's Precision Projector Laboratory at Caltech. This code was used to produce the results in Plazas et al. 2018 (arXiv:1712.06642)



## Code:  bf_ppl.py 

## Usage of the code: 

% python bf\_ppl.py <name\_of\_output\_directory>

The code is run in ‘lucius’. 

## Output: 

The following ASCII files will be created in the directory “out\_dir/<name\_of\_output\_directory>”. Currently, “out_dir” is set to “/projector/aplazas/” in lucius. Use the code "plot_fn.py" to read them and plot them. The figures for the paper come mainly from the PDF file produced after running that code ("plot_fn.py").  


jay_relative_size.dat:  Data for Fig. 8
jay_metric.dat: Data for Fig. 2
jay_metric_surrounding.dat: Data for Fig. 4

jay\_NORM\_spots.dat
jay\_NORM\_flats.dat
jay\_B.dat
skip.dat
jay\_c2\_flats.dat
jay\_c2\_spots.dat
jay\_residual\_center\_pixel\_flat.dat
jay\_diff\_fluxes\_center\_pixel.dat
jay\_ratio\_fluxes\_center\_pixel.dat
jay\_median\_flux\_flats\_center\_pixel.dat
jay\_median\_flux\_spots\_center\_pixel.dat
jay\_median\_deficit\_flux\_flats\_center_pixel.dat

selected\_positions\_centroid.dat:  positions of the selected spots after centroid condition has been satisfied 


for string in ["center", "n1", "n2", "n3", "n4"]:
    cmd="rm "+out_dir+"flat_calibration_"+string+".dat"; run_shell_cmd(cmd)


## Parts of the code (from top to bottom): 

#### 1. Plotting options 

### 2. Function definitions 

### 3. Parameters: 

dir = <name\_of\_output\_directory> 

pp = PdfPages(out\_dir+"bf\_ppl\_out.pdf") :this is the name of the output PDF file with plots 

sigma_cut : threshold for sigma clipping 

gain:  gain in e/ADU

ysize : y-size of the detector

nchan: number of channels 

nref=number of reference columns 

stamp_string:  Size of the postage stamp around point sources. It is either “three” or “five”. However, the “five” version of the code is not fully tested yet, so use “three”. 

correct_NL: Should we correct for NL? Boolean. 

correct_IPC: Should we correct for IPC? Boolean. 

simulation: should we use simulated data? Boolean. Set to False to use PPL data. 
root\_sim= Substring to identify the directory where the simulations are. Example: 
"/projector/aplazas/"+root\_sim+"\_OFFSET00\_LOW\_NOISE\_NO\_NL/\*BACKGROUND\*\_00[1-2]\*.fits"
See the parameters “list\_of\_???\_sims” below to modify the path as you wish. 


examine_ramps: Boolean. If true, discard individual ramps that look like outliers. 

MASK:  Bad pixel mask of the detector derived by Eric. Will be use by Sextractor (in the parameter fie) when finding the sources.  it will also be used during the Big Loop of step 10 when going to the postage stamp of each point source, discarding those stats that do not have all good pixels. 

centroid\_type = Either ‘center' or ‘corner’.  Chose sources close to the center of the pixel or to the corner. 
centroid\_threshold=distance in pixel from center or corner. Usually set to 0.1.

x\_cut, y\_cut = Exclude these pixels from the edge

SEXTRACTOR: location of extractor executable 


list\_of\_darks\_ppl: list of dark FITS files from PPL 

list\_of\_spots\_ppl\_1: list of spots FITS files from PPL, part 1 

list\_of\_spots\_ppl\_2: list of spots FITS files from PPL, part 2

list\_of\_flats\_ppl\_1: list of flats FITS files from PPL, part 1 

list\_of\_flats\_ppl\_2: list of flats FITS files from PPL, part 2 

list\_of\_darks\_sims: list of dark FITS files from simulations 

list\_of\_flats\_sims: list of flats FITS files from simulations 

list\_of\_spot\_sims: list of spots FITS files from simulations 

Note: For the list above, double check with the PPL log book in Dropbox. The first ten ramps of the 100 flats and spots ramps are not used, for example (to mitigate burn-in effect).  Also, I split the list of files into two lists due to the way I list them through the use of regular expressions, but this is not necessary. 

Make sure that the spots and darks have the same exposure time (same number of samples per ramp). 


### 4. Load Data 

Here the code uses “badger.getRampsFromFiles”  so make sure you can import “badger” or that you at least have the required files in your directory. 

### 5. Stack data 

If the number of files is less than 40, take the median. If it is larger, split the list in 3, take the median of each part, and then the mean of the last 3 medians. This is to avoid running out of memory. 


### 6. Switch te sign of ADU (ADU-> 2^16 -1 -ADU), subtract mean of reference pixels, convert ADU to electrons

Note:  A couple of parameters are defined here: 

start\_sample\_spots, end\_sample\_spots=1, shapes\_spots[0] 
start\_sample\_flats, end\_sample\_flats=1, shapes_flats[0] 

The starting frame for flats and spots is “1”, meaning that we are discarding frame “0”.  For the darks, use the same as for the spots. 


 ### 7. Discard frame 0 and correct for IPC with kernel K 

GLOBAL\_SPOTS, GLOBAL\_FLATS, and darks are vectors that contain the frames for the median ramps. 


### 8. Run SEXtractor on last frame of median ramp GLOBAL_SPOTS if not running simulations 

Convert back to ADU theist frame, use  “daofind\_sex\_detection.config”  as configuration file for SEXtractor. 

Output catalog will be placed in: 

out=out\_dir + '/' + prefix + '\_sextractor\_out\_last\_sample\_ramp\_100.param'


If you are using simulations,  the code does not run SExtractor and uses a catalog of positions created by the user when making the simulations. 


### 9. Centroid calculation from last frame of spots ramp, in electrons.

Subtract the bias (B_spots), and then, within the loop, calculate the unweighted centroid after subtracting the local background. . If centroid type is ‘corner’, select only those sources in a given Cartesian quadrant (Region 1 to 4). I changed this by hand and ran the code 4 times to get the data to produce Figure 7 of the paper. 


### 10. Big loop over sources to correct for NL (from flat fields), and calculate f_N

to_plot: ID of sources to plot in the final PDF  

Loop over sources: 

-Discard if it has at least one bad pixel. 
-For each source, loop over pixels in postage stamp 
- Use function “fit\_pixel\_ramp” to fit a quadratic function to the ramps of the spots, flats, and darks. 
- After fitting the ramps, calculate model residuals. Figure 2 in paper. 
- Correct stamps of darks, spots, and flats for NL by using the quadratic formula. Subtract darks after correcting for NL. 
- Calculate size of corrected stamp, save in a vector. 
-Loop over each pixel the corrected spot stamp: 
           -Calculate signal and time difference between consecutive frames 
	-Turn electrons into electrons per time for each difference: 
            rates\_vec\_jay=delta\_sig/delta\_time
            - In the process, calculate  F_i - <F_i> to eventually produce Figure 6. 
            - Calculate difference in rate with respect to first frame, and then normalize to produce the f_N metric: “jay metric”. The vector is  s\_vec\_jay/=NORM
            
             - For the central pixel, calculate the coefficient B : 
                           new_B = (m/fc)*(NORM/(val0*delta\_t/1000)) 
                           The parameters used are derived from the fit:  m, m_err=linear_fit_m (samples, s\_vec\_jay, err)
                            Save that “new_B” in a vector; use those numbers to produce histogram of B in paper. 
                             Note that new_B and b = 2*(c1)*(c2)/fc  are consistent with each other. 


### 11. After big loop, save files with output data 

### 12. Calculate the mean of the size of the postage stamp in each frame; then calculate relative size to first frame

Figure 8 of paper

### 13. Plots: 
These won’t be the final plots in the paper. Those are produced by another code (plot_fn.py), using the output ASCII files listed above. 






