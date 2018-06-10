# bf_hxrg_ppl
Code to measure the BF effect in HXRG detectors using data from JPL's Precision Projector Laboratory at Caltech. This code was used to produce the results in Plazas et al. 2018 (arXiv:1712.06642)

## Main code:  bf_ppl.py 

### Usage of the code: 

% python bf\_ppl.py 

The code is run in `lucius`, and it reads parameters a configuration filed named `config_bf_ppl.ini`. The code analyzes PPL data (`SOLO` or `FILMTRIP` mode) and simulated data (darks, flats, and spots) created by using a combination of the codes `sim.py` and `hxrg_simulator.py` (see bellow).  

### Files that should be in the same directory where "bf\_ppl.py" is:

- `badger.py`: Slightly modified version. In particular, the lines `fileStr['FRAMTIME'][i] = thisHdr['FRAMTIME']` and `fileStr['LODFILE'][i] = thisHdr['LODFILE']` were added in the function ` organize_ramps`. 
- `sextractor_engine.py`
- `pixel_rejector.py`
- `moments.py`
- `sigma_clip.py` : From Dr. E. Sheldon
- `config_bf_ppl.ini`

### Parameters in `config_bf_ppl.ini` for `bf_ppl.py`.
The parameters for running the code should be specified in a configuration file named `config_bf_ppl.ini`. Their names and 
meanings are as follows: 

- `OutDirRoot`: Output directory path. 

- `OutDirName`: Name of output directory. Will be created if it does not exists already. All the output files (see below)
will be placed in this directory, whose path is given by `OutDirRoot`.

- `OutPDFName`: Name of the output PDF file where diagnostic and preliminary plots will be produced. 

- `SigmaCut`: Cut for sigma-clipping when averaging ramps over spots.

- `Gain`: Mean gain of the detector, in electrons per ADU.

- `YSize`: Size of the (squared) detector, in pixels (e.g., 2048)

- `NChan`: Number of channels in the detector.

- `NRef`: Number of reference pixels to use (e.g., 3).

- `CorrectNL`: Should the code correct for detector nonlinearity? Boolean (`True` or `False`).

- `PolyOrder`: Order of the polynomial to correct for nonlinearity. Should be `2` or `3`.

- `CorrectIPC`: Should the code correct for IPC? The constant kernel is denoted by `K` in the code. Boolean (`True` or `False`).

- `SubtractDark`: Should the code subtract dark images? Boolean (`True` or `False`)

- `Simulation`: Should the code use simulated ramps (`True`) or PPL data (`False`)?

- `ExamineRamps`: Should the code look for outlier ramps (`True`) or PPL data (`False`). The results will be plotted in the first three pages of `OutPDFName`.

- `DiscardRampsSpots`: List of spot ramps number, separated by a white space, that should be discarded. E.g.: `4 20 22`. 
Set to `-1` if no ramps should be discarded. 

- `DiscardRampsDarks`: List of dark ramps number, separated by a white space, that should be discarded. E.g.: `1 3 6`. 
Set to `-1` if no ramps should be discarded. 

- `DiscardRampsFlats`:  List of flats ramps number, separated by a white space, that should be discarded. E.g.: `-1`. 
Set to `-1` if no ramps should be discarded. 

- `StartFrameFlats`: Number of starting frame for ramps of flats. If you don't wish to discard any frame, set it to `0`. If you want to discard the first frame, set it to `1`. 

- `StartFrameSpots`: Number of starting frame for ramps of spots and darks. If you don't wish to discard any frame, set it to `0`. If you want to discard the first frame, set it to `1`.

- `XBorderCut`: Do not use spots whose `x` coordinate is within this number of pixels rom the detector border. E.g., (`50`)

- `YBorderCut`: Do not use spots whose `y` coordinate is within this number of pixels rom the detector border. E.g., (`50`)

- `EndSpotVector`: When doing the big loop over all spots in the stacked image, this value is the final number of spots that should be considered. E.g., `-1` represents that the last entry of the spots vector is `-1`, i.e., the whole vector. 

- `BadPixelMask`: Bad pixel mask by Dr. Eric Huff. E.g., `/projector/aplazas/master-euclid-mask.fits`. Will be use by `SExtractor` (in the parameter fie) when finding the sources.  It will also be used to discard the postage stamps of the spots that have at least 1 pixel whose value is different from `0` (good) in the mask. 

- `CentroidType`: Either `center` or `corner`.  Chose sources close to the center of the pixel or to the corner.

- `CentroidThreshold`: Distance, in pixels, from center or corner of a pixel. E.g., `0.1`.

- `RegionCorner`: If `CentroidType` is set to `corner`, this value should be set to `1`, `2`, `3`, or `4`, according to the quadrant of the `3x3` postage stamp that you wish to stack.  

- `SextractorPath`: Path to `Sextractor` executable. E.g., `/usr/local/optical/sextractor/bin/sex`.

- `ListDarksPPL`: List of PPL files with darks. E.g., `/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-000[0-9]*.fits`.

- `ListSpotsPPL`: List of PPL files with spots. E.g., `/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-02[3-9][0-9]*.fits`

- `ListFlatsPPL`: List of PPL files with flats. E.g., `/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-00[3-9][0-9]*.fits`. If you want to mitigate or get rid of the `burn-in` effect, discard the first few files at this point. 

- `ListDarksSimulation`: List of simulated files with darks. E.g., `/projector/aplazas/TESTJULY21_90_V9_OFFSET00_LOW_NOISE_NO_NL/*BACKGROUND*_00[1-2]*.fits`. 

- `ListSpotsSimulation`: List of simulated files with spots. E.g., `/projector/aplazas/TESTJULY21_90_V9_OFFSET00_LOW_NOISE_NO_NL/*OBJECT*.fits`

- `ListFlatsSimulation`: List of simulated files with flats. E.g.,`/projector/aplazas/TESTJULY21_90_V9_OFFSET00_LOW_NOISE_NO_NL/*FLAT*.fits`

### Output: 

The following ASCII files will be created by the code and placed in the directory `OutDirRoot`+`OutDirName`. Use the code `plot_fn.py` to read and plot them. The figures for the paper come mainly from the PDF file produced after running that code (`plot_fn.py`). 

- `jay_relative_size.dat`:  Data for Fig. 8 of the paper. 
   - 3 columns, one row per ramp frame. 
   - | mean size (returned by `moments.py`) |  error on mean | signal in central pixel of average between consecutive frames |

- `jay_metric.dat`: Data for Fig. 3  of the paper.
   - (1 + 3\*Nframes) columns
   - |pixel number| Nframes columns: `f_N` for each frame | Nframes columns: error on `f_N` per frame | Nframes columns: mean signal per frame |

- `jay_metric_surrounding.dat`: Data for Fig. 4 of the paper.
   -  Row 1: sum of neighbors; Row 2: central pixel. 
   - |row ID | Nframes columns: `f_N` per frame | Nframes columns: error on `f_N` per frame |

-  `jay_median_flux_flats_pixel_[1-9].dat` and `jay_residual_pixel\_[1-9]\_flat.dat`: Data for Fig. 2 of the paper.
    - One file per pixel, in each case. 
    - Each of the `jay_median_flux_flats_pixel_[1-9].dat` files has one column with the mean signal in all frames of the ramp.
    - Each of the `jay_residual_pixel\_[1-9]\_flat.dat` files has Nframes columns and Nspots rows. The green spots in Fig. 2 correspond to the data in this file, and the red curve, the median per frame (column). 

-  `jay_B.dat`: Data for Fig. 9 of the paper.
   - 8 columns:  |`fc` | `b` | `c1` | `c2` | `m` | `m_err` | `m/fc` | `new_B`; Nspots rows. 
   - The column used to create the histogram in Fig. 9 is `new_B`, which is the same as Eqn. 11 of the paper.
   - `b = 2*(c1)*(c2)/fc`, an alternative to `new_B`. The two give consisten mean values. 

- `jay_metric_no_norm.dat`: Data for Fig. 6 of the paper. Eqn. 4 of the paper with the normalization set to 1. 

- `selected_positions_centroid.dat`:  `x` and `y` positions in teh detector of the selected spots after centroid condition has been satisfied. 

- `OutPDFName`: PDF file with some preliminary and diagnostic plots. The plots for the paper are produced from the flies described above and by running `plot_fn.py`.  

- Note: To produce Fig. 7, run the code 4 times as described in part 9 below, setting the parameters `CentroidType` to `corner` and 
`CornerRegion` to `1`, `2`, `3`, or `4`, respectively. Don't forget to change the `OutDirName` every time. 


### Parts of the code (from top to bottom): 

#### 1. Plotting options 

#### 2. Function definitions 

#### 3. Parameters: 
Read in from `config_bf_ppl.py`. 

Note: For the list above, double check with the PPL log book in Dropbox. The first ten ramps of the 100 flats and spots ramps are not used, for example (to mitigate burn-in effect).  Also, I split the list of files into two lists due to the way I list them through the use of regular expressions, but this is not necessary. 

Make sure that the spots and darks have the same exposure time (same number of samples per ramp). 


#### 4. Load Data 

Here the code uses `badger.getRampsFromFiles`  so make sure you can import `badger` or that you at least have the required files in your directory. 

Time along the code is in miliseconds, because it is read from the files produced by the PPL, which record time in tehse units. Thus, sometimes we need to multiply by 1000 (c.f., `plot_fn.py`) to report time in seconds (as is done in the paper). 

#### 5. Stack data 

If the number of files is less than 40, take the median. If it is larger, split the list in 3, take the median of each part, and then the mean of the last 3 medians. This is to avoid running out of memory. 


#### 6. Switch te sign of ADU (ADU-> 2^16 -1 -ADU), subtract mean of reference pixels, convert ADU to electrons
 
#### 7. Correct for IPC with kernel K 

`GLOBAL_SPOTS`, `GLOBAL_FLATS`, and darks are vectors that contain the frames for the median ramps. 

#### 8. Run SEXtractor on last frame of median ramp GLOBAL_SPOTS if not running simulations 

Convert back to ADU the last frame, use  `daofind_sex_detection.config`  as configuration file for `SExtractor`. 

The output catalog will be placed in: 

- `out=out_dir + '/' + prefix + '_sextractor_out_last_sample_ramp_100.param'`

If you are using simulations,  the code does not run `SExtractor` and uses a catalog of positions created by the user when making the simulations. 


#### 9. Centroid calculation from last frame of spots ramp, in electrons.

Subtract the bias (`B_spots`), and then, within the loop, calculate the unweighted centroid after subtracting the local background. . If `CentroidType` is `corner`, select only those sources in a given Cartesian quadrant (Region 1 to 4). I changed this by hand and ran the code 4 times to get the data to produce Fig. 7 of the paper. 


#### 10. Big loop over sources to correct for NL (from flat fields), and calculate f_N

Loop over sources: 

- Discard if it has at least one bad pixel. 

- For each source, loop over pixels in postage stamp 

- Use function `fit_pixel_ramp` to fit a quadratic or cubic function to the ramps of the spots, flats, and darks. 

- After fitting the ramps, calculate model residuals. Fig. 2 in paper. 

- Correct stamps of darks, spots, and flats for NL by using `np.root` for the quadratic and cubic cases. Subtract darks after correcting for NL. 

- Calculate size of corrected stamp, save in a vector. 

- Loop over each pixel the corrected spot stamp: 

- Calculate signal and time difference between consecutive frames 
  
- Turn electrons into electrons per time for each difference: `rates_vec_jay=delta_sig/delta_time`.
  
- In the process, calculate  `F_i - <F_i>` to eventually produce Figure 6. 
	    
- Calculate difference in rate with respect to first frame, and then normalize to produce the `f_N` metric: `jay_metric`. The vector is  `s_vec_jay/=NORM`.
	    
- For the central pixel, calculate the coefficient B : 
    - `new_B` = (m/fc)\*(NORM/(val0\*delta\_t/1000))
		 
    - The parameters used are derived from the fit:  `m, m_err=linear_fit_m (samples, s_vec_jay, err)`
		 
    - Save that `new_B` in a vector; use those numbers to produce histogram of B in paper. 
		 
     - Note that `new_B` and `b = 2*(c1)*(c2)/fc` are consistent with each other. 


#### 11. After big loop, save files with output data 

#### 12. Calculate the mean of the size of the postage stamp in each frame; then calculate relative size to first frame

Fig. 8 of paper

#### 13. Plots: 
These won’t be the final plots in the paper. Those are produced by another code (`plot_fn.py`), using the output ASCII files listed above. 

## Code: plot_fn.py

After running `bf_ppl.py` for different configurations (e.g., simulations, PPL data center, PPL data corner in each Cartesian quadrant), a set of ASCII files is produced in the outpur directory that was especified in the configuration file. Then `plot_fn.py` reads those files to produce most of the plots that ended up in the paper. In addition, `plot_fn.py` reads files from simulations. 

### Configuration file `config_plot_fn.ini` for `plot_fn.py`

- `DoSigmaClipping`: Boolean. Should the code do sigma clipping? 

- `SigmaClippingCut`: Integer. Number of sigmas if doing sigma-clipping. 

- `FramTimePPL`: Float. `FRAMTIME` keyword in the header of the PPL files. E.g., `0.837632`.

- `NFramesPPL`: Integer. Number of final frames in the `f_n` plot. You can use `OutPDFName` to figure it out. It is not the initial number of frames in the raw data/simulated data, because the process we subtract consecutive frames (and sometimes we discard the first frame). 

- `NFramesSim`: Integer. Same as `NFramesPPL`, but for simulations. 

- `NFramesCornerPPL`: Integer. Same as `NFramesPPL`, but for the 4 cases where the centroid of the spots are close to the corner of a pixel.  

- `PPLDataDirCenter`: Location of `OutDirName` from `bf_ppl.py`, after running the code with `CentroidType` set to `center`. E.g.,  `/Users/amalagon/NL_plots/ASCII_FILES_TO_PLOT/H_FILTER_CENTER_PPL/MAR21_H_BAND_F11_CUBIC`

- `PPLDataDirCornerR1`: Location of `OutDirName` from `bf_ppl.py`, after running the code with `CentroidType` set to `corner` and `RegionCorner` to `1`. E.g., `/Users/amalagon/NL_plots/ASCII_FILES_TO_PLOT/CORNER_PPL/SECOND_RUN/REGION1_xc_lt_0_yc_lt_0`

- `PPLDataDirCornerR2`: Location of `OutDirName` from `bf_ppl.py`, after running the code with `CentroidType` set to `corner` and `RegionCorner` to `2`. E.g.,`/Users/amalagon/NL_plots/ASCII_FILES_TO_PLOT/CORNER_PPL/SECOND_RUN/REGION2_xc_gt_0_yc_lt_0`

- `PPLDataDirCornerR3`: Location of `OutDirName` from `bf_ppl.py`, after running the code with `CentroidType` set to `corner` and `RegionCorner` to `3`. E.g.,`/Users/amalagon/NL_plots/ASCII_FILES_TO_PLOT/CORNER_PPL/SECOND_RUN/REGION3_xc_lt_0_yc_gt_0`

- `PPLDataDirCornerR4`: Location of `OutDirName` from `bf_ppl.py`, after running the code with `CentroidType` set to `corner` and `RegionCorner` to `4`. E.g.,`/Users/amalagon/NL_plots/ASCII_FILES_TO_PLOT/CORNER_PPL/SECOND_RUN/REGION3_xc_lt_0_yc_gt_0`

The following parameters specify the location of the directories with files from simulations, used to produced Fig. 5 and the green histogram in Fig. 9. 

- `SimDirCenterNLCorrected`: `/Users/amalagon/NL_plots/ASCII_FILES_TO_PLOT/CENTER_SIM_NL_CORRECTED`

- `SimDirCenterNLNotCorrected`: `/Users/amalagon/NL_plots/ASCII_FILES_TO_PLOT/CENTER_SIM_NL_NOT_CORRECTED`

- `SimDirCenterNothing`: `/Users/amalagon/NL_plots/ASCII_FILES_TO_PLOT/SIMS_NOTHING`

- `SimDirCenterBHistogram`: `/Users/amalagon/NL_plots/ASCII_FILES_TO_PLOT/CENTER_SIM_BF_90RAMPS_V7`


- Note: If you use different simulations, you need to especify the fluxes per pixel for a simulated spot. This is given by the following dictionary in the code: 

`simulations_flux={'1':[402.957, 677.249, 955.99, 1240.43], '2': [1825.97 ,3078.59 ,4363.04, 5672.92], '3':[402.915, 677.282, 957.607, 1241.39], \
'4': [1827.36, 3079.58, 4361.75, 5672.8], '5':[28285.2, 46977, 65505.4, 83830] , '6': [1832.58 ,3088.91, 4374.23 ,5687.13], \
'7': [402.496, 676.444, 954.201, 1238.62], '8': [1826.36, 3079.46, 4365.73, 5677.67],  '9': [403.374 ,677.778, 956.593, 1240.58]}`


## Code: sim.py 

Uses GalSim to produce a simulated 2k by 2k scene with a grid of point sources. The number of spots depends on the size of their individual postage stamps; this can be chaged at the beginning of the code. As input, the code reads the PPL PSF model file provided by Chaz (`chazPSF_lamda1_cd3_f11_pix1_noboxcar.fits`). The FITS image will be saved in a directory called `output`. You can change this in the variable `file_name`. To change the placement of the sources, modify the variable offset as needed (e.g., offset=(ud(), ud()) for random offsets, or offset=(0.0, 0.0) for sources perfectly located at the center of the pixel). The simulated scene will be used by the code `hxrg_simulator.py` to produce simulated ramps. 

## Code: hxrg_simulator.py

Originally written by Dr. Chaz Shapiro. This version has small modifications to add BF (from the Power Law model in GalSim) and IPC. Uses as input the image created with `sim.py`.

