#coding=utf-8
import numpy as np
import matplotlib
matplotlib.use('Pdf')
import subprocess as S
import sys
import os
import astropy.io.fits as pf
import collections
from scipy import optimize
from scipy import ndimage
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf  import PdfPages
import matplotlib.font_manager as fm
from scipy import stats

import matplotlib.patches as patches
import matplotlib.cm as cm


import badger
import glob

import sigma_clip
import moments 

import gc



############ PLOTTING PARAMETERS #############

plt.minorticks_on()
import matplotlib as mpl
mpl.rc('lines', linewidth=1, color='black', linestyle='-')
mpl.rc('font', family='serif',weight='normal', size=10.0 )
mpl.rc('text',  color='black', usetex=False)
mpl.rc('axes',  edgecolor='black', linewidth=1, grid=False, titlesize='x-large', labelsize='x-large', labelweight='normal',labelcolor='black')
mpl.rc('axes.formatter', limits=[-4,4])

mpl.rcParams['xtick.major.size']=7
mpl.rcParams['xtick.minor.size']=4
mpl.rcParams['xtick.major.pad']=5
mpl.rcParams['xtick.minor.pad']=5
mpl.rcParams['xtick.labelsize']= 10 #'x-large'
mpl.rcParams['xtick.minor.width']= 1.0
mpl.rcParams['xtick.major.width']= 1.0

mpl.rcParams['ytick.major.size']=7
mpl.rcParams['ytick.minor.size']=4
mpl.rcParams['ytick.major.pad']=5
mpl.rcParams['ytick.minor.pad']=5
mpl.rcParams['ytick.labelsize']= 10 #'x-large'
mpl.rcParams['ytick.minor.width']= 1.0
mpl.rcParams['ytick.major.width']= 1.0
mpl.rc ('legend', numpoints=1, fontsize='x-large', shadow=False, frameon=False)


prop = fm.FontProperties(size=5)
loc_label='upper right'


################################## FUNCTION DEFINITIONS ################################


def fitfunc_m (x, m):
    x=np.array(x)
    return m*x
#pinit=[0.1]

def linear_fit_m (x,y,y_err):
    pinit=[-0.001]
    pfinal, covar=optimize.curve_fit(fitfunc_m,x, y, p0=pinit, sigma=y_err,  maxfev=100000)
    return pfinal[0], np.sqrt(covar[0])



### QUADRATIC POLYNOMIAL
fit_order=2
def fitfunc (x, p0, p1, p2):
    x=np.array(x)
    return p0 + p1*x + p2*(p1*x)**2
#pinit=np.repeat([1.0], fit_order + 1)
pinit=[10000, 1000, -7.76e-7]
pinit2=[10000,1000,-7e-10]


#0 0 1 1
#0 1 2 2
#0 2 3 3
#1 0 4 4
#1 1 5 5
#1 2 6 6
#2 0 7 7
#2 1 8 8
#2 2 9 9

dict_3x3={(0,0):(-1,1), (0,1):(0,1), (0,2): (1,1), \
          (1,0):(-1,0), (1,1):(0,0), (1,2): (1,0), \
          (2,0):(-1,-1),(2,1):(0,-1), (2,2): (1,-1) }

def get_centroid_3x3 (stamp):
    s1,s2=stamp.shape
    if not s1 == 3 and s2 == 3: 
	print "Error: stamp must be 3X3"
        sys.exit()
    xcent,ycent=0,0
    sum_stamp=np.sum(stamp)
    for (j,i), val in np.ndenumerate(stamp):
        xcent+=dict_3x3[(j,i)][0]*val
        ycent+=dict_3x3[(j,i)][1]*val
    return xcent/sum_stamp, ycent/sum_stamp

#(0, 0) 1
#(0, 1) 2
#(0, 2) 3
#(0, 3) 4
#(0, 4) 5
#(1, 0) 6
#(1, 1) 7
#(1, 2) 8
#(1, 3) 9
#(1, 4) 10
#(2, 0) 11
#(2, 1) 12
#(2, 2) 13
#(2, 3) 14
#(2, 4) 15
#(3, 0) 16
#(3, 1) 17
#(3, 2) 18
#(3, 3) 19
#(3, 4) 20
#(4, 0) 21
#(4, 1) 22
#(4, 2) 23
#(4, 3) 24
#(4, 4) 25

dict_5x5={(0,0):(-2,2), (0,1):(-1,2), (0,2):(0,2), (0,3):(1,2), (0,4):(2,2), \
      (1,0):(-2,1), (1,1):(-1,1), (1,2):(0,1), (1,3):(1,1), (1,4):(2,1), \
      (2,0):(-2,0), (2,1):(-1,0), (2,2):(0,0), (2,3):(1,0), (2,4):(2,0), \
      (3,0):(-2,-1), (3,1):(-1,-1), (3,2):(0,-1), (3,3):(1,-1), (3,4):(2,-1), \
      (4,0):(-2,-2), (4,1):(-1,-2), (4,2):(0,-2), (4,3):(1,-2), (4,4):(2,-2) }

def get_centroid_5x5 (stamp):
    s1,s2=stamp.shape
    if not s1 == 5 and s2 == 5:
        print "Error: stamp must be 5X5"
        sys.exit()
    xcent,ycent=0,0
    sum_stamp=np.sum(stamp)
    for (j,i), val in np.ndenumerate(stamp):
        xcent+=dict_5x5[(j,i)][0]*val
        ycent+=dict_5x5[(j,i)][1]*val
    return xcent/sum_stamp, ycent/sum_stamp


a=np.array([[3,10,3],[10,50,10],[3,10,3]])
print a
s,t=get_centroid_3x3(a)
print s, t

a=np.array([[7.5,10,3],[10,50,10],[3,10,3]])
print a
s,t=get_centroid_3x3(a)
print s, t, np.sqrt (s**2 + t**2)

a=np.array([[1,10,4.3],[6.8,50,50],[7.6,50,50]])
print a
s,t=get_centroid_3x3(a)
print s, t, np.sqrt (s**2 + t**2)


#sys.exit()




def fit_pixel_ramp (ramp='', time='', i=0, j=0, pinit=pinit):
     if not len (ramp) == len (time):
         print "inside function ;fit_pixel_ramp': len(ramp) not equal to len (time)."
         print "len (ramp) == len (time): ", len (ramp), len (time)
         print "ramp: ", ramp, ramp.shape
         sys.exit()
         print "time: ", time
         sys.exit()   

 
     flag=False
    

     time_vec, signal_vec, signal_vec_err=[],[],[]
     #a=[]
     #b=[]
     #counter=0
     for t, sample_array in zip(time,ramp):

        index_x, index_y = j, i  #np.unravel_index (counter, sample_array.shape)
        s=sample_array[index_x, index_y]

        time_vec.append(t)
        #signal_vec.append( 2**16-1-s )
        signal_vec.append(s) # s=ADU_dark - ADU_data
        #print "np.sqrt(s): ", np.sqrt(s)
        signal_vec_err.append(np.sqrt(1))    ##### TOMATO
        #if not counter == 0:
        #    a.append(time_vec[counter]-time_vec[counter-1])
        #    b.append(signal_vec[counter]-signal_vec[counter-1])
        #counter+=1
     time_vec=np.array(time_vec)
     signal_vec=np.array(signal_vec)
     signal_vec_err=np.array(signal_vec_err)

     #a=np.array(a)
     #b=np.array(b)
     print "time_vec, signal_vec: ", time_vec, signal_vec
     #print "diff time, diff signal, ratio: ", a, b, b/a
     if signal_vec[0] < 0: 
	print "Inside fit_pixel_ramp, negative signal: ", signal_vec
        signal_vec=np.fabs(signal_vec)

     if np.isnan(signal_vec).any() or np.isinf(signal_vec_err).any():
        flag=True
        print "FLAG TRUE 1"
        return [0.,0.,0.], np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]), 0., flag, signal_vec
     else:
        pfinal, covar, chi2_red = get_final_parameters_first_fit (x=time_vec, y=signal_vec, y_err=signal_vec_err, pinit=pinit)
        if np.isinf(pfinal).any(): #or np.isinf(covar).any():
            flag=True
            print "FLAG TRUE 2"
            print "pfinal, covar: ", pfinal, covar
            print "time_vec, y=signal_vec, y_err=signal_vec_err", time_vec, signal_vec, signal_vec_err
            sys.exit()
            return [0.,0.,0.],  np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]), 0., flag, signal_vec
        else:
           
            p0, p1, p2 = pfinal[0], pfinal[1], pfinal[2]
            p0_err=np.sqrt(np.diag(covar))[0]
            p1_err=np.sqrt(np.diag(covar))[1]
            p2_err=np.sqrt(np.diag(covar))[2]
            print "Todo bien in fit pixel ramp: pfinal, covar, chi2_red, flag: ", pfinal, covar, chi2_red, flag
            return pfinal, covar, chi2_red, flag, signal_vec


def correct_and_fit_pixel_ramp (ramp='', dark_ramp='', time='', i=0, j=0, pinit=pinit, c0_spot=0, c0_dark=0, c1_spot=0, c1_dark=0, c2=1,):
    #c2=-7.76e-7 
    if not len (ramp) == len (time):
         print "inside function 'correct_and_fit_pixel_ramp': len(ramp) not equal to len (time)."
         print "len (ramp), len (time): ", len (ramp), len (time)
         sys.exit()
    
    flag=False
    time_vec, signal_vec, signal_vec_err, sum_stamp=[],[],[],[]
    if not len(time)==len(ramp):
        print "len(time)==len(ramp) not satisfied"
        sys.exit(1)
    a,b, delta_sum_stamp, counter=[],[],[],0
    for t, sample_array, sample_dark in zip(time,ramp, dark_ramp):
        index_x, index_y =  j, i #np.unravel_index (counter, sample_array.shape)
        ## SPOT
  
        s_temp=sample_array[index_x, index_y]
        #print "s before NL correction: ", s_temp
        #print "c0_spot, c0_dark, c2: ", c0_spot, c2_dark, c2
        
        #s=c0_spot + c1_spot*t
        s=(1/(2*c2))*(-1+np.sqrt(1-4*c2*(c0_spot-s_temp)))         ### Apply the correction here!
        #print "signal: %g, corrected signal: %g" %(s_temp, s)
        #print "s/s2: ", s/s2, s, s2
        #sys.exit(1)
        ## Same pixels, in dark fields
        d_temp=sample_dark[index_x,index_y]
        d=(1/(2*c2))*(-1+np.sqrt(1-4*c2*(c0_dark-d_temp))) 
        #d=c0_dark + c1_dark*t
         
        ##### TEMP, turn off NL corection 
        if correct_NL == False:
            s=s_temp
            d=d_temp
            
        print " "
        print "c0_spot, c0_dark, c2 (from flat) ", c0_spot, c0_dark, c2
        print "signal: %g, corrected_signal: %g " %(s_temp, s)
        print "dark: %g, corrected_dark: %g " %(d_temp, d) 

        ## Subtract dark from data
        s-=d
        
        print "dark subtracted signal (corrected and not corrected): %g, %g" %(s, s_temp-d_temp)
        #print bias_frame[j,i]
        #sys.exit(1)
        time_vec.append(t)
        #signal_vec.append( 2**16-1-s )
        signal_vec.append(s) # s=ADU_dark - ADU_data
        signal_vec_err.append(np.sqrt(1))   ## TOMATO2
        sum_stamp.append(np.sum(sample_array))

        if not counter == 0:
            a.append(time_vec[counter]-time_vec[counter-1])
            b.append(signal_vec[counter]-signal_vec[counter-1])
            delta_sum_stamp.append(sum_stamp[counter]-sum_stamp[counter-1])
        delta_sum_stamp.append(sum_stamp)
        counter+=1
    print " "
    print "j, i: ", j, i
    print " "
    #if j == 2 and i == 2: 
    #    sys.exit()

    time_vec=np.array(time_vec)
    signal_vec=np.array(signal_vec)
    signal_vec_err=np.array(signal_vec_err)

    delta_sum_stamp=np.array(delta_sum_stamp)
    a=np.array(a)
    b=np.array(b)
    r=b/a
   
    #r=signal_vec
 

    #print "time_vec, corrected signal_vec, signal/time: ", time_vec, signal_vec, signal_vec/time_vec
    print "diff time, diff signal, ratio of diffs ", a, b, b/a


    if np.isnan(signal_vec).any() or np.isinf(signal_vec_err).any():
        flag=True
        print "inside corract_and fit function, NAN ", signal_vec, signal_vec_err
        #sys.exit()
        return [0.,0.,0.], np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]) , 0., flag, r, b, delta_sum_stamp, signal_vec
    else:
        pfinal, covar, chi2_red = get_final_parameters_first_fit (x=time_vec, y=signal_vec, y_err=signal_vec_err, pinit=pinit)
        if np.isinf(pfinal).any() or np.isinf(covar).any():
            c = np.polyfit (time_vec, signal_vec, 2)
            pfinal=c
            pfinal[0]=c[2]
            pfinal[1]=c[1]
            pfinal[2]=c[0]/c[1]**2
            return pfinal, np.array([[1.,1.,1.],[1.,1.,1.],[1.,1.,1.]]), 1, flag, r, b, delta_sum_stamp, signal_vec  
            #flag=True
            #print "inside corract_and fit function, INF ", pfinal, covar
            #sys.exit() 
            #return [0.,0.,0.], np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]), 0., flag, r, b, delta_sum_stamp, signal_vec
        else:
            p0, p1, p2 = pfinal[0], pfinal[1], pfinal[2]
            p0_err=np.sqrt(np.diag(covar))[0]
            p1_err=np.sqrt(np.diag(covar))[1]
            p2_err=np.sqrt(np.diag(covar))[2]
            return pfinal, covar, chi2_red, flag, r , b, delta_sum_stamp, signal_vec


def get_final_parameters_first_fit (x=[], y=[], y_err=[], pinit=pinit):
    ###  this function gets all the parameters according to the polynomial order: alpha, beta, gamma, delta, etc
    x,y,y_err = np.array(x), np.array(y), np.array(y_err)
    
    print "Inside get_final_parameters_first_fit (x, y, y_err): ", x,y,y_err
    for i in range(len(y)):
	if y[i] < 0: 
            y[i] = 0
            y_err[i] = 1

    pfinal, covar=optimize.curve_fit(fitfunc,x, y, p0=[1,1,1], sigma=y_err,  maxfev=1000000000)
    chi2= np.power( ( fitfunc ( x, *pfinal) - y)/y_err,  2).sum()
    n,d=len(x), len(pinit)
    nu=n-d
    chi2_red=chi2*1./nu
    #c=np.polyfit (x, y, 2)
    #pfinal=c
    #pfinal[0]=c[2]
    #pfinal[1]=c[1]
    #pfinal[2]=c[0]/c[1]**2 
    #chi2_red=1
    #covar=np.array([[1,1,1],[1,1,1],[1,1,1]])
    return pfinal, covar, chi2_red


def aux_quadratic_fit (t_vec, s_vec, s_e_vec,pinit=pinit, label="_"):
    p, c, chi2 = get_final_parameters_first_fit (x=t_vec, y=s_vec, y_err=s_e_vec, pinit=pinit)
    pe=np.sqrt(np.diag(c))
    print "p, pe, chi2: ", p, pe, chi2
    p_all=np.array([p,pe])
    #np.savetxt (outdir+"parameters_quadratic_%s.txt"%label, p_all)
    return p, pe, chi2, c

def get_residual_error (l='', s='', varc0='', varc1='', covc0c1='', t=''):
    #t=t/1000
    #return 100*np.sqrt((s/l**2)+ s**2*(varc0+t**2*varc1 + 2*t*covc0c1)/l**2)
    return 100*np.sqrt(s)/l


def run_shell_cmd (cmd):
    print >>sys.stderr, cmd
    S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()



def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()


################################## PARAMETERS ################################


dir=sys.argv[1]
out_dir='/projector/aplazas/'+dir+"/"
cmd="mkdir -v %s"%(out_dir)
S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
print "OUTPUT DIRECTORY: ", out_dir
pp=PdfPages(out_dir+"bf_ppl_out.pdf")

sigma_cut=3.0
#prop = fm.FontProperties(size=5)
#loc_label='upper right'
gain=2.7 #e/ADU


ysize=2048
nchan=32
colsperchan=ysize/nchan
nref=3
gain=2.7

#ytest, xtest=481-1 ,1953-1
ytest, xtest = 10, 10

stamp_string='three' 

correct_NL=True
correct_IPC=True
diff_spots=True
simulation=True
examine_ramps=True

SUM_STAMP_2D=[]
SUM_STAMP_2D_FLATS=[]
SUM_STAMP_2D_DARKS=[]

if simulation == True:
    root="TESTJULY21_90_V9"

    ##### Read in bias frame used in the simulation, to compare with the C0 coefficient
    #bias_file='/home/aplazas/simulations/bias_frame_fixed.fits'
    bias_file='/home/aplazas/simulations/bias_frame.fits'
    bias_frame_temp = (2**16 - 1 - pf.getdata(bias_file).astype('float64'))*gain
   

    #t_bias=bias_frame_temp
    #diffref= t_bias[-nref:].reshape(nref,nchan,colsperchan)
    #diffref = diffref.mean(2).mean(0)  #average last 4 rows in each channel
    #diffref = np.tile(diffref,(colsperchan,1)).T.flatten() #[c0,c0,...,c1,c1,...]        
    #t_bias-=diffref
    #bias_frame = t_bias*gain
    bias_frame=bias_frame_temp    
    
    ytest, xtest=481-1 ,1953-1
    ytest, xtest= 32-1-2,32-1 -2
    print "bias shape: "
    print "bias at test pixel: ", bias_frame[ytest, xtest]

### Andres's data, 03-02-17

### DATA
if simulation == False:
    files_darks=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-000[0-9]*.fits")  # 10 ramps of darks 
    ### PPL data for SURF project, see log
    #files_darks=glob.glob("/projector/cshapiro/data-raw/WFIRST/2017-06-08/raw/andres/ROIREAD_400/img-000[1-5]*.fits")


else: ## SIMULATIONS:
    #files_darks=glob.glob("/projector/aplazas/"+root+"_OFFSET00_LOW_NOISE_NO_NL/*BACKGROUND*.fits")
    files_darks=glob.glob("/projector/aplazas/"+root+"_OFFSET00_LOW_NOISE_NO_NL/*BACKGROUND*_00[1-2]*.fits")  ## Only 10 ramps/out of 90
    
    #files_darks=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-000[0-9]*.fits") 
    #files_darks=glob.glob("/projector/cshapiro/data-raw/WFIRST/2017-06-08/raw/andres/ROIREAD_400/img-000[1-5]*.fits")

if simulation == False: 
    ## PPL DATA
    #Spots:
    files1=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-02[1-9][0-9]*.fits")  # andres-02[1-9][0-9]*.fits, 80 ramps of spots, wlamp=180 SPOTS
    files2=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-030[0-9]*.fits")  # remaining 10 (initially 100 ramps, discarded first 10).   


    ## TEMP. Only a few ramps
    #files1=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-02[8-9][0-9]*.fits")




    #files2=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-030[0-9]*.fits")    ### TOTAL: ramps 210 to 309


    #### 
    #files1=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-03[6-9][0-9]*.fits")  # andres-02[1-9][0-9]*.fits, 100 ramps of spots, wlamp=180
    #files2=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-040[0-9]*.fits")
    
    #### 6-10, spots, ROI 400, SURF data, see log
    #files1 = glob.glob("/projector/cshapiro/data-raw/WFIRST/2017-06-08/raw/andres/ROIREAD_400/img-00[01][67890]*.fits")

    #files1=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-01[3-9][0-9]_000*.fits") 

    ## Flats at 100W, taking the place of SPOTS. Temp: 
    #files1=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-01[2-9][0-9]_000*.fits") 


    files_spots=sorted(files1+files2)


    #files1=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-01[12][0-9]*.fits")  # andres-01[1-9][0-9]*.fits100 ramps of flats, wlamp=100
    #files2=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-020[0-9]*.fits")
    #files_flats=sorted(files1+files2)


    files1=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-00[1-9][0-9]_000*.fits")   ## Pick flats away from first ramps to avoid "burn-in", flats at 120W. Discard first 10.
    files2=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-010[0-9]_000*.fits")  # remaining 10, 120W  

    ### TEMP
    #files1=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-00[8-9][0-9]_000*.fits")

    files_flats=sorted(files1+files2)

    #files_spots=files_flats
    ## Use Stacked (100 files) file
    #files_spots=glob.glob ("/projector/aplazas/stacked/stacked100Spots-*.fits")
    #files_flats=glob.glob ("/projector/aplazas/stacked/stacked100FlatsWLAMP120W-*.fits")

else:
    ## Use simulations 
    files_spots=glob.glob("/projector/aplazas/"+root+"_OFFSET00_LOW_NOISE_NO_NL/*OBJECT*.fits")
    files_flats=glob.glob("/projector/aplazas/"+root+"_OFFSET00_LOW_NOISE_NO_NL/*FLAT*.fits")

    
    #files_flats=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-00[4-9][0-9]_000*.fits")  #10-109: 120W 
    #files_spots=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-01[4-9][0-9]_000*.fits")  #110-209: 100W
    #files_spots=files_flats
   

allDarks, infoDarks = badger.getRampsFromFiles(sorted(files_darks))
#temp_darks=np.median(allDarks, axis=0)
files_darks=0


### Initial reductions of spots, flats, and darks: BIAS and IPC
#temp_spots=np.array(data_spots)
#temp_flats=np.array(data_flats)

allSpots, infoSpots = badger.getRampsFromFiles(sorted(files_spots))
files_spots=0

allFlats, infoFlats = badger.getRampsFromFiles(sorted(files_flats))
files_flats=0


#allSpots=allDarks


def plot_ratio_ramps (all, title='', discard=[]):
  
  num_plots= len(all)
  colormap = plt.cm.gist_ncar
  plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])


  fig=plt.figure()
  ax=fig.add_subplot(111)
  first_vec=[]
  labels=[]
  new_ramp_list=[]
  for i, ramp in enumerate((all)): 
    if not len(discard) == 0:
        #discard=len(all)-np.array(discard)
        if i in discard: continue

    temp_vec=[]
    for sample in ramp: 
        sample= 2**16 - 1 - sample
        mean_region= np.mean (sample[200:1900, 200:1900])
        if i == 0:
            first_vec.append(mean_region)    
        else:
            temp_vec.append(mean_region)
    if i == 0: 
        first_vec=np.array(first_vec)
        continue
    else:
        temp_vec=np.array(temp_vec)
    #ratio = 1- temp_vec/first_vec 
    ratio = temp_vec


    new_ramp_list.append(ramp)

    #if ratio[-2] < ratio[1]:
    #    print "ESTE ES!!!!!! : "
    #    print i, ratio, title

    plt.plot(ratio, '.-' )
    ax.annotate('%g' %i, xy=(5, ratio[5]), xytext=(5, ratio[5]), size=4)
    labels.append("%g" %(i))
    ax.set_title(title)
  #plt.legend(labels, ncol=7, loc='upper center', bbox_to_anchor=[0.5, 1.1], columnspacing=1.0, labelspacing=0.0, handletextpad=0.0, handlelength=1.5, fancybox=True, shadow=True)
  pp.savefig()
  

  return new_ramp_list


### Discard first ramps for burn-in, 
if examine_ramps == True: 
  allDarks=plot_ratio_ramps (allDarks, title='Darks', discard=[6])
  allFlats=plot_ratio_ramps (allFlats, title='Flats', discard=range(1,15)) #[1,2,3,4,5,6,7,8,9,10])
  allSpots=plot_ratio_ramps (allSpots, title='Spots', discard=range(1,15)+[72]) #1,2,3,4,5,6,7,8,9,10,21,24,25]+[26,27,28,30,32,33,34,36,37,38,42,46,54])

  #pp.close()
  #sys.exit()



#files_spots.close()
#files_flats.close()
#files_darks.close()

print "len all Spots, allFlats, allDarks:", len(allSpots), len(allFlats), len(allDarks)
print "infoSpots: ", infoSpots





"""
midpoint=50

a=np.sum(allSpots[:midpoint], axis=0)/len(allSpots)
b=np.sum(allSpots[midpoint:], axis=0)/len(allSpots)
temp_spots=a+b


midpoint=len(allFlats)/2
a=np.sum(allFlats[:midpoint], axis=0)/len(allFlats)
b=np.sum(allFlats[midpoint:], axis=0)/len(allFlats)
temp_flats=a+b
"""

### Try mean of medians


if len(allSpots) < 40:
    temp_spots=np.median(allSpots, axis=0)
    temp_flats=np.median(allFlats, axis=0)
else:
    midpoint=len(allSpots)/3
    a= np.median(allSpots[:midpoint], axis=0)
    b= np.median(allSpots[midpoint:2*midpoint], axis=0)
    c= np.median(allSpots[2*midpoint:], axis=0)
    temp_spots= (a+b+c)*1./3

    a= np.median(allFlats[:midpoint], axis=0)
    b= np.median(allFlats[midpoint:2*midpoint], axis=0)
    c= np.median(allFlats[2*midpoint:], axis=0)
    temp_flats=(a+b+c)*1./3



if len(allDarks) < 40: 
    temp_darks=np.median(allDarks, axis=0)
else:
    a= np.median(allDarks[:midpoint], axis=0)
    b= np.median(allDarks[midpoint:2*midpoint], axis=0)
    c= np.median(allDarks[2*midpoint:], axis=0)
    temp_darks=(a+b+c)*1./3    



#print temp_spots[-1][100,100]

#for i in range(len(allSpots[:midpoint])):
#    print temp_spots[-1][100,100], allSpots[i][-1][100,100], temp_flats[-1][100,100], allFlats[i][-1][100,100] 
    

shapes_spots = temp_spots.shape   # (6, 2048, 2048) --> (sample # in ramp, y, x)
shapes_flats = temp_flats.shape  
shapes_darks = temp_darks.shape


start_sample_spots, end_sample_spots=1, shapes_spots[0] 
start_sample_flats, end_sample_flats=1, shapes_flats[0]  


data_spots=np.zeros(shapes_spots)
data_darks=np.zeros(shapes_darks)
data_flats=np.zeros(shapes_flats)

print "shapes_spots, shapes_darks, shapes_flats: ", shapes_spots, shapes_darks, shapes_flats


if not shapes_spots [0] == shapes_darks[0]:
    print "Mean spot ramp and mean dark ramps don't have same number of samples/frames. "
    sys.exit()


### CHAZ" 
#diffref = diff[-args.nref:].reshape(args.nref,nchan,colsperchan)
#diffref = diffref.mean(2).mean(0)  #average last 4 rows in each channel
#diffref = np.tile(diffref,(colsperchan,1)).T.flatten() #[c0,c0,...,c1,c1,...]
#for i in range(ysize): diff[i]-=diffref

ysize=2048

if not shapes_spots[2] == ysize:
    print "ysize in spots ramp not 2048"
    sys.exit(1)

if not shapes_darks[2] == ysize:
    print "ysize in darks ramp not 2048"
    sys.exit(1)

if not shapes_flats[2] == ysize:
    print "ysize in flats ramp not 2048"
    sys.exit(1)



nchan=32
colsperchan=ysize/nchan
nref=3

### Swicth sign and correct for reference region 
for sample in range(shapes_spots[0]):
        t_spots = (2**16-1-temp_spots[sample,:,])
        #Spots
        diffref= t_spots[-nref:].reshape(nref,nchan,colsperchan)
        diffref = diffref.mean(2).mean(0)  #average last 4 rows in each channel
        diffref = np.tile(diffref,(colsperchan,1)).T.flatten() #[c0,c0,...,c1,c1,...]        
        for i in range(ysize): 
            t_spots[i]-=diffref
        data_spots[sample,:,:] = t_spots
        
for sample in range(shapes_darks[0]):
        t_darks = (2**16-1-temp_darks[sample,:,])
        diffref= t_darks[-nref:].reshape(nref,nchan,colsperchan)
        diffref = diffref.mean(2).mean(0)  #average last 4 rows in each channel
        diffref = np.tile(diffref,(colsperchan,1)).T.flatten() #[c0,c0,...,c1,c1,...]        
        for i in range(ysize): 
            t_darks[i]-=diffref        
        data_darks[sample,:,:] = t_darks

for sample in range(shapes_flats[0]):
        t_flats = (2**16-1-temp_flats[sample,:,]) 
        #Flats
        diffref= t_flats[-nref:].reshape(nref,nchan,colsperchan)
        diffref = diffref.mean(2).mean(0)  #average last 4 rows in each channel
        diffref = np.tile(diffref,(colsperchan,1)).T.flatten() #[c0,c0,...,c1,c1,...]        
        for i in range(ysize): 
            t_flats[i]-=diffref
        data_flats[sample,:,:] = t_flats

#print "data_darks after switching sign and correcting for ref: ",
#for sample in data_darks:
#        print sample[ytest, xtest]


print "data_spots after switching sign: ",
for sample in data_spots:
        print sample[ytest, xtest]


print "data_flats after switching sign: ",
for sample in data_flats:
        print sample[ytest, xtest]


#Convert to electrons, and cut the domain in samples  (e.g., to avoid transisnts at the begining and saturation at the end) 
data_spots= (data_spots)*gain
data_flats= (data_flats)*gain
data_darks= (data_darks)*gain



print "data_darks after converting to electrons: ",
for sample in data_darks:
        print sample[ytest, xtest]


print "data_spots after converting to electrons: ",
for sample in data_spots:
        print sample[ytest, xtest]


print "data_flats after converting to electrons: ",
for sample in data_flats:
        print sample[ytest, xtest]



## Reference frame. Discard frame 0. 
## For NL, don't subtract: the C0 will take care of it. 
## For detecting spots through unweighted centroid: you need to subtract it 
B_spots=data_spots[start_sample_spots]
B_flats=data_flats[start_sample_flats]
B_darks=data_darks[start_sample_spots]

print np.mean(B_spots), np.mean(B_flats), np.mean(B_darks)
print B_spots.shape
print B_spots



#Uncomment bias if not IPC

data_spots= data_spots[start_sample_spots:end_sample_spots]  #- B_spots  # subtract bias here only to select objects. When adding IPC, see below
data_flats= data_flats[start_sample_flats:end_sample_flats]  #- B_flats
data_darks= data_darks[start_sample_spots:end_sample_spots]  #- B_darks


print "data_darks: ",
for sample in data_darks:
        print sample[ytest, xtest]


print "data_spots"
for sample in data_spots:
        print sample[ytest, xtest]


print "data_flats"
for sample in data_flats:
        print sample[ytest, xtest]




M=np.array( [[0,0.007,0],[0.009,-0.032,0.009],[0,0.007,0]] )
I=np.array([[0,0,0],[0,1,0],[0,0,0]])
K=I-M
print "K: "
print K


## Loop over each sample, correct for IPC, bias, and Reference pixels
GLOBAL_SPOTS, GLOBAL_FLATS= [], []
darks=[]

if not len(data_spots) == len(data_darks):
	print "len of data_spots and data_darks is different"
        sys.exit()

for sample_s, sample_d in zip(data_spots, data_darks):
    print " "
    con_s = ndimage.filters.convolve (sample_s, K, mode='constant', cval=0.0)
    con_d = ndimage.filters.convolve (sample_d, K, mode='constant', cval=0.0)
    if correct_IPC==False: 
        con_s, con_d = sample_s, sample_d
    GLOBAL_SPOTS.append(con_s) 
    darks.append(con_d) 

for sample_f in data_flats:
    print " "
    con_f = ndimage.filters.convolve (sample_f, K, mode='constant', cval=0.0)
    if correct_IPC==False:
        con_f= sample_f
    GLOBAL_FLATS.append(con_f) 


GLOBAL_SPOTS=np.array(GLOBAL_SPOTS) #discard last
GLOBAL_FLATS=np.array(GLOBAL_FLATS)
darks=np.array(darks)
    

### Should be one ramps
#allFlats, infoFlats = badger.getRampsFromFiles(files_flats)
#allSpots, infoSpots = badger.getRampsFromFiles(files_spots)
#allDarks, infoDarks = badger.getRampsFromFiles(files_darks)


allSpots=0
allDarks=0
allFlats=0

#print files_flats
#print files_spots
#print "Number of flat ramps: ", len(allFlats)
#print "Number of spots ramps: ", len(allSpots)
#print "Number of dark ramps: ", len(allDarks)

if not infoDarks['EXPTIME'][0] == infoSpots['EXPTIME'][0]: 
    print "Exposure time for spots and darks is not the same."
    sys.exit(1)


print "infoDarks['sample'][0]: ", infoDarks['sample'][0] 
print "infoFlats['sample'][0]: ", infoFlats['sample'][0]

print "shapes_darks[0]: ", shapes_darks[0]
print "shapes_flats[0]: ", shapes_spots[0]

#sys.exit()

time_darks=np.linspace(0.0, infoDarks['EXPTIME'][0], shapes_darks[0])
time_darks=time_darks[start_sample_spots:end_sample_spots]

time_flats=np.linspace(0.0, infoFlats['EXPTIME'][0], shapes_flats[0])
time_flats=time_flats[start_sample_flats:end_sample_flats]


print "len(time_darks), len(time_flats): ", len(time_darks), len(time_flats)
print "len(GLOBAL_SPOTS): ", len(GLOBAL_SPOTS)


if not len(time_flats) == len(GLOBAL_FLATS):
    print "len(time_flats) not the same as len(GLOBAL_FLATS)"
    print len(time_flats), len(GLOBAL_FLATS)
    sys.exit()


if not len(time_darks) == len(GLOBAL_SPOTS):
    print "len(time_darks) not the same as len(GLOBAL_SPOTS)"
    print len(time_darks), len(GLOBAL_SPOTS)
    sys.exit()

if not len(time_darks) == len(darks):
    print "len(time_darks) not the same as len(GLOBAL_darks)"
    print len(time_darks), len(darks)
    sys.exit()





print "len(time_darks), len(time_flats), len(GLOBAL_SPOTS), len(GLOBAL_FLATS), len(darks)", len(time_darks), len(time_flats), len(GLOBAL_SPOTS), len(GLOBAL_FLATS), len(darks)


#print allSpots

#sys.exit(1)

#title_string=infoRamps['fileName'][0]
#title_string2=infoRamps2['fileName'][0]

#dir=sys.argv[1]
#out_dir='/projector/aplazas/'+dir+"/"
#cmd="mkdir -v %s"%(out_dir)
#S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()

#print "OUTPUT DIRECTORY: ", out_dir

cmd="rm "+out_dir+"jay_NORM_spots.dat"; run_shell_cmd(cmd)
cmd="rm "+out_dir+"jay_NORM_flats.dat"; run_shell_cmd(cmd)
cmd="rm "+out_dir+"jay_B.dat"; run_shell_cmd(cmd)
cmd="rm "+out_dir+"skip.dat"; run_shell_cmd(cmd)
cmd="rm "+out_dir+"jay_c2_flats.dat"; run_shell_cmd(cmd)
cmd="rm "+out_dir+"jay_c2_spots.dat"; run_shell_cmd(cmd)
cmd="rm "+out_dir+"jay_residual_center_pixel_flat.dat"; run_shell_cmd(cmd)
cmd="rm "+out_dir+"jay_diff_fluxes_center_pixel.dat"; run_shell_cmd(cmd)
cmd="rm "+out_dir+"jay_ratio_fluxes_center_pixel.dat"; run_shell_cmd(cmd)
cmd="rm "+out_dir+"jay_median_flux_flats_center_pixel.dat"; run_shell_cmd(cmd)
cmd="rm "+out_dir+"jay_median_flux_spots_center_pixel.dat"; run_shell_cmd(cmd)
cmd="rm "+out_dir+"jay_median_deficit_flux_flats_center_pixel.dat"; run_shell_cmd(cmd)



for string in ["center", "n1", "n2", "n3", "n4"]:
    cmd="rm "+out_dir+"flat_calibration_"+string+".dat"; run_shell_cmd(cmd)




#data_flats=allFlats[0]  ## allRamps should have len == 1; this script is for one ramp. 
#data_spots=allSpots[0]  ## 

#last=int(sys.argv[1])
#data_flats=np.mean(allFlats[0:last], axis=0)
#allFlats=0
#data_spots=np.mean(allSpots[0:last], axis=0)
#allSpots=0


pp=PdfPages(out_dir+"nl_spots_ramps_stacked_100.pdf")



fig=plt.figure()
ax=fig.add_subplot(121)

ax=fig.add_subplot(121)
ax.imshow(GLOBAL_FLATS[-1], cmap=cm.seismic, origin="lower", interpolation='nearest')#, vmin=mean_science - mean_science/factor, vmax=mean_science + mean_science/factor)
ax.set_title("100 stacked images: Flats, last frame")
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
plt.colorbar(PCM, ax=ax)

ax=fig.add_subplot(122)
ax.imshow(GLOBAL_SPOTS[-1], cmap=cm.seismic, origin="lower", interpolation='nearest')#, vmin=mean_science - mean_science/factor, vmax=mean_science + mean_science/factor)
ax.set_title("100 stacked images: Spots, last frame")
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
plt.colorbar(PCM, ax=ax)

plt.tight_layout()
pp.savefig()


MASK=pf.open("/projector/aplazas/master-euclid-mask.fits")[0].data
fill_value=-999999

### Use only pixels for which x > 1000
#GLOBAL_RAMP=GLOBAL_RAMP[:,:,1000:]
#GLOBAL_RAMP2=GLOBAL_RAMP2[:,:,1000:]
#MASK=MASK[:,1000:]

npix_total=MASK.shape[0]*MASK.shape[1]
centroid_threshold=0.1
nl_threshold_f, nl_threshold_s= 3.0, 3.0

x_cut, y_cut = 10, 10
centroid_type = 'center'

if centroid_type == 'center':
    MIN_CENTROID_THRESHOLD, MAX_CENTROID_THRESHOLD = 0.001, 0.0 + centroid_threshold
elif centroid_type == 'corner':
   MIN_CENTROID_THRESHOLD, MAX_CENTROID_THRESHOLD = np.sqrt(2)*0.5 - centroid_threshold, np.sqrt(2)*0.5
else: 
    print "Enter a vaild centroid type: center or corner."
    sys.exit(1)

SEXTRACTOR="/usr/local/optical/sextractor/bin/sex"

#pp=PdfPages(out_dir+"nl_spots_ramps_stacked_%g.pdf" %ramp_number)

#### Step 1. Run sextractor on spots. 

gc.collect()

#last=GLOBAL_SPOTS[-1]/gain   ## Return to ADU to do the detecttion in SExtractor

cmd="rm "+out_dir+"/flat_calibration_*.dat"; run_shell_cmd(cmd)


last=data_spots[-1]/gain
print "last.shape: ", last.shape
cmd="rm science_andres.fits"; run_shell_cmd(cmd)


pf.writeto ("science_andres.fits", last, clobber=True)
prefix="hola"


"""
cmd="%s science_andres.fits, science_andres.fits -c daofind_sex_detection.config"%(SEXTRACTOR); run_shell_cmd (cmd)
#cmd="mkdir -v %s" %out_dir; run_shell_cmd(cmd)
out=out_dir + '/' + prefix + '_sextractor_out_last_sample_ramp_100.param'
cmd="mv check.fits %s/%s_check.fits" %(out_dir, prefix); run_shell_cmd(cmd)
cmd="mv output.cat %s" %(out); run_shell_cmd(cmd)
#cmd="rm science_andres.fits"; run_shell_cmd(cmd)
"""


out=out_dir + '/' + prefix + '_sextractor_out_last_sample_ramp_100.param'

if simulation == False:    
    ## REad position from detected objects catalog
    data=np.genfromtxt(out)
    flux_d=data[:,0]
    x_d=data[:,5]
    y_d=data[:,6]

    print len(x_d), len(y_d)
    positions_file='Source Extractor'
else:
    ### SIMULATIONS
    #positions_file="/projector/aplazas/MAY30/sextractor_positions_ppl_data_03mar17.txt" 
    #positions_file="/home/aplazas/simulations/sim_spots_position_random_stamp_size_16.txt"
    #positions_file="/home/aplazas/simulations/sim_spots_position_00_05.txt"
    #positions_file="/home/aplazas/simulations/sim_spots_position_00.txt"
    positions_file="/home/aplazas/simulations/sim_spots_position_random_norm_0_point_1.txt"
    data_HAND=np.genfromtxt(positions_file)
    x_d=data_HAND[:,0]
    y_d=data_HAND[:,1]

print len(x_d), len (y_d)
#sys.exit()

print "NUMBER OF DETECTED OBJECTS: ", len(x_d)

x_int, y_int = np.rint(x_d), np.rint(y_d)

#### Calculate centroid by hand from last image


LAST_TEST= GLOBAL_SPOTS[-1]
pf.writeto(out_dir+"LAST.fits", LAST_TEST, clobber=True)
#sys.exit()



base_pix=7
temp_hist2=[]

DETECTION_GLOBAL_SPOTS=GLOBAL_SPOTS - B_spots  # use this to calculate unweighted centroids 

print GLOBAL_SPOTS
print B_spots
#sys.exit()


x_int_filtered, y_int_filtered, flux_filtered, central_filtered =[],[],[],[] 
x_centroid_filtered, y_centroid_filtered=[],[]


#print x_int, y_int

if stamp_string == 'three':
    stamp_end = 1
elif stamp_string == 'five':
    stamp_end=2
else:
    print "Invalid stamp string"
    sys.exit(1)



for (xc, yc) in zip (x_int, y_int):
    if xc > 2043 or xc < 5: continue
    if yc > 2000 or yc < 5: continue
    xc-=1
    yc-=1

     
    ### To select centred and cornered spots, use the method of unweighted centroid: answers are; center: 0.0; corner: 0.707 (sqrt(2)/2).
    stamp=DETECTION_GLOBAL_SPOTS[-1][yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end] 
    ## TEMP
    #if stamp[2,2] > 6e4: continue    

    print " "  
    ## get local background through median 
    #print DETECTION_GLOBAL_SPOTS[-1][yc + 10:yc+30, xc - 20:xc + 20]
    bck = np.median (DETECTION_GLOBAL_SPOTS[-1][yc + 10:yc+30, xc - 20:xc + 20])
    #bck=0
    print "local background: ", bck
    
    if np.isnan(bck): 
    	continue

    stamp-=bck

    flux=np.sum(GLOBAL_SPOTS[-1][yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end])

    flux_central=GLOBAL_SPOTS[-1][yc,xc]

    #print np.max(stamp), np.sum(stamp), np.max(stamp)/np.sum(stamp)
    #temp_hist.append(np.max(stamp)/np.sum(stamp))
    print stamp
    if stamp_string == 'three':
        x_centroid, y_centroid = get_centroid_3x3 (stamp)    
    elif stamp_string == 'five':
        x_centroid, y_centroid = get_centroid_5x5 (stamp)
    else:
        print "invalid stamp string"
        sys.exit(1)
    norm_centroid=np.sqrt (x_centroid**2 + y_centroid**2)
    if np.isnan(norm_centroid) or (norm_centroid > 1.0): continue
    temp_hist2.append(norm_centroid)
    print "x centroid, y centroid, norm: ", x_centroid, y_centroid, norm_centroid
    
    if centroid_type == 'center':
        if (norm_centroid >= MIN_CENTROID_THRESHOLD) and (norm_centroid <= MAX_CENTROID_THRESHOLD):
            x_int_filtered.append(xc)
            y_int_filtered.append(yc)
            flux_filtered.append(flux)
            central_filtered.append(flux_central)
            x_centroid_filtered.append(x_centroid)
            y_centroid_filtered.append(y_centroid)
            print "HOLA", xc, yc, norm_centroid, flux

    if  centroid_type == 'corner':
        if (norm_centroid >= MIN_CENTROID_THRESHOLD) and (norm_centroid <= MAX_CENTROID_THRESHOLD) and (x_centroid > 0) and (y_centroid > 0):    #Region 4
        #if (norm_centroid >= MIN_CENTROID_THRESHOLD) and (norm_centroid <= MAX_CENTROID_THRESHOLD) and (x_centroid < 0) and (y_centroid > 0):    #Region 3
        #if (norm_centroid >= MIN_CENTROID_THRESHOLD) and (norm_centroid <= MAX_CENTROID_THRESHOLD) and (x_centroid > 0) and (y_centroid < 0):    #Region 2
        #if (norm_centroid >= MIN_CENTROID_THRESHOLD) and (norm_centroid <= MAX_CENTROID_THRESHOLD) and (x_centroid < 0) and (y_centroid < 0):     #Region 1
            x_int_filtered.append(xc)
            y_int_filtered.append(yc)
            flux_filtered.append(flux)
            central_filtered.append(flux_central)
            x_centroid_filtered.append(x_centroid)
            y_centroid_filtered.append(y_centroid)
            print "HOLA CORNER", xc, yc, norm_centroid, flux  

## TEMP
#x_int_filtered=x_int
#y_int_filtered=y_int

x_int_filtered=np.array(x_int_filtered)
y_int_filtered=np.array(y_int_filtered)

##### print the positions of teh selected spots (centroid < 0.1)

f=open(out_dir+"selected_positions_centroid.dat", 'w')
for a, b in zip(x_int_filtered, y_int_filtered):
    line="%g %g \n" %(a,b)
    f.write(line)
f.close()


## Do a scatter plot of centroids 
fig=plt.figure()
ax=fig.add_subplot(111)
ax.scatter (x_centroid_filtered, y_centroid_filtered, color='r', s=10, marker='.', alpha=0.3)
ax.set_xlabel("x centroid")
ax.set_ylabel("y centroid")
pp.savefig()


end=-1

temp_hist2=np.array(temp_hist2)
m,dev,indixec = sigma_clip.sigma_clip(temp_hist2,niter=6, nsig=sigma_cut, get_indices=True, verbose=False)

fig = plt.figure()
ax = fig.add_subplot(211)
n, bins, patches_out = ax.hist(temp_hist2, 1000, normed=False, color='green', histtype='step', alpha=0.75 , label='N, mean (after sc), std (after sc): \n  %g, %g, %g' %(len(temp_hist2), m,dev ) )
#ax.set_title('norm of unweighted centroid in stamp \n %s' %(positions_file), size=8)
ax.legend(loc=loc_label , fancybox=True, ncol=1, numpoints=1, prop = prop)
ax.tick_params(axis='both', which='major', labelsize=7)
ax.set_yscale('log')
ax.set_xlabel('norm of centroid vector')
ax.set_ylabel('Counts')
#plt.xlim([-0.2,1])

central_filtered=np.array(central_filtered)
flux_filtered=np.array(flux_filtered)
m,dev,indixec = sigma_clip.sigma_clip(central_filtered,niter=6, nsig=sigma_cut, get_indices=True, verbose=False)

ax = fig.add_subplot(212)
n, bins, patches_out = ax.hist(central_filtered, 1000, normed=False, color='blue', histtype='step', alpha=0.75 , label='N, mean (after sc), std (after sc): \n  %g, %g, %g' %(len(central_filtered), m,dev ) )
#ax.set_title('Flux of stamp (back subtracted) \n %s' %(positions_file), size=8)
ax.legend(loc=loc_label , fancybox=True, ncol=1, numpoints=1, prop = prop)
ax.tick_params(axis='both', which='major', labelsize=7)
ax.set_xlabel('Flux in central pixel (electrons, from last sample; back. subtracted)')
ax.set_ylabel('Counts')

#ax.set_yscale('log')

pp.savefig()
#pp.close()

#sys.exit()

#For real data, seems to be 1e5-2.3e5
min_flux=1.5e5 #8e4
max_flux=4.0e5 #2.5e5

#### Cuts for binning data by flux measured in sample with longest time
cut1, cut2, cut3, cut4 = min_flux, 2.75e5, max_flux, max_flux
index= (cut1 <= flux_filtered) & (flux_filtered < cut2)
flux_fbin1=flux_filtered[index]
x_int_fbin1=x_int_filtered[index]
y_int_fbin1=y_int_filtered[index]

index= (cut2 <= flux_filtered) & (flux_filtered < cut3)
flux_fbin2=flux_filtered[index]
x_int_fbin2=x_int_filtered[index]
y_int_fbin2=y_int_filtered[index]
#r_fbin2=r[index]

#index= (cut3 <= flux_filtered) & (flux_filtered <= cut4)
#flux_fbin3=flux_filtered[index]
#x_int_fbin3=x_int_filtered[index]
#y_int_fbin3=y_int_filtered[index]
#r_fbin3=r[index]

#print "Len(bin1), len(bin2), len(bin3): ", len (flux_fbin1), len(flux_fbin2), len(flux_fbin3)
print "len (all filtered): ", len(flux_filtered)

#"""



#sys.exit()

### Step 2: 
### Get NL coefficients for pixels with flats
### Apply correction to 9 pixels within spot. 

def plot_pixel_ramp(ramp='', time='', fig='',i=0, j=0, counter='', fmt='k-o', plot_flag=False):
    if not len (ramp) == len (time):
	print "inside function 'plot_pixel_ramp': len(ramp) not equal to len (time)."
        print "len(ramp), len(time): ", len(ramp), len(time)
        sys.exit()

    time_vec, signal_vec, signal_vec_err=[],[],[]
    a,b=[],[]
    counter2=0
    for t, sample_array in zip(time[:],ramp[:]):
        
        index_x, index_y =  j, i #np.unravel_index (counter, sample_array.shape)
        
        s=sample_array[index_x, index_y]
        time_vec.append(t)
        #signal_vec.append( 2**16-1-s )
        signal_vec.append(s) # s=ADU_dark - ADU_data
        signal_vec_err.append(np.sqrt(s))
        if not counter2 == 0:
            a.append(time_vec[counter2]-time_vec[counter2-1])
            b.append(signal_vec[counter2]-signal_vec[counter2-1])
        counter2+=1

    
    time_vec=np.array(time_vec)
    signal_vec=np.array(signal_vec)
    signal_vec_err=np.array(signal_vec_err)
    
    if plot_flag: 
        ax=fig.add_subplot(3,3,counter)
        ax.errorbar(time_vec, signal_vec, yerr=None, fmt=fmt, alpha=1.0, label="(%g, %g, %g)" %(i,j,counter), markersize=5) 
    #ax.errorbar(time_vec[1:],b, yerr=None, fmt=fmt, alpha=1.0, label="(%g, %g, %g)" %(i,j,counter), markersize=5)         

        ax.legend(loc='upper left', fancybox=True, ncol=1, numpoints=1, prop = prop)
        #ax.set_xlabel ('Time (seconds)')
        #ax.set_ylabel ('Signal (e-)')
        ax.tick_params(axis='both', which='major', labelsize=5)
        ax.tick_params(axis='both', which='minor', labelsize=5)    
    
    return np.array(b)/np.array(a)
    #return signal_vec

#initialize dictionaries that will hold ramps for each of the 25 pixels
#start_x, end_x=1,1
#start_y, end_y=1,1

### For 9 pixels

if stamp_string == 'three':
    #stamp_range_global_x, stamp_range_global_y, stamp_range_local = [xc-1,xc,xc+1],[yc-1,yc,yc+1], [0,1,2]
    STAMP_SIZE=3
    get_pixel_index={(0,0):7, (0,1):4, (0,2):1, (1,0):8, (1,1):5, (1,2): 2, (2,0):9, (2,1):6, (2,2):3 }
    get_pair_index={7:(0,0), 4:(0,1), 1:(0,2), 8:(1,0), 5:(1,1), 2:(1,2), 9:(2,0), 6:(2,1), 3:(2,2) }
elif stamp_string == 'five': 
    #stamp_range_global_x, stamp_range_global_y, stamp_range_local=[ xc-2,xc-1,xc,xc+1,xc+2],[ yc-2,yc-1,yc,yc+1,yc+2], [0,1,2,3,4]
    STAMP_SIZE=5
    #for 25 pixels
    get_pixel_index={(0,4):1, (1,4):2, (2,4):3, (3,4):4, (4,4):5, \
                 (0,3):6, (1,3):7, (2,3):8, (3,3):9, (4,3):10, \
                 (0,2):11, (1,2):12, (2,2):13, (3,2):14, (4,2):15,\
                 (0,1):16, (1,1):17, (2,1):18, (3,1):19, (4,1):20,\
                 (0,0):21, (1,0):22, (2,0):23, (3,0):24, (4,0):25}

    get_pair_index={ 1:(0,4),   2:(1,4),   3:(2,4), 4:(3,4),   5:(4,4), \
                 6:(0,3),  7:(1,3), 8:(2,3), 9:(3,3), 10:(4,3), \
                 11:(0,2), 12:(1,2), 13:(2,2), 14:(3,2), 15:(4,2),\
                 16:(0,1), 17:(1,1), 18:(2,1), 19:(3,1), 20:(4,1),\
                 21:(0,0), 22:(1,0), 23:(2,0), 24:(3,0), 25:(4,0)}
else:
    print "invalid string stamp"
    sys.exit(1)



def initialize_ramps_dict (ramps_dict):
    for j in range(STAMP_SIZE):
        for i in range(STAMP_SIZE):
            ramps_dict[(j,i)]=[]

ramps_dict_all, ramps_dict_all_jay={},{}
ramps_dict_fbin1={}
ramps_dict_fbin2={}
ramps_dict_fbin3={}

ramps_dict_all_jay_mean_signal={}


initialize_ramps_dict (ramps_dict_all_jay_mean_signal)
initialize_ramps_dict (ramps_dict_all)
initialize_ramps_dict (ramps_dict_all_jay)
initialize_ramps_dict (ramps_dict_fbin1)
initialize_ramps_dict (ramps_dict_fbin2)
initialize_ramps_dict (ramps_dict_fbin3)

#start_x, end_x=1,1
#start_y, end_y=1,1

print len(x_int), len(y_int)
super_counter=1

max_vec=[]
sigma_x_vec, sigma_y_vec=[],[]


#dict_c2_flat={"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[],"10":[],\
             #"11":[],"12":[],"13":[],"14":[],"15":[],"16":[],"17":[],"18":[],"19":[],"20":[],\
             #"21":[],"22":[],"23":[],"24":[],"25":[]}
#dict_c2_spot={"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[],"10":[],\
             #"11":[],"12":[],"13":[],"14":[],"15":[],"16":[],"17":[],"18":[],"19":[],"20":[],\
             #"21":[],"22":[],"23":[],"24":[],"25":[]}
#dict_c2_dark={"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[],"10":[],\
             #"11":[],"12":[],"13":[],"14":[],"15":[],"16":[],"17":[],"18":[],"19":[],"20":[],\
             #"21":[],"22":[],"23":[],"24":[],"25":[]}

#dict_c2_flat_corr={"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[],"10":[],\
             #"11":[],"12":[],"13":[],"14":[],"15":[],"16":[],"17":[],"18":[],"19":[],"20":[],\
             #"21":[],"22":[],"23":[],"24":[],"25":[]}
#dict_c2_spot_corr={"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[],"10":[],\
             #"11":[],"12":[],"13":[],"14":[],"15":[],"16":[],"17":[],"18":[],"19":[],"20":[],\
             #"21":[],"22":[],"23":[],"24":[],"25":[]}
#dict_c2_dark_corr={"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[],"10":[],\
             #"11":[],"12":[],"13":[],"14":[],"15":[],"16":[],"17":[],"18":[],"19":[],"20":[],\
             #"21":[],"22":[],"23":[],"24":[],"25":[]}


final=0

size_final=[]
central_signal_size_final=[]
#siga=False
#siga2=False

def plot_near_nl_corr (fig, frame, p_spot_corr, signal_corrected_spots, p_flat_corr, signal_corrected_flats, pos_x, pos_y):
    linear_s = p_spot_corr[1]*time_darks
    #linear_s=[12889.7802830, 25779.5605659, 38669.3408489, 51559.1211319, 64448.9014149]
    #linear_s=[ 6444.8901415, 12889.7802830, 19334.6704245, 25779.5605659, 32224.4507074]
    #linear_s=[ 19334.6704245, 38669.3408489, 58004.0112734, 77338.6816978, 96673.3521223 ]
    res_s=100*(signal_corrected_spots-linear_s)/linear_s
  
    print "linear_s: ", linear_s
    print "signal_corrected_spots ", signal_corrected_spots
    
    linear_f = p_flat_corr[1]*time_flats
    #linear_f=[18000.0000000, 36000.0000000, 54000.0000000, 72000.0000000, 90000.0000000]
    #linear_f=[9000.0000000, 18000.0000000, 27000.0000000, 36000.0000000, 45000.0000000]
    #linear_f=[27000.0000000, 54000.0000000, 81000.0000000, 108000.0000000, 135000.0000000] 
    res_f=100*(signal_corrected_flats-linear_f)/linear_f
                    
    print "linear_f: ", linear_f
    print "signal_corrected_flats ", signal_corrected_flats

    flag=False
    if  (np.max(np.abs(res_f))) > nl_threshold_f or (np.max(np.abs(res_s))) > nl_threshold_s:
        flag=True
        return flag 

    ax=fig.add_subplot(3,3, frame)
    plt.plot(time_darks, res_s, 'r.-', label='spots')
    plt.plot(time_flats, res_f, 'b.-', label='flats')
    ax.legend(loc=loc_label , fancybox=True, ncol=1, numpoints=1, prop = prop)
    
    ax.set_ylabel('Fractional NL\n(after NL correction)', size=8)
    ax.set_xlabel('Time (mili seconds)', size=8)
    ax.set_title("%g electrons.\n Position: (%g,%g)" %(signal_corrected_spots[-1],pos_x,pos_y), size=8)
    plt.tight_layout() 

    print "plot_near_nl_corr residuals spots in %: ", res_s
    print "plot_near_nl_corr residuals flats in %: ", res_f 
    print  "plot_near_nl_corr: ", frame, p_spot_corr, signal_corrected_spots, p_flat_corr, signal_corrected_flats, pos_x, pos_y 

    ### PRINT 
    if frame == 5: string="center"
    if frame == 2: string="n1"
    if frame == 4: string="n2"
    if frame == 6: string="n3"
    if frame == 8: string="n4"
    f=open(out_dir+"flat_calibration_"+string+".dat", 'a+')
    #line_array=np.hstack ((linear_s, res_s, linear_f, res_f))
    #line=""
    #for x in line_array:
    #    line+="%g "%x
    #line+="\n"
    if flag == False: 
        for (a,b,c,d) in zip(linear_s, res_s, linear_f, res_f):
            line="%g %g %g %g \n"%(a,b,c,d)
            f.write(line)
        f.close()
    return flag


c0_diff_spots_vec=[]
c0_diff_flats_vec=[]
c0_diff_darks_vec=[]

to_plot=[1,2,3,4,5,10,15,20,25,30,35,40,45,50,51,52,60,63,68,70,72,73,75,79]
#to_plot=[15, 20]


to_plot=[1,2,3,4,5,10,15,20,25,30,35,40,45,50,51,52,60,63,68,70,72,73,75,79, 74,115,180,181,206,226,228,274,306,317,318,321,340,344,345]


print np.min(flux_filtered), np.max(flux_filtered)
#sys.exit()


NORM_flats_big_vec=[]
NORM_spots_big_vec=[]


big_vec_c2_flats=[]
big_vec_c2_spots=[]
residual_vec=[]
ratio_fluxes_vec=[]
diff_fluxes_vec=[]


signal_flat_SUPER_VEC=[]
signal_spot_SUPER_VEC=[]

signal_deficit_central_pixel_spots=[]
signal_deficit_central_pixel_flats=[]


dict_residual={(0,0):[],(0,1):[],(0,2):[],(1,0):[],(1,1):[],(1,2):[],(2,0):[],(2,1):[],(2,2):[]}
dict_median_flux_flats={(0,0):[],(0,1):[],(0,2):[],(1,0):[],(1,1):[],(1,2):[],(2,0):[],(2,1):[],(2,2):[]}
dict_median_flux_spots={(0,0):[],(0,1):[],(0,2):[],(1,0):[],(1,1):[],(1,2):[],(2,0):[],(2,1):[],(2,2):[]}




for L, (xc, yc, f, central) in enumerate( zip (x_int_filtered[:end], y_int_filtered[:end], flux_filtered[:end], central_filtered[:end]) ):
    

    stamp=GLOBAL_SPOTS[-1][yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end]
    stamp_mask = MASK [yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end]
    if not np.sum(stamp_mask) == 0:
	print "Discarding spot because it has at least one bad pixel"
        continue



    if L in to_plot:
        plot_flag=True
        st=stamp
    else:
        plot_flag=False   

 
    #if central < 1.2e5 and central > 2.1e5: 
    #    continue    
  

    #print "L: ", L 
    #continue

 
    siga=False
    siga2=False
    #xc-=1
    #yc-=1
 

    print " "
    print "OBJECT AT: (%g,%g)" %(xc, yc)
    
    max=np.max(stamp)    
    max_vec.append(max)
    #### Fit a 2D Gaussian to stamp
    #x,y = np.array([0,1,2]), np.array([0,1,2])
    #x, y = np.meshgrid (x,y)
    #initial_guess=(max,1.0,1.0,2.0,2.0,0,1e2)
    #popt, pcov = optimize.curve_fit(twoD_Gaussian, (x, y), stamp.ravel(), p0=initial_guess, maxfev=100000000)
    #print "GAUSSIAN parameters and errors: ", popt, np.diag(pcov)
    #sx, sy = popt[3], popt[4]
    #sigma_x_vec.append(sx)
    #sigma_y_vec.append(sy)

    #continue

    print "stamp shape:: ", stamp.shape
    
    c0_mat_spots=0*stamp
    c0_mat_darks=0*stamp
    c0_mat_flats=0*stamp

    c2_mat_flats=-7.76e-7*np.ones_like(stamp)
    c2_mat_spots=0*stamp
    c2_mat_spots_corr=0*stamp
 
    c2_mat_flats_err=0*stamp
    c2_mat_spots_err=0*stamp
    c2_mat_spots_corr_err=0*stamp

    c1_mat_spots=0*stamp

    counter=1
    fig2=plt.figure()
    fig3=plt.figure()
 

    #if True:
    #if L in to_plot:
    #    plot_flag=True
    #    st=stamp
    #else:
    #    plot_flag=False


    counter=1
    #siga=False
    #TEMP_FLAG=False
    #siga2=False

    corrected_signal_dict={}
    sample_corrected_ramp_per_pixel={1:[],2:[],3:[],4:[],5:[]}
    initialize_ramps_dict (corrected_signal_dict)

    if L in to_plot:
    #if True:
	fig_near=plt.figure()



    ### Check NL in central pixel
    #stamp_spots=GLOBAL_SPOTS[:, yc-2:yc+3, xc-2:xc+3]
    #stamp_flats=GLOBAL_FLATS[:, yc-2:yc+3, xc-2:xc+3]
    #stamp_darks=darks[:, yc-2:yc+3, xc-2:xc+3]


    #p_spot, cov_spot, chi2_spot,flag2 =fit_pixel_ramp (ramp=stamp_spots,time=time, i=2, j=2, pinit=pinit)
    #p_flat, cov_flat, chi2_flat,flag1 =fit_pixel_ramp (ramp=stamp_flats,time=time, i=2, j=2, pinit=pinit)
    #p_dark, cov_dark, chi2_dark,flag3 =fit_pixel_ramp (ramp=stamp_darks,time=time, i=2, j=2, pinit=pinit)

    #p_spot_corr, cov_spot_corr, chi2_spot_corr, flag4, flux_rate_spots, delta_signal_spots, delta_sums_spots, signal_corrected_spots = correct_and_fit_pixel_ramp (ramp=stamp_spots, dark_ramp=stamp_darks, time=time, i=2, j=2, pinit=pinit, c0_spot=p_spot[0], c0_dark=p_dark[0], c1_spot=p_spot[1], c1_dark=p_dark[1], c2=p_flat[2])


    #linear_s = p_spot_corr[0] + p_spot_corr[1]*time
    #res_s=100*(signal_corrected_spots-linear_s)/linear_s

    #a=np.fabs(np.max(res_s))
    #if a > 5: continue
    #####
    if stamp_string == 'three':
        stamp_range_global_x, stamp_range_global_y, stamp_range_local=[xc-1,xc,xc+1],[yc-1,yc,yc+1], [0,1,2]
    elif stamp_string == 'five':
        stamp_range_global_x, stamp_range_global_y, stamp_range_local=[xc-2,xc-1,xc,xc+1,xc+2],[ yc-2,yc-1,yc,yc+1,yc+2], [0,1,2,3,4]    
    else:
        print "invalid stamp_string"
        sys.exit(1)

    dict_B={(1,1):[[],[]], (1,0):[[],[]], (0,1):[[],[]], (2,1):[[],[]], (1,2):[[],[]]}


 
    for i,k in zip(stamp_range_global_x, stamp_range_local):
        #if i == 0 or i == 4: continue
        for j,l in zip(stamp_range_global_y, stamp_range_local):
            #if j == 0 or j == 4: continue
            i,j=int(i), int(j)
            #i,j,k,l=xc,yc,1,1
            stamp_spots=GLOBAL_SPOTS[:, yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end]
            stamp_flats=GLOBAL_FLATS[:, yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end]
            stamp_darks=darks[:, yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end]


            print "Last sample of stamp spots: ", stamp_spots[-1]
            #sys.exit(1)


            #plot_pixel_ramp(ramp=stamp_flats, time=time, i=k,j=l, fig=fig2, counter=counter, plot_flag=False)
            #plot_pixel_ramp(ramp=stamp_spots, time=time, fig=fig3, i=k,j=l, counter=counter, fmt='b-o', plot_flag=False)
            
            p_flat, cov_flat, chi2_flat,flag1, signal_flat =fit_pixel_ramp (ramp=stamp_flats,time=time_flats, i=k, j=l, pinit=pinit)
            if flag1 == True:
                #siga=True
                print "flat fit"
                sys.exit()
                #break
     
 
            print p_flat[0]
            time_darks=np.array(time_darks) 
            p_spot, cov_spot, chi2_spot,flag2, signal_spot =fit_pixel_ramp (ramp=stamp_spots,time=time_darks, i=k, j=l, pinit=pinit)         


            if flag2 == True:
                #siga=True
                print "spot fit"
                sys.exit()
                #break
            print p_spot[0]           

            
            print time_flats, p_flat
            print time_darks, p_spot
 
            p_dark, cov_dark, chi2_dark,flag3, signal =fit_pixel_ramp (ramp=stamp_darks,time=time_darks, i=k, j=l, pinit=pinit)
            if flag3 == True: 
                #siga=True
                print "dark fit"
                sys.exit()
                #break
            
            #if i == xc and j == yc:
            ### To test NL model, print data - model / model for centra pixel in flats
            model_flat = fitfunc (time_flats, p_flat[0], p_flat[1], p_flat[2])
            if not len (model_flat) == len (signal_flat):
                    print "Lengths of model and signal in flat not the same"
                    sys.exit(1)
            residual=100*(signal_flat - model_flat)/ model_flat
            print "MODEL RESIDUAL IN  PIXEL: ", k,l,  residual
            #abs_res=np.max(np.abs(residual))
                #if abs_res > 0.2:
                #   FLAG_RESIDUAL_FLAT_CENTRAL_PIXEL=True
                #else:
                   #FLAG_RESIDUAL_FLAT_CENTRAL_PIXEL=False
            dict_residual[(k,l)].append(residual)
                
            dict_median_flux_flats[(k,l)].append(signal_flat)
            #print "ANDRES:::: ", dict_median_flux_flats[(k,l)]
            dict_median_flux_spots[(k,l)].append(signal_spot)

            print "p_flat[0], p_dark[0], p_spot[0]: ", p_flat[0], p_dark[0], p_spot[0]
            #p_flat[0]=0.0
            #p_dark[0]=0.0
            #p_spot[0]=0.0
            print "p_flat[0], p_dark[0], p_spot[0]: ", p_flat[0], p_dark[0], p_spot[0]
            #if i == xc and j == yc:
            #    sys.exit() 


            #if p_spot[0] < 0: 
		#siga2=True
                #print "p_spot[0]< 0 ", p_spot
                #p_spot=np.fabs(p_spot)
                #sys.exit()
		#break

            #print time, stamp_darks
           

            #print "p_dark, cov_dark, chi2_dark,flag3: ", p_dark, cov_dark, chi2_dark,flag3 
            #sys.exit(1)
            p_spot_corr, cov_spot_corr, chi2_spot_corr, flag4, flux_rate_spots, delta_signal_spots, delta_sums_spots, signal_corrected_spots = correct_and_fit_pixel_ramp (ramp=stamp_spots, dark_ramp=stamp_darks, time=time_darks, i=k, j=l, pinit=pinit, c0_spot=p_spot[0], c0_dark=p_dark[0], c1_spot=p_spot[1], c1_dark=p_dark[1], c2=p_spot[2])
            
            if (k,l) in [(1,1), (1,0), (0,1), (1,2), (2,1)]:
                dict_B[(k,l)][0]=p_spot_corr
                dict_B[(k,l)][1]=signal_corrected_spots 
             
 
            #if flag4 == True: continue
            #print "xc, xy, i, j, p flats , pe flats, p spots, pe spots: ", xc, yc, i, j, p_flat, np.diag(cov_flat), p_spot, np.diag(cov_spot)
            
            #corrected_signal_dict[(l,k)].append(signal_corrected_spots)
            #for sample in range(1,len(signal_corrected_spots)+1):
		#sample_corrected_ramp_per_pixel[sample].append(signal_corrected_spots[sample-1])
            
            print "Last sampe of stamp flats: ", stamp_flats[-1]

            p_flat_corr, cov_flat_corr, chi2_flat_corr, flag4, flux_rate_flats, delta_signal_flats, delta_sum_flats, signal_corrected_flats = correct_and_fit_pixel_ramp (ramp=stamp_flats, dark_ramp=stamp_darks, time=time_flats, i=k, j=l, pinit=pinit, c0_spot=p_flat[0], c0_dark=p_dark[0], c1_spot=p_flat[1], c1_dark=p_dark[1], c2=p_flat[2])

            #if l == 1 and k == 2:
            #    print "Spots: ", p_spot, p_spot_corr
            #    print "Flats: ", p_flat, p_flat_corr
            #    sys.exit()

  
            if L in to_plot:
            #if True: 
                # Check for linearity in nearest four neighbours 
                if (i == xc-1) and (j == yc): 
                    use1=plot_near_nl_corr (fig_near, 4, p_spot_corr, signal_corrected_spots, p_flat_corr, signal_corrected_flats, xc-1, yc) 
                if (i == xc+1) and (j == yc):
                    use2=plot_near_nl_corr (fig_near, 6, p_spot_corr, signal_corrected_spots, p_flat_corr, signal_corrected_flats, xc+1, yc)
                if (i == xc) and (j == yc+1):
                    use3=plot_near_nl_corr (fig_near, 2, p_spot_corr, signal_corrected_spots, p_flat_corr, signal_corrected_flats, xc, yc+1)
                
                if (i == xc) and (j == yc-1):
                    use4=plot_near_nl_corr (fig_near, 8, p_spot_corr, signal_corrected_spots, p_flat_corr, signal_corrected_flats, xc, yc-1)
                
                if (i == xc) and (j == yc):
                    print "XC, YC: ", xc, yc 
                    use5=plot_near_nl_corr (fig_near, 5, p_spot_corr, signal_corrected_spots, p_flat_corr, signal_corrected_flats, xc, yc)
   

            print " "
            print "p flats , pe flats: ", p_flat, np.diag(cov_flat)
            print "p spots, pe spots: ", p_spot, np.diag(cov_spot)
            print "p dark, pe darks: ", p_dark, np.diag(cov_dark)
            print "p spots corr, pe spots corr:  ", p_spot_corr, np.diag(cov_spot_corr)
            print "flag1, flag2, flag3: ", flag1, flag2, flag3
            print "counter, i, j, xc, yc: ", counter, i, j, xc, yc
            print " "


            c2_mat_flats[l,k], c2_mat_flats_err[l,k]=p_flat[2], np.diag(cov_flat)[2] 
            c2_mat_spots[l,k], c2_mat_spots_err[l,k]=p_spot[2], np.diag(cov_spot)[2]
            c2_mat_spots_corr[l,k], c2_mat_spots_corr_err[l,k] =p_spot_corr[2], np.diag(cov_spot_corr)[2]
            c0_mat_spots[l,k]=p_spot[0]
            c1_mat_spots[l,k]=p_spot[1]           
            c0_mat_flats[l,k]=p_flat[0]
            c0_mat_darks[l,k]=p_dark[0]


            """
            diff_last_first = GLOBAL_SPOTS[-1][yc-2:yc+3, xc-2:xc+3] - GLOBAL_SPOTS[0][yc-2:yc+3, xc-2:xc+3]
            delta_t=time[-1] - time[-2]
            NORM = np.sum(diff_last_first)/ (  (len(GLOBAL_SPOTS)-1)   *delta_t)            

            s_vec=[]
            val0=flux_rate_spots[0]
            for val in flux_rate_spots:
                val=(val-val0)/NORM
                print "VAL: ", val
                s_vec.append(val)
            """
       
            counter+=1


    if L in to_plot:
        #plt.suptitle("Fractional NL after correction")
        pp.savefig(fig_near)
    #if correct_NL == True:
    #    if (use1 and use2 and use3 and use4 and use5): continue  ### If the NL correction is not < than 2%
    #if correct_NL == True:
        #if FLAG_RESIDUAL_FLAT_CENTRAL_PIXEL==True: continue 
        
    # Save fitted C2 coefficients from flat-field images 
    for i in c2_mat_flats.flatten(): big_vec_c2_flats.append(i)
    #for i in c2_mat_spots.flatten(): big_vec_c2_spots.append(i)
    big_vec_c2_spots.append(c2_mat_spots[1,1])


    # Produce a NL-corrected stamp
    GLOBAL_SPOTS_STAMPS_CORR, corr_stamps_sum=[],[]
    for t, stamp, dark_stamp in zip(time_darks, GLOBAL_SPOTS[:, yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end], darks [:,yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end]):
	print "stamp: ", stamp
	print " "
   	print "dark_stamp: ", dark_stamp
        print " "
        print "c2_mat_flats: ", c2_mat_flats
        print " "
        print "c0_mat_spots: ", c0_mat_spots
        print " "
        print "c0_mat_darks: ", c0_mat_darks
        print " "
        corr_stamp_spots = (1/(2*c2_mat_flats))*(-1+np.sqrt(1-4*c2_mat_flats*(c0_mat_spots-stamp)) )
        corr_stamp_darks = (1/(2*c2_mat_flats))*(-1+np.sqrt(1-4*c2_mat_flats*(c0_mat_darks-dark_stamp)) )
        print "corr_stamp_darks: ", corr_stamp_darks
        print " "
        
        if correct_NL == False:
            corr_stamp = stamp - dark_stamp
        else:
            corr_stamp = corr_stamp_spots - corr_stamp_darks

        print "corr_stamp_spots: ", corr_stamp_spots
	#print "corr_stamp_darks: ", corr_stamp_darks
        print "corr_stamp: ", corr_stamp


        GLOBAL_SPOTS_STAMPS_CORR.append (corr_stamp)
        corr_stamps_sum.append(np.sum(corr_stamp))
   


    #### CORRECT STAMP OF FLATS, ALSO SAVE CORRECTED DARKS
    GLOBAL_FLATS_STAMPS_CORR, corr_stamps_flats_sum=[],[]
    GLOBAL_DARKS_STAMPS_CORR, corr_stamps_darks_sum=[],[]

    for t, stamp, dark_stamp in zip(time_flats, GLOBAL_FLATS[:, yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end], darks [:,yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end]):
        print "FLAT stamp: ", stamp
        print " "
        print "dark_stamp: ", dark_stamp
        print " "
        print "c2_mat_flats: ", c2_mat_flats
        print " "
        print "c0_mat_spots: ", c0_mat_spots
        print " "
        print "c0_mat_darks: ", c0_mat_darks
        print " "
        corr_stamp_flats = (1/(2*c2_mat_flats))*(-1+np.sqrt(1-4*c2_mat_flats*(c0_mat_flats-stamp)) )
        corr_stamp_darks = (1/(2*c2_mat_flats))*(-1+np.sqrt(1-4*c2_mat_flats*(c0_mat_darks-dark_stamp)) )
        print "corr_stamp_darks: ", corr_stamp_darks
        print " "

        if correct_NL == False:
            corr_stamp = stamp - dark_stamp
        else:
            corr_stamp = corr_stamp_flats #- corr_stamp_darks

        print "corr_stamp_flats: ", corr_stamp_flats
        #print "corr_stamp_darks: ", corr_stamp_darks
        print "corr_stamp: ", corr_stamp

        GLOBAL_FLATS_STAMPS_CORR.append (corr_stamp)
        corr_stamps_flats_sum.append(np.sum(corr_stamp))

        GLOBAL_DARKS_STAMPS_CORR.append(corr_stamp_darks) 
        corr_stamps_darks_sum.append(np.sum(corr_stamp_darks)) 



    ##### Calculate the size of the corrected stamp
    temp_size, temp_central_signal=[],[]
    for i in range(0,len(GLOBAL_SPOTS_STAMPS_CORR)-1):
        print i
        #bck_f = np.median (DETECTION_GLOBAL_SPOTS[i+1][yc + 10:yc+30, xc - 20:xc + 20]) 
        #bck_i = np.median (DETECTION_GLOBAL_SPOTS[i][yc + 10:yc+30, xc - 20:xc + 20]) 
        stamp_f = GLOBAL_SPOTS_STAMPS_CORR[i+1]#[yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end] #- bck_f
        stamp_i = GLOBAL_SPOTS_STAMPS_CORR[i]#[#yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end] #- bck_i
        diff=stamp_f-stamp_i
        ref_e1, ref_e2, ref_s, centroid_x, centroid_y, dxc, dyc = moments.measure_e1e2R2(diff, skysub=False, sigma=1)
        print i,  ref_e1, ref_e2, ref_s, centroid_x, centroid_y
        temp_size.append(ref_s)
        temp_central_signal.append(  (stamp_f[1,1] + stamp_i[1,1])*0.5  )

    size_final.append(temp_size)
    central_signal_size_final.append( temp_central_signal )


    #### End of calculating size of corrected stamp 
 
    #sys.exit()
    #### Calculate Chaz's metric and Jay's metric 
    siga=False
    for i,k in zip(stamp_range_global_x, stamp_range_local):
        
        for j,l in zip(stamp_range_global_y, stamp_range_local):
            sig=[]
            for sample_s in GLOBAL_SPOTS_STAMPS_CORR:
            #for sample_s in GLOBAL_FLATS_STAMPS_CORR:
                sig.append(sample_s[l,k])

            delta_sig, delta_time, delta_sum =[],[],[]
            charge_previous_pixel = []
        
            delta_mean_signal=[]

            for counter in range(len(sig)-1):
               delta_sig.append(sig[counter+1]-sig[counter])
               delta_mean_signal.append( (sig[counter+1]+sig[counter])*0.5   )

               delta_sum.append(corr_stamps_sum[counter+1]-corr_stamps_sum[counter])
               delta_time.append(time_darks[counter+1]-time_darks[counter])
               charge_previous_pixel.append(sig[counter])

            delta_sig=np.array(delta_sig)
            delta_sum=np.array(delta_sum)
            delta_time=np.array(delta_time) 
            charge_previous_pixel=np.array(charge_previous_pixel)
            delta_mean_signal=np.array(delta_mean_signal)

            print delta_mean_signal
            #sys.exit()

            s_vec=delta_sig/delta_sum
            rates_vec_jay=delta_sig/delta_time ## Rates
	    print "JAY's METRIC: delta_sig, delta_time: ratio is rates_vec_jay ", delta_sig, delta_time, rates_vec_jay 
            #if k == 1 and l == 1:
            #    print p_spot, p_flat
            #    print "central pixel, rates jay / charge previous pixel: ", rates_vec_jay/charge_previous_pixel
            #    print "charge_previous_pixel: ", charge_previous_pixel
            #   print "signal: ", sig
            #    sys.exit() 

            ### Second 'derivative' for Jay's metric
	    #diff_last_first = GLOBAL_SPOTS_STAMPS_CORR[-1] - GLOBAL_SPOTS_STAMPS_CORR[0]
            #delta_t=time_darks[-1] - time_darks[-2]
           
            if diff_spots == False: 
                diff_last_first = GLOBAL_FLATS_STAMPS_CORR[-1] - GLOBAL_FLATS_STAMPS_CORR[0]
                delta_t=time_flats[-1] - time_flats[-2]
                NORM = np.sum(diff_last_first)/ (  (len(GLOBAL_FLATS_STAMPS_CORR) -1 ) *delta_t)   # M-2 = 4
            else:
                diff_last_first = GLOBAL_SPOTS_STAMPS_CORR[-1] - GLOBAL_SPOTS_STAMPS_CORR[0]
                delta_t=time_darks[-1] - time_darks[-2]
                NORM = np.sum(diff_last_first)/ (  (len(GLOBAL_SPOTS_STAMPS_CORR)-1)*delta_t)

            print "NORM:", NORM 
            print "len(GLOBAL_SPOTS_STAMPS_CORR): ", len(GLOBAL_SPOTS_STAMPS_CORR)
            print "delta_t: ", delta_t
            #sys.exit()
 
            ##### Save norm of flats and spots at the same time 
            diff_last_first_flats = GLOBAL_FLATS_STAMPS_CORR[-1] - GLOBAL_FLATS_STAMPS_CORR[0]
            delta_t_flats=time_flats[-1] - time_flats[-2]
            NORM_flats = np.sum(diff_last_first_flats)/ (  ( len(GLOBAL_FLATS_STAMPS_CORR) -1) *delta_t_flats)


            diff_last_first = GLOBAL_SPOTS_STAMPS_CORR[-1] - GLOBAL_SPOTS_STAMPS_CORR[0]
            delta_t=time_darks[-1] - time_darks[-2]
            NORM_spots = np.sum(diff_last_first)/ (  (len(GLOBAL_SPOTS_STAMPS_CORR) -1 ) *delta_t)
         

            NORM_spots_big_vec.append(NORM_spots)
            NORM_flats_big_vec.append(NORM_flats)
            ######## 

            NORM=NORM_spots

            s_vec_jay=[]
            val0=rates_vec_jay[0]
            for val in rates_vec_jay:
                val=(val-val0)
                print "VAL: ", val
                s_vec_jay.append(val)    

            s_vec_jay=np.array(s_vec_jay)        
           
            ### For central pixel, record F_i/<F_i>  - 1 and F_i - <F_i>
            if i == xc and j == yc:
                r=(rates_vec_jay/np.mean(rates_vec_jay)) - 1
                ratio_fluxes_vec.append(r)
                
                diff=rates_vec_jay - np.mean(rates_vec_jay)
                diff_fluxes_vec.append(diff)
            

                def_spots=100*s_vec_jay/val0
                signal_deficit_central_pixel_spots.append(def_spots)
                



            #NORM=1
            #### Actual f_N metric
            s_vec_jay/=NORM    



            print "CHAZ's metric: delta_sig, delta_sum, s_vec ", delta_sig, delta_sum, s_vec
            print "JAY's metric: NORM, rates_vec_jay, s_vec_jay ", NORM, rates_vec_jay, s_vec_jay

            if np.isnan(s_vec).any() == True or np.isnan(s_vec_jay).any() == True:
               print "Breaking bad"
               siga=True
               break
          

          


            if k == 1 and l == 1:
                c=dict_B[(1,1)][0][1]
                print "central c1: ", c
                sum=0
                for key in [(1,0),(0,1),(1,2),(2,1)]:
                    sum+=dict_B[key][0][1]
                sum/=4.0 # mean contrast in neighbors 
                print "c1*1000, <c1_neighbors>*1000: ", c*1000, sum*1000
                fc=1000*(c-sum)
                c1=dict_B[(1,1)][0][1]*1000
                c2=dict_B[(1,1)][0][2]
                b = 2*(c1)*(c2)/fc

                samples=range(1, len(s_vec_jay) + 1)
                err=np.repeat(1, len(s_vec_jay))
		m, m_err=linear_fit_m (samples, s_vec_jay, err)
                print " "
                print " "
                print " "
                print "F_C, B, c1, c2, m, m_err, m/F_c: ", fc,  b, c1, c2, m , m_err, m/fc
                ###### DEBUGGING: seems there are some constants messing to get B
                new_B = (m/fc)*(NORM/(val0*delta_t/1000))   # delta_s in miliseconds            
                print "old B, new B, NORM, extra factor: ", m/fc, new_B,  NORM, (NORM/(val0*delta_t/1000))
                print " "
                print " "
                print " "

                f=open(out_dir+"jay_B.dat", 'a+')
                line="%g %g %g %g %g %g %g %g\n" %(fc,b, c1, c2, m, m_err, m/fc, new_B)
                f.write(line)
                f.close()                

  
            ramps_dict_all[(l,k)].append(s_vec)
            ramps_dict_all_jay[(l,k)].append(s_vec_jay)
    
            ##
            ramps_dict_all_jay_mean_signal[(l,k)].append(delta_mean_signal)
            


    if plot_flag == True:   
        stamp=GLOBAL_SPOTS_STAMPS_CORR[-1] 
        #fig2.suptitle('FLATS.' + " Center: (%g,%g) " %(xc,yc))
        #pp.savefig(fig2)
        

        #fig3.suptitle('SPOTS.'+" Center: (%g,%g) " %(xc,yc))
        #pp.savefig(fig3)

        #from matplotlib.colors import LogNorm 
        fig1=plt.figure()
        ax=fig1.add_subplot (2,2,1)
        cax=ax.imshow(stamp, cmap=cm.gray, origin="lower", interpolation='nearest')#, vmin=mean_science - mean_science/factor, vmax=mean_science + mean_science/factor)
        for (j,i),label in np.ndenumerate(stamp):
            print j,i, label
            label="%g e/s \n (%g e)" %(label/(time_darks[-1]/1e3), label)
            ax.text(i,j,label,ha='center',va='center', color='red', size=4)
        ax.set_title("Center: (%g,%g) " %(xc,yc))
        #PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
        #plt.colorbar(PCM, ax=ax)
        #fig.colorbar(cax)
        fig.colorbar(cax)
       
        ax=fig1.add_subplot (2,2,2)
        cax=ax.imshow(c2_mat_flats, cmap=cm.gray, origin="lower", interpolation='nearest' )
        ax.set_title("C2 flat")
        #PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
        #plt.colorbar(PCM, ax=ax)
        for (j,i), label_val in np.ndenumerate(c2_mat_flats):
            label_err=c2_mat_flats_err[j,i]
            label="%.3g +/- \n %.3g" %(label_val, label_err)
            print j,i, label
            ax.text(i,j,label,ha='center',va='center', size=5, color='red')
        fig.colorbar(cax)


        ax=fig1.add_subplot (2,2,3)
        cax=ax.imshow(c2_mat_spots, cmap=cm.gray, origin="lower", interpolation='nearest')
        ax.set_title("C2 spots")
        #PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
        #plt.colorbar(PCM, ax=ax)
        for (j,i),label_val in np.ndenumerate(c2_mat_spots):
            label_err=c2_mat_spots_err[j,i]
            label="%.3g +/- \n %.3g" %(label_val, label_err)
            print j,i, label
            ax.text(i,j,label,ha='center',va='center', size=5, color='red')
        fig.colorbar(cax)

        ax=fig1.add_subplot (2,2,4)
        cax=ax.imshow(c2_mat_spots_corr, cmap=cm.gray, origin="lower", interpolation='nearest')#, vmin=mean_science - mean_science/factor, vmax=mean_science + mean_science/factor)
        ax.set_title("C2 spots after NL correction")
        #PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
        #plt.colorbar(PCM, ax=ax)
        for (j,i),label_val in np.ndenumerate(c2_mat_spots_corr):
            label_err=c2_mat_spots_corr_err[j,i]
            label="%.3g +/- \n %.3g" %(label_val, label_err)
            print j,i, label
            ax.text(i,j,label,ha='center',va='center', size=5, color='red')
        fig.colorbar(cax)

        plt.tight_layout()

        pp.savefig(fig1)


    SUM_STAMP_2D.append(GLOBAL_SPOTS_STAMPS_CORR[-1])
    SUM_STAMP_2D_FLATS.append(GLOBAL_FLATS_STAMPS_CORR[-1])   
    SUM_STAMP_2D_DARKS.append(GLOBAL_DARKS_STAMPS_CORR[-1])

    super_counter+=1

SUM_STAMP_2D=np.array(SUM_STAMP_2D)
MEDIAN_STAMP=np.median(SUM_STAMP_2D, axis=0)

SUM_STAMP_2D_FLATS=np.array(SUM_STAMP_2D_FLATS)
MEDIAN_STAMP_FLATS=np.median(SUM_STAMP_2D_FLATS, axis=0)


SUM_STAMP_2D_DARKS=np.array(SUM_STAMP_2D_DARKS)
MEDIAN_STAMP_DARKS=np.median(SUM_STAMP_2D_DARKS, axis=0)



######## Plot average stamp with numbers


### FLATS
fig1=plt.figure()
ax=fig1.add_subplot (1,1,1)
cax=ax.imshow(MEDIAN_STAMP, cmap=cm.gray, origin="lower", interpolation='nearest')#, vmin=mean_science - mean_science/factor, vmax=mean_science + mean_science/factor)
for (j,i),label in np.ndenumerate(MEDIAN_STAMP):
            print j,i, label
            label="%g e$^{-}$/s \n (%g e$^{-}$)" %(label/(time_darks[-1]/1e3), label)
            if j == 1 and i == 1:
                ax.text(i,j,label,ha='center',va='center', color='black', size=11)
            else:
                ax.text(i,j,label,ha='center',va='center', color='white', size=11) 

pp.savefig()

#DARKS
fig1=plt.figure()
ax=fig1.add_subplot (1,1,1)
cax=ax.imshow(MEDIAN_STAMP_DARKS, cmap=cm.gray, origin="lower", interpolation='nearest')#, vmin=mean_science - mean_science/factor, vmax=mean_science + mean_science/factor)
for (j,i),label in np.ndenumerate(MEDIAN_STAMP_DARKS):
            print j,i, label
            label="%g e$^{-}$/s \n (%g e$^{-}$)" %(label/(time_darks[-1]/1e3), label)
            if j == 1 and i == 1:
                ax.text(i,j,label,ha='center',va='center', color='black', size=11)
            else:
                ax.text(i,j,label,ha='center',va='center', color='white', size=11)

pp.savefig()


#SPOTS
fig1=plt.figure()
ax=fig1.add_subplot (1,2,1)
cax=ax.imshow(MEDIAN_STAMP, cmap=cm.gray, origin="lower", interpolation='nearest')#, vmin=mean_science - mean_science/factor, vmax=mean_science + mean_science/factor)
for (j,i),label in np.ndenumerate(MEDIAN_STAMP):
            print j,i, label
            label="%g e$^{-}$/s \n (%g e$^{-}$)" %(label/(time_darks[-1]/1e3), label)
            if j == 1 and i == 1:
                ax.text(i,j,label,ha='center',va='center', color='black', size=9)
            else:
                ax.text(i,j,label,ha='center',va='center', color='white', size=9)

#ax.set_title("Center: (%g,%g) " %(xc,yc))
#PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
#plt.colorbar(PCM, ax=ax)
#fig.colorbar(cax)
#fig.colorbar(cax)

ax=fig1.add_subplot (1,2,2)
cax=ax.imshow(MEDIAN_STAMP_FLATS, cmap=cm.gray, origin="lower", interpolation='nearest')#, vmin=mean_science - mean_science/factor, vmax=mean_science + mean_science/factor)
for (j,i),label in np.ndenumerate(MEDIAN_STAMP_FLATS):
            print j,i, label
            label="%g e$^{-}$/s \n (%g e$^{-}$)" %(label/(time_darks[-1]/1e3), label)
            if j == 1 and i == 1:
                ax.text(i,j,label,ha='center',va='center', color='black', size=9)
            else:
                ax.text(i,j,label,ha='center',va='center', color='white', size=9)


pp.savefig()





######





big_vec_c2_flats=np.array(big_vec_c2_flats)
big_vec_c2_spots=np.array(big_vec_c2_spots)
#residual_vec=np.array(residual_vec)

ratio_fluxes_vec=np.array(ratio_fluxes_vec)
diff_fluxes_vec=np.array(diff_fluxes_vec)


#signal_flat_SUPER_VEC=np.array(signal_flat_SUPER_VEC)
#signal_spot_SUPER_VEC=np.array(signal_spot_SUPER_VEC)

#median_flat_signal = np.median(signal_flat_SUPER_VEC, axis=0)
#median_spot_signal = np.median(signal_spot_SUPER_VEC, axis=0)



signal_deficit_central_pixel_spots=np.array(signal_deficit_central_pixel_spots)
#signal_deficit_central_pixel_flats=np.array(signal_deficit_central_pixel_flats)

median_def_spot_signal=np.median(signal_deficit_central_pixel_spots, axis=0)
#median_def_flat_signal=np.median(signal_deficit_central_pixel_flats, axis=0)

NORM_flats_big_vec=np.array(NORM_flats_big_vec)
NORM_spots_big_vec=np.array(NORM_spots_big_vec)


np.savetxt (out_dir+"jay_NORM_spots.dat", NORM_spots_big_vec)
np.savetxt (out_dir+"jay_NORM_flats.dat", NORM_flats_big_vec)
np.savetxt (out_dir+"jay_c2_flats.dat", big_vec_c2_flats)
np.savetxt (out_dir+"jay_c2_spots.dat", big_vec_c2_spots)
#np.savetxt (out_dir+"jay_residual_center_pixel_flat.dat", residual_vec)
np.savetxt (out_dir+"jay_ratio_fluxes_center_pixel.dat", ratio_fluxes_vec)
np.savetxt (out_dir+"jay_diff_fluxes_center_pixel.dat", diff_fluxes_vec)
#np.savetxt (out_dir+"jay_median_flux_flats_center_pixel.dat", median_flat_signal)
#np.savetxt (out_dir+"jay_median_flux_spots_center_pixel.dat", median_spot_signal)

#np.savetxt (out_dir+"jay_median_deficit_flux_flats_center_pixel.dat", median_def_flat_signal)
np.savetxt (out_dir+"jay_median_deficit_flux_spots_center_pixel.dat", median_def_spot_signal)


##### Residual of NL model for each pixel. Print a separate file for each pixel. 

for key in dict_residual:
    vec=np.array(np.array(dict_residual[key]))
    pixel_number=get_pixel_index[key]
    np.savetxt (out_dir+"jay_residual_pixel_%g_flat.dat" %pixel_number, vec)

    print "vec 1: ", vec


    vec=np.array(np.array(dict_median_flux_flats[key]))
    vec=np.median(vec,axis=0)
    np.savetxt (out_dir+"jay_median_flux_flats_pixel_%g.dat" %pixel_number, vec)

    print "vec 2: ", vec


    vec=np.array(np.array(dict_median_flux_spots[key]))
    vec=np.median(vec,axis=0)
    np.savetxt (out_dir+"jay_median_flux_spots_pixel_%g.dat" %pixel_number, vec)

    print "vec 3: ", vec


print "FINAL: ", final
print "len(x_int[0:end]): ", len(x_int[0:end])
size_final=np.array(size_final)
central_signal_size_final=np.array(central_signal_size_final)
print len(size_final), size_final.shape
size_final_mean=np.nanmean(size_final, axis=0)
size_final_err=np.nanstd(size_final, axis=0)/np.sqrt(len(size_final[~np.isnan(size_final)]))

central_signal_size_final_mean=np.nanmean(central_signal_size_final, axis=0)



f=open(out_dir+"jay_relative_size.dat", 'w')
for a,b,c in zip (size_final_mean, size_final_err, central_signal_size_final_mean):
    line="%g %g %g \n"%(a,b,c) 
    f.write(line)
f.close()

print "size_final_mean, size_final_err: ", size_final_mean, size_final_err
print len(size_final_mean)



size_final_mean=(size_final_mean - size_final_mean[0])/size_final_mean[0]

fig=plt.figure()
ax=fig.add_subplot(111)
plt.errorbar (range(1,len(size_final_mean)+1), size_final_mean, yerr=size_final_err, markersize=7)
ax.set_xlabel('Sample number')
ax.set_ylabel('Normalized change in size')
#plt.ylim([-0.001,0.15])
plt.xlim([0,len(size_final_mean)+1])
pp.savefig()

##### Plot histogram of measured C2 
#print "dict_c2_flat: ", dict_c2_flat
#print "dict_c2_flat_corr: ", dict_c2_flat_corr


true_NL=-7.76e-7


"""
fig=plt.figure()
for key in dict_c2_flat:
    ax = fig.add_subplot(3,3,int(key))
    #n, bins, patches_out = ax.hist(dict_c2_flat[key], 50, normed=False, color='blue', histtype='step', alpha=0.75 , label='flat, N: %g' %(len(dict_c2_flat[key])))
    #n, bins, patches_out = ax.hist(dict_c2_spot[key], 50, normed=False, color='red', histtype='step', alpha=0.75 , label='spot, N: %g' %(len(dict_c2_spot[key])))
    #n, bins, patches_out = ax.hist(dict_c2_dark[key], 50, normed=False, color='green', histtype='step', alpha=0.75 , label='dark, N: %g' %(len(dict_c2_dark[key])))


    # Do sigma clipping first 
    temp_mean, temp_sigma, indices = sigma_clip.sigma_clip(dict_c2_flat_corr[key],niter=6, nsig=sigma_cut, get_indices=True, verbose=False)
    vector=np.array(dict_c2_flat[key])[indices]
    print "flat, max, min, mean: ", np.max(vector), np.min(vector), np.mean(vector)
    #true_NL_dark=dict_c2_dark_corr[key]
    #true_NL_flat=dict_c2_flat_corr[key]
    #true_NL_spot=dict_c2_spot_corr[key]

    #plt.plot (np.array(dict_c2_dark[key]) - true_NL_dark, 'g.', alpha=0.2)
    plt.plot (np.array(dict_c2_flat[key]) - true_NL, 'b.', alpha=0.2)
    #plt.plot (np.array(dict_c2_spot[key]) - true_NL_spot, 'r.', alpha=0.2)
    
    #plt.plot (np.array(dict_c2_dark_corr[key]), 'g.', alpha=0.2)
    #plt.plot (np.array(vector) - true_NL, 'b.', alpha=0.2)
    #plt.plot (np.array(dict_c2_spot_corr[key]), 'r.', alpha=0.2)
    #plt.ylim([np.min(vector), np.max(vector)])
    ax.set_title('Fitted C2, flat', size=8)
    ax.legend(loc=loc_label , fancybox=True, ncol=1, numpoints=1, prop = prop)
    ax.tick_params(axis='both', which='major', labelsize=7)
pp.savefig()


fig=plt.figure()
for key in dict_c2_flat:
    ax = fig.add_subplot(3,3,int(key))
    #n, bins, patches_out = ax.hist(dict_c2_flat[key], 50, normed=False, color='blue', histtype='step', alpha=0.75 , label='flat, N: %g' %(len(dict_c2_flat[key])))
    #n, bins, patches_out = ax.hist(dict_c2_spot[key], 50, normed=False, color='red', histtype='step', alpha=0.75 , label='spot, N: %g' %(len(dict_c2_spot[key])))
    #n, bins, patches_out = ax.hist(dict_c2_dark[key], 50, normed=False, color='green', histtype='step', alpha=0.75 , label='dark, N: %g' %(len(dict_c2_dark[key])))

    temp_mean, temp_sigma, indices = sigma_clip.sigma_clip(dict_c2_spot_corr[key],niter=6, nsig=sigma_cut, get_indices=True, verbose=False)
    vector=np.array(dict_c2_spot[key])[indices]
    print "spot, max, min, mean: ", np.max(vector), np.min(vector), np.mean(vector)
    #true_NL_dark=dict_c2_dark_corr[key]
    #true_NL_flat=dict_c2_flat_corr[key]
    #true_NL_spot=dict_c2_spot_corr[key]

    #plt.plot (np.array(dict_c2_dark[key]) - true_NL_dark, 'g.', alpha=0.2)
    #plt.plot (np.array(dict_c2_flat[key]) - true_NL_flat, 'b.', alpha=0.2)
    plt.plot (np.array(dict_c2_spot[key]) - true_NL, 'r.', alpha=0.2)

    #plt.plot (np.array(dict_c2_dark_corr[key]), 'g.', alpha=0.2)
    #plt.plot (np.array(dict_c2_flat_corr[key]), 'b.', alpha=0.2)
    #plt.plot (np.array(vector) - true_NL, 'r.', alpha=0.2)
    #plt.ylim([-0.5e-5, 0.5e-5])
    #plt.ylim([np.min(vector), np.max(vector)])
    #ax.set_title('Fitted C2, spot', size=8)
    ax.legend(loc=loc_label , fancybox=True, ncol=1, numpoints=1, prop = prop)
    ax.tick_params(axis='both', which='major', labelsize=7)
pp.savefig()
"""



#pp.close()
#sys.exit(1)



#cmd="rm science_andres.fits"; run_shell_cmd(cmd)


max_vec, sigma_x_vec = np.array(max_vec), np.array(sigma_x_vec)
print max_vec, sigma_x_vec
 
flux_filtered=np.array(flux_filtered)

fig = plt.figure()
ax = fig.add_subplot(111)
n, bins, patches_out = ax.hist(flux_filtered, 50, normed=False, facecolor='red', histtype='step', alpha=0.75) # , label='Mean: %g \n Scatter:%g'%(mean, scatter))
ax.set_title('Histogram of flux (sum in each stamp, last fream)', size=8)
ax.legend(loc=loc_label , fancybox=True, ncol=1, numpoints=1, prop = prop)
ax.tick_params(axis='both', which='major', labelsize=7)

#ax = fig.add_subplot(212)
#ax.scatter (max_vec, sigma_x_vec, color='r')
#ax.scatter (max_vec, sigma_y_vec, color='b')
pp.savefig()
#pp.close()
#sys.exit()

print "END OF FOR LOOP"
print len(ramps_dict_all[(1,1)]), len(ramps_dict_all[(2,2)]), len(ramps_dict_all[(0,1)]), len(ramps_dict_all_jay[(0,1)])
print "SUPER counter: ", super_counter

###### Get the ramps for the bins in flux 


def get_ramps (x_int, y_int, flux, ramps_dict):

    for xc,yc,f in zip(x_int,y_int,flux):    
        #xc,yc=int(xc)-1, int(yc)-1
        #if xc + 2 > 2047: continue
        #if yc + 2 > 2047: continue
        stamp=GLOBAL_SPOTS[-1][yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end]
        print "Inside get_ramps: FLUX_AUTO vs stamp.sum(): ", f, stamp.sum()
        f=stamp.sum()

        print "stamp shape:: ", stamp.shape

        c0_mat_spots=0*stamp

        c2_mat_flats=0*stamp
        c2_mat_spots=0*stamp
        c2_mat_spots_corr=0*stamp

        c2_mat_flats_err=0*stamp
        c2_mat_spots_err=0*stamp
        c2_mat_spots_corr_err=0*stamp

        c1_mat_spots=0*stamp



        counter=1
        #siga=False
        #TEMP_FLAG=False
        for i,k in zip(stamp_range_global_x, stamp_range_local):
            siga=False
            for j,l in zip(stamp_range_global_y, stamp_range_local):
               
                stamp_spots=GLOBAL_SPOTS[:, yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end]
                stamp_flats=GLOBAL_FLATS[:, yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end]
                stamp_darks=darks[:, yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end]
                
                plot_pixel_ramp(ramp=stamp_flats, time=time_flats, i=k,j=l, fig=fig2, counter=counter, plot_flag=False)
                plot_pixel_ramp(ramp=stamp_spots, time=time_darks, fig=fig3, i=k,j=l, counter=counter, fmt='b-o', plot_flag=False)

                p_flat, cov_flat, chi2_flat, flag1 =fit_pixel_ramp (ramp=stamp_flats,time=time_flats, i=k, j=l, pinit=pinit)
                #if flag1 == True: continue
                p_spot, cov_spot, chi2_spot, flag2 =fit_pixel_ramp (ramp=stamp_spots,time=time_darks, i=k, j=l, pinit=pinit)
                #if flag2 == True: continue
                p_dark, cov_dark, chi2_dark, flag3 =fit_pixel_ramp (ramp=stamp_darks, time=time_darks, i=k, j=l, pinit=pinit)
                #if flag3 == True: continue


                p_spot_corr, cov_spot_corr, chi2_spot_corr, flag4, flux_rate_spots, delta_signal_spots, delta_sum_spots, signal_corrected_spots = correct_and_fit_pixel_ramp (ramp=stamp_spots, dark_ramp=stamp_darks, time=time_darks, i=k, j=l, pinit=pinit, c0_spot=p_spot[0], c0_dark=p_dark[0], c1_spot=p_spot[1], c1_dark=p_dark[1], c2=p_flat[2])
                #if flag4 == True: continue


                if p_spot[0] == 0: 
                   print k,l
                   print p_spot
                   sys.exit(1)

                
                c2_mat_spots[l,k], c2_mat_spots_err[l,k]=p_spot[2], np.diag(cov_spot)[2]
                c2_mat_spots_corr[l,k], c2_mat_spots_corr_err[l,k] =p_spot_corr[2], np.diag(cov_spot_corr)[2]
                c0_mat_spots[l,k]=p_spot[0]
                c1_mat_spots[l,k]=p_spot[1]



                ##### Spots
                #flux_rate_spots-=p_spot_corr[0]
                #flux_rate_spots=flux_rate_spots[1:]/time[1:]
                diff_last_first = GLOBAL_SPOTS[-1][yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end] - GLOBAL_SPOTS[0][yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end]
                delta_t=time_darks[-1] - time_darks[-2]
                NORM = np.sum(diff_last_first)/ ( ( len(GLOBAL_SPOTS) -1) *delta_t)                

                #NORM=1
 
                """ 
                s_vec=[]
                val0=flux_rate_spots[0]
                for val in flux_rate_spots:
                    #val=(val-val0)/ (f/ (time[-1]-time[-2]))
                    val=(val-val0)/NORM
                    print "VAL: ", val
                    s_vec.append(val)
                ramps_dict[(l,k)].append(s_vec)
                counter+=1
                print "pixel counter, i (along x), j (along y): ", counter, i, j 
                """

                #s_vec=[]
                #for p,s in zip(delta_signal_spots, delta_sums_spots):
                #    s_vec.append(p/s)
                #    print "RATIO P/S: ", p/s
                #s_vec=np.array(s_vec)
                #ramps_dict[(l,k)].append(s_vec)
                #counter+=1
                #if np.isnan(s_vec).any() == True:
                #    siga=True
                #    break
                #### TEMP
                #if (i==xc) and (j==yc) and (s_vec[0]) < 0.3:
                #    TEMP_FLAG=True
                #    break



        # Produce a NL-corrected stamp
        GLOBAL_SPOTS_STAMPS_CORR, corr_stamps_sum=[],[]
        for t, stamp, dark_stamp in zip(time_darks, GLOBAL_SPOTS[:,yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end], darks [:,yc-stamp_end:yc+1+stamp_end, xc-stamp_end:xc+1+stamp_end]):
            corr_stamp_spots = (1/(2*c2_mat_flats))*(-1+np.sqrt(1-4*c2_mat_flats*(c0_mat_spots-stamp)) )
            corr_stamp_darks = (1/(2*c2_mat_flats))*(-1+np.sqrt(1-4*c2_mat_flats*(c0_mat_darks-dark_stamp)) )
            #corr_stamp = corr_stamp_spots - corr_stamp_darks
            #corr_stamp = c0_mat_spots + c1_mat_spots*t  - dark_stamp
            corr_stamp = stamp - dark_stamp
            GLOBAL_SPOTS_STAMPS_CORR.append (corr_stamp)
            corr_stamps_sum.append(np.sum(corr_stamp))


        siga=False
        for i,k in zip(stamp_range_global_x, stamp_range_local):

            for j,l in zip(stamp_range_global_y, stamp_range_local):
                sig=[]
                for sample_s in GLOBAL_SPOTS_STAMPS_CORR:
                    sig.append(sample_s[l,k])

                delta_sig, delta_sum =[],[]
                for counter in range(len(sig)-1):
                    delta_sig.append(sig[counter+1]-sig[counter])
                    delta_sum.append(corr_stamps_sum[counter+1]-corr_stamps_sum[counter])

                delta_sig=np.array(delta_sig)
                delta_sum=np.array(delta_sum)

                s_vec=delta_sig/delta_sum

                if np.isnan(s_vec).any() == True:
                    print "Breaking bad"
                    #siga=True
                    #break
                print "delta_sig, delta_sum, s_vec: ", delta_sig, delta_sum, s_vec
                ramps_dict[(l,k)].append(s_vec)


            if siga ==True: break
        if siga==True: continue        
        #if TEMP_FLAG==True: continue


#get_ramps (x_int_fbin1, y_int_fbin1, flux_fbin1, ramps_dict_fbin1)
#get_ramps (x_int_fbin2, y_int_fbin2, flux_fbin2, ramps_dict_fbin2)
#get_ramps (x_int_fbin3, y_int_fbin3, flux_fbin3, ramps_dict_fbin3)




## Helper function to stack ramps in each pixel ad plot them. To be used in the nested loop immediatly below.
def stack_ramps_and_plot (ax, ramps_dict_mean_signal, ramps_dict, fmt_string, counter_pixel, label=" ", title=" "):
    #print "HOLA "
    (j,i)=get_pair_index[counter_pixel]

    signal_vector=ramps_dict_mean_signal[(j,i)]
    ramps_vector=ramps_dict[(j,i)]
    
    signal_vector=np.array(signal_vector)
    ramps_vector=np.array(ramps_vector)
        
    if not len (ramps_vector) == len (signal_vector): 
        print "ERROR"
        sys.exit()

    #print "HOLA 2"
    #print "counter_pixel: ", counter_pixel
    ##### Stack vectors doing sigma clipping in each component
    #print "len(ramps_vector): ", len(ramps_vector)
    #print "len(ramps_vector[:,1]): ", len(ramps_vector[:,1])
    #print "ramps_vector[:,k]: ", ramps_vector[:,1]
    #print ramps_vector.shape
    #sys.exit()

    stacked_signal, stacked, stacked_err, stacked_numbers=[],[], [], []
    #for k in range(len(ramps_vector[:,1])):
    #for k in [0,1,2,3]:

    #print range(shapes_darks[0] -start_sample_spots -1)
    #sys.exit()

    for k in range(ramps_vector.shape[1]):
        #if len(ramps_vector) == 0: 
        #    stacked_signal.append(0)
        #    stacked.append(0)
        #    stacked_err.append(0)
        #    stacked_numbers.append(1)
        #    continue
        print "Bin: %g, number of objects: %g" %(k, len(ramps_vector[:,k]))
     	#if len(ramps_vector) == 0: continue
        print "ramps_vector[:,k]:", ramps_vector[:,k]
        print "signal_vector[:,k]:", signal_vector[:,k]

        if k == -1:
            print "hola" 
            #stacked_signal.append(0)
            #stacked.append(0)
            #stacked_err.append(0)
            #stacked_numbers.append(1)
        else:

            temp_mean, temp_sigma, indices = sigma_clip.sigma_clip(ramps_vector[:,k],niter=6, nsig=sigma_cut, get_indices=True, verbose=True)
            temp_mean_signal, temp_sigma_signal, indices_signal = sigma_clip.sigma_clip(signal_vector[:,k],niter=6, nsig=sigma_cut, get_indices=True, verbose=True)


    	    if (indices is None)  and (indices_signal is not None):
                n=1
                stacked_signal.append(temp_mean_signal)
                stacked.append(0.0)
                stacked_err.append(0.0)
    	    else:
                stacked_signal.append(temp_mean_signal)
                stacked.append(temp_mean)
                stacked_err.append(temp_sigma)
       	        n=len(indices)
    	    
            print "Number of points in bin %g: %g " %(k,n)
    	    stacked_numbers.append(n)

    stacked_signal = np.array(stacked_signal)
    stacked=np.array(stacked)
    stacked_err=np.array(stacked_err)
    stacked_numbers=np.array(stacked_numbers)

    #stacked=np.mean(ramps_vector, axis=0)   #Old way: just mean, no sigma clipping
    #stacked_err=np.std(ramps_vector, axis=0)

    samples=(range(1, len(stacked)+1))
    #print "Samples, STACKED: ", samples, stacked
    print stacked_signal, stacked_err, stacked_numbers

    plt.errorbar (stacked_signal ,stacked, yerr=stacked_err/np.sqrt(stacked_numbers), fmt=fmt_string, markersize=4, label=label)
    #plt.errorbar (samples,stacked, yerr=stacked_err/np.sqrt(len(ramps_vector[:,1])), fmt=fmt_string, markersize=4, label=label)
    #return (stacked, stacked_err/np.sqrt(stacked_numbers))

   
    if STAMP_SIZE == 5 and counter_pixel==13:  #5
        plt.legend(loc='upper right', fancybox=True, ncol=1, numpoints=1, prop = prop)
        ax.set_yticklabels(ax.get_yticks(), size=5, visible=True)
        #plt.ylim([-3e-2,3e-2])
    elif STAMP_SIZE == 3 and counter_pixel==5:
        plt.legend(loc='upper right', fancybox=True, ncol=1, numpoints=1, prop = prop)
        ax.set_yticklabels(ax.get_yticks(), size=5, visible=True)
        #plt.ylim([-0.03, 0.005])
    else:    
        print "hola nada!!"
        #plt.ylim([-2e-3,5e-3])
    print "range(shapes_darks[0]): ", range(shapes_darks[0])
    #plt.xlim([0, np.max(range(shapes_darks[0]))])  # number of samples
    #plt.ylim([0, -1e-2])

    #ax.fill_between(range(-100,nsamples+100), -2e-3, 2e-3, facecolor='gray', alpha=0.3)
    if STAMP_SIZE == 5 and counter_pixel in [1,6,11,16,21]: #in [1,6,11,16,21]: #in [1,4,7]:
        ax.set_ylabel(r"$f_{N}$", size =9)
    
    if STAMP_SIZE == 5 and counter_pixel in [21,22,23,24,25]: #[21,22,23,24,25]: #[7,8,9]:
        ax.set_xlabel("Frame number (time)", size=7)
    
    if STAMP_SIZE == 5 and counter_pixel in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]:
        ax.set_xticklabels(ax.get_xticks(), size=9, visible=False)
    #if counter_pixel in [2,3,4,5,7,8,9,10,12,14,15,17,18,19,20,22,23,24,25]:
    #    ax.set_yticklabels(ax.get_yticks(), size=9, visible=True)
        

    if STAMP_SIZE == 3 and counter_pixel in [1,4,7]: #in [1,6,11,16,21]: #in [1,4,7]:
        ax.set_ylabel(r"$f_{N}$", size =16)

    if STAMP_SIZE == 3 and counter_pixel in [7,8,9]: #[21,22,23,24,25]: #[7,8,9]:
        ax.set_xlabel("Frame number (time)", size=11)

    if STAMP_SIZE == 3 and counter_pixel in [1,2,3,4,5,6]:
        ax.set_xticklabels(ax.get_xticks(), size=9, visible=False)

    ax.tick_params(labelsize=7)
    #ax.set_title (title + " \n pixel#: %g"%(counter_pixel), size=6)
    print "HELLO "

    return (stacked_signal, stacked, stacked_err/np.sqrt(stacked_numbers))





def stack_ramps_and_plot2 (ax, ramps_dict, fmt_string, counter_pixel, label=" ", title=" "):
    #print "HOLA "
    (j,i)=get_pair_index[counter_pixel]
    ramps_vector=ramps_dict[(j,i)]
    ramps_vector=np.array(ramps_vector)
    #print "HOLA 2"
    #print "counter_pixel: ", counter_pixel
    ##### Stack vectors doing sigma clipping in each component
    #print "len(ramps_vector): ", len(ramps_vector)
    #print "ramps_vector[:,k]: ", ramps_vector[:,1]
    print ramps_vector.shape


    for r in ramps_vector:
        if len(ramps_vector) == 0:
            r=np.zeros(len(GLOBAL_SPOTS))
        plt.errorbar (range(len(r)),r, yerr=None, fmt=fmt_string, markersize=2)


    #plt.errorbar (samples,stacked, yerr=stacked_err/np.sqrt(len(ramps_vector[:,1])), fmt=fmt_string, markersize=4, label=label)
    if STAMP_SIZE == 5 and counter_pixel==13:  #5
        plt.legend(loc='upper right', fancybox=True, ncol=1, numpoints=1, prop = prop)
    elif STAMP_SIZE == 3 and counter_pixel==5:  #5
        plt.legend(loc='upper right', fancybox=True, ncol=1, numpoints=1, prop = prop)
    else:
        print "Incorrect STAMP_SIZE"
        sys.exit(1)

    #plt.xlim([0,nsamples+1])
    plt.xlim([0,range(shapes_darks[0])])
    #plt.ylim([0, -1e-2])

    #ax.fill_between(range(-100,nsamples+100), -2e-3, 2e-3, facecolor='gray', alpha=0.3)
    if STAMP_SIZE == 5 and counter_pixel in [1,6,11,16,21]: #in [1,6,11,16,21]: #in [1,4,7]:
        ax.set_ylabel(r"$f_{N}$", size =8)
    elif STAMP_SIZE == 3 and counter_pixel in [1,4,7]: #in [1,6,11,16,21]: #in [1,4,7]:
        ax.set_ylabel(r"$f_{N}$", size =8)    
    else:
        print "Incorrect STAMP_SIZE"
        sys.exit(1)  


    if STAMP_SIZE == 5 and counter_pixel in [21,22,23,24,25]: #[21,22,23,24,25]: #[7,8,9]:
        ax.set_xlabel("Frame number (time)", size =7)
    elif STAMP_SIZE == 3 and counter_pixel in [7,8,9]: #[21,22,23,24,25]: #[7,8,9]:
        ax.set_xlabel("Frame number (time)", size =11)
    else: 
        print "Incorrect STAMP_SIZE"
        sys.exit(1)

    #if counter_pixel == 5:
        #plt.ylim([-3e-2,3e-2])
    #else:
        #plt.ylim([-6e-3,6e-3])
        #plt.ylim([-7e-6,7e-6])
    ax.tick_params(labelsize=7)
    ax.set_title (title + "  \n pixel#: %g"%(counter_pixel), size=7)




#print "ramps_dict_all: ", ramps_dict_all

"""
fig=plt.figure()
counter_pixel=1
st_vec=np.zeros(25)

#get_pixel_index={(0,0):7, (0,1):4, (0,2):1, (1,0):8, (1,1):5, (1,2): 2, (2,0):9, (2,1):6, (2,2):3 }
for j in range(-2,3):
    for i in range(-2,3):
        for (l,k), val in np.ndenumerate(st):
            print l, k, st[k,l]
            pixel_index=get_pixel_index[(l,k)]
	    st_vec[pixel_index-1]=st[k,l]
        
        ax1=fig.add_subplot(5,5,counter_pixel)
        stack_ramps_and_plot (ax1, ramps_dict_all, 'k-o', counter_pixel, label="All", title='') #'%g'%st_vec[counter_pixel-1])
        #stack_ramps_and_plot (ax1, ramps_dict_fbin1, 'r-o', counter_pixel, label='')# label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut1, cut2))
        #stack_ramps_and_plot (ax1, ramps_dict_fbin2, 'b-o', counter_pixel, label='')#label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut2,cut3))
        #stack_ramps_and_plot (ax1, ramps_dict_fbin3, 'g-o', counter_pixel, label='')#label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut3, cut4), title='%g'%st_vec[counter_pixel-1])

        counter_pixel+=1
plt.tight_layout()
plt.suptitle("Chaz's metric")
pp.savefig(fig)
"""

###PLOT JAY's METRIC
fig=plt.figure()
counter_pixel=1
if stamp_string == 'five': 
    st_vec=np.zeros(25)
    dict={'1':[],'2':[],'3':[],'4':[],'5':[],'6':[],'7':[],'8':[],'9':[],'10':[],\
    '11':[],'12':[],'13':[],'14':[],'15':[],'16':[],'17':[],'18':[],'19':[],'20':[],'21':[],'22':[],\
    '23':[],'24':[],'25':[]}
elif stamp_string == 'three':
    dict={'1':[],'2':[],'3':[],'4':[],'5':[],'6':[],'7':[],'8':[],'9':[]}
    st_vec=np.zeros(9)
else:
    print "error!!! " 
    sys.exit()


for j in range(-stamp_end,1+stamp_end):
    for i in range(-stamp_end,1+stamp_end):
        for (l,k), val in np.ndenumerate(st):
            print l, k, st[k,l]
            pixel_index=get_pixel_index[(l,k)]
            st_vec[pixel_index-1]=st[k,l]

        if stamp_string == 'five': 
            ax1=fig.add_subplot(5,5,counter_pixel)
        elif stamp_string == 'three':
            ax1=fig.add_subplot(3,3,counter_pixel)
        else:
            print "error! "
            sys.exit()
        print i, j
        print " "
        print " "
        print " "
        print " "
        print " "
        x, y,y_err=stack_ramps_and_plot (ax1, ramps_dict_all_jay_mean_signal, ramps_dict_all_jay, 'k-o', counter_pixel, label="All", title='Second metric (Jay)') #'%g'%st_vec[counter_pixel-1])
        #stack_ramps_and_plot (ax1, ramps_dict_fbin1, 'r-o', counter_pixel, label='')# label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut1, cut2))
        #stack_ramps_and_plot (ax1, ramps_dict_fbin2, 'b-o', counter_pixel, label='')#label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut2,cut3))
        #stack_ramps_and_plot (ax1, ramps_dict_fbin3, 'g-o', counter_pixel, label='')#label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut3, cut4), title='%g'%st_vec[counter_pixel-1])
        dict["%g"%counter_pixel]=(y,y_err, x)
        counter_pixel+=1
#plt.suptitle("Jay's metric")
#plt.tight_layout()
pp.savefig(fig)
#pp.close()
#sys.exit()


print dict

f=open(out_dir+'jay_metric.dat', 'w')
for key in dict:
    line="%s "%key
    a=np.concatenate((dict[key][0], dict[key][1], dict[key][2]))
    for b in a:
        line+="%g "%b
    line+="\n"
    f.write(line)
f.close()


#print len(ramps_dict_fbin1[(1,3)]), len(ramps_dict_fbin2[(1,3)])
#sys.exit()


"""
fig=plt.figure()
counter_pixel=1
st_vec=[0,0,0,0,0,0,0,0,0]
#get_pixel_index={(0,0):7, (0,1):4, (0,2):1, (1,0):8, (1,1):5, (1,2): 2, (2,0):9, (2,1):6, (2,2):3 }
for j in range(-start_y,end_y+1):
    for i in range(-start_x,end_x+1):
        for (l,k), val in np.ndenumerate(st):
            print l, k, st[k,l]
            pixel_index=get_pixel_index[(l,k)]
            st_vec[pixel_index-1]=st[k,l]

        ax1=fig.add_subplot(3,3,counter_pixel)
        stack_ramps_and_plot2 (ax1, ramps_dict_all, 'k-o', counter_pixel, label="All", title='%g'%st_vec[counter_pixel-1])
        #stack_ramps_and_plot2 (ax1, ramps_dict_fbin1, 'r-o', counter_pixel, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut1, cut2))
        #stack_ramps_and_plot2 (ax1, ramps_dict_fbin2, 'b-o', counter_pixel, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut2,cut3))
        #stack_ramps_and_plot2 (ax1, ramps_dict_fbin3, 'g-o', counter_pixel, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut3, cut4), title='%g'%st_vec[counter_pixel-1])

        counter_pixel+=1
plt.tight_layout()
pp.savefig(fig)
"""

"""
fig=plt.figure()
counter_pixel=1
st_vec=[0,0,0,0,0,0,0,0,0]
#get_pixel_index={(0,0):7, (0,1):4, (0,2):1, (1,0):8, (1,1):5, (1,2): 2, (2,0):9, (2,1):6, (2,2):3 }
for j in range(-start_y,end_y+1):
    for i in range(-start_x,end_x+1):
        for (l,k), val in np.ndenumerate(st):
            print l, k, st[k,l]
            pixel_index=get_pixel_index[(l,k)]
            st_vec[pixel_index-1]=st[k,l]

        ax1=fig.add_subplot(3,3,counter_pixel)
        stack_ramps_and_plot2 (ax1, ramps_dict_all, 'k-o', counter_pixel, label="All", title='%g'%st_vec[counter_pixel-1])
        #stack_ramps_and_plot2 (ax1, ramps_dict_fbin1, 'r-o', counter_pixel, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut1, cut2))
        #stack_ramps_and_plot2 (ax1, ramps_dict_fbin2, 'b-o', counter_pixel, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut2,cut3))
        #stack_ramps_and_plot2 (ax1, ramps_dict_fbin3, 'g-o', counter_pixel, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut3, cut4), title='%g'%st_vec[counter_pixel-1])

        counter_pixel+=1
plt.tight_layout()
pp.savefig(fig)
"""
"""
fig=plt.figure()
counter_pixel=1
st_vec=[0,0,0,0,0,0,0,0,0]
#get_pixel_index={(0,0):7, (0,1):4, (0,2):1, (1,0):8, (1,1):5, (1,2): 2, (2,0):9, (2,1):6, (2,2):3 }
for j in range(-start_y,end_y+1):
    for i in range(-start_x,end_x+1):
        for (l,k), val in np.ndenumerate(st):
            print l, k, st[k,l]
            pixel_index=get_pixel_index[(l,k)]
            st_vec[pixel_index-1]=st[k,l]

        ax1=fig.add_subplot(3,3,counter_pixel)
        stack_ramps_and_plot2 (ax1, ramps_dict_all_jay, 'k-o', counter_pixel, label="All", title='%g'%st_vec[counter_pixel-1])
        #stack_ramps_and_plot2 (ax1, ramps_dict_fbin1, 'r-o', counter_pixel, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut1, cut2))
        #stack_ramps_and_plot2 (ax1, ramps_dict_fbin2, 'b-o', counter_pixel, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut2,cut3))
        #stack_ramps_and_plot2 (ax1, ramps_dict_fbin3, 'g-o', counter_pixel, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut3, cut4), title='%g'%st_vec[counter_pixel-1])

        counter_pixel+=1

plt.suptitle("Jay's metric")
plt.tight_layout()
pp.savefig(fig)
"""

"""
fig=plt.figure()
counter_pixel=1
st_vec=np.zeros(25)
#get_pixel_index={(0,0):7, (0,1):4, (0,2):1, (1,0):8, (1,1):5, (1,2): 2, (2,0):9, (2,1):6, (2,2):3 }
for j in range(-2,3):
    for i in range(-2,3):
        for (l,k), val in np.ndenumerate(st):
            print l, k, st[k,l]
            pixel_index=get_pixel_index[(l,k)]
            st_vec[pixel_index-1]=st[k,l]

        ax1=fig.add_subplot(5,5,counter_pixel)
        stack_ramps_and_plot2 (ax1, ramps_dict_all, 'k-o', counter_pixel, label="All", title='')#title='%g'%st_vec[counter_pixel-1])
        #stack_ramps_and_plot2 (ax1, ramps_dict_fbin1, 'r-o', counter_pixel, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut1, cut2))
        #stack_ramps_and_plot2 (ax1, ramps_dict_fbin2, 'b-o', counter_pixel, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut2,cut3))
        #stack_ramps_and_plot2 (ax1, ramps_dict_fbin3, 'g-o', counter_pixel, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut3, cut4), title='%g'%st_vec[counter_pixel-1])

        counter_pixel+=1
plt.tight_layout()
pp.savefig(fig)

"""



###### Sum the 8 neighboring pixels 

def plot_surrounding_pixels (fig, ramps_dict, fmt, label):
    surrounding_pixels=[]
    surrounding_pixels_err=[]
    if stamp_string == 'five':
        final = 26
    elif stamp_string == 'three':
        final = 10
    else:
        print "ERROR!!!!!"
        sys.exit()

    for i in range(1,final):
        if STAMP_SIZE == 5 and i==13: 
            continue   #5
        if STAMP_SIZE == 3 and i == 5:
            continue
        (j,i)=get_pair_index[i]
        ramps_vector=ramps_dict[(j,i)]
        ramps_vector=np.array(ramps_vector)
        #stacked=np.mean(ramps_vector, axis=0)
        
        stacked, stacked_err, stacked_numbers=[],[],[]
        #for k in range(len(ramps_vector[0])):
        for k in range(ramps_vector.shape[1]):
            if len(ramps_vector) == 0:
                stacked.append(0)
                stacked_err.append(0)
                stacked_numbers.append(1)
                continue

            if k == -1:
                stacked.append(0)
                stacked_err.append(0)
                stacked_numbers.append(1) 
            else:  
                temp_mean, temp_sigma, indices = sigma_clip.sigma_clip(ramps_vector[:,k],niter=6, nsig=sigma_cut, get_indices=True, verbose=True)
                stacked.append(temp_mean)
                stacked_err.append(temp_sigma)
                if indices is None:
                    n=0
                else:
                    n=len(indices)
                stacked_numbers.append(n)
        
        stacked=np.array(stacked)
        stacked_err=np.array(stacked_err)
        stacked_numbers=np.array(stacked_numbers)
        
        surrounding_pixels.append(stacked)
        surrounding_pixels_err.append (stacked_err/np.sqrt(stacked_numbers))

    surrounding_pixels=np.array(surrounding_pixels)
    surrounding_pixels_err=np.array(surrounding_pixels_err)

    stacked_surrounding_pixels=np.sum(surrounding_pixels, axis=0)  # sum of the means in the 8 pixels
    stacked_surrounding_pixels_err=np.sqrt ( np.sum(surrounding_pixels_err**2, axis=0)  ) ## add in cuadrature


    ##### Stack vectors doing sigma clipping in each component
    if STAMP_SIZE == 5:
        ramps_vector=ramps_dict[(2,2)]  # Number five or (1,1) is central pixel / 13 or 2,2
    elif STAMP_SIZE == 3:
        ramps_vector=ramps_dict[(1,1)]
    else:
        print "error!!!!!"
        sys.exit()
    ramps_vector=np.array(ramps_vector)
    stacked_central_pixel, stacked_central_pixel_err, stacked_central_pixel_numbers =[],[], []
    #for k in range(len(ramps_vector[0])):
    for k in range(ramps_vector.shape[1]):
        if len(ramps_vector) == 0:
            stacked_central_pixel.append(0)
            stacked_central_pixel_err.append(0)
            stacked_central_pixel_numbers.append(1)
            continue

        if k == -1:
            stacked_central_pixel.append(0)
            stacked_central_pixel_err.append(0)
            stacked_central_pixel_numbers.append(1)
        else:
            temp_mean, temp_sigma, indices = sigma_clip.sigma_clip(ramps_vector[:,k],niter=6, nsig=sigma_cut, get_indices=True, verbose=True)
            if indices is None:
                n=1
                stacked_central_pixel.append(0.0)
                stacked_central_pixel_err.append(0.0)
            else:
                stacked_central_pixel.append(temp_mean)
                stacked_central_pixel_err.append(temp_sigma)
                n=len(indices)
            stacked_central_pixel_numbers.append(n) #number of objects

    stacked_central_pixel=np.array(stacked_central_pixel)
    stacked_central_pixel_err=np.array(stacked_central_pixel_err)
    stacked_central_pixel_numbers=np.array(stacked_central_pixel_numbers)
    #stacked_central_pixel=np.mean(ramps_vector, axis=0)

    print "stacked_central_pixel_numbers", stacked_central_pixel_numbers
    print "stacked_central_pixel_numbers", stacked_central_pixel_numbers


    samples=(range(1, len(stacked_central_pixel)+1))
    print len(samples), len(stacked_surrounding_pixels), len (stacked_central_pixel), range(ramps_vector.shape[1])

    dict={'1':[], '2':[]}   
          
    #fig=plt.figure()
    ax=fig.add_subplot(211)
    #ax.fill_between(range(-100,nsamples+100), -2e-3, 2e-3, facecolor='gray', alpha=0.3)

    plt.errorbar (samples,stacked_surrounding_pixels, yerr=stacked_surrounding_pixels_err, fmt=fmt, markersize=9, label=label)
    ax.set_ylabel(r"$f_{N}$", size =14)
    #ax.set_xlabel("Sample number (time)", size =8)
    ax.set_title ("Sum of %g surrounding pixels" %(final-2), size=11)
    plt.legend(loc='upper left', fancybox=True, ncol=1, numpoints=1, prop = prop)
    plt.xlim([0, np.max(range(shapes_darks[0])) +1  ])
    #plt.ylim([0.24,0.31])

    #plt.ylim([-0.8e-7,0.8e-7])

    ax=fig.add_subplot(212)
    #ax.fill_between(range(-100,nsamples+100), -2e-3, 2e-3, facecolor='gray', alpha=0.3)

    plt.errorbar (samples,stacked_central_pixel, yerr=stacked_central_pixel_err/np.sqrt(stacked_central_pixel_numbers), fmt=fmt, markersize=8, label=label)
    ax.set_ylabel(r"$f_{N}$", size =14)
    ax.set_xlabel("Frame number (time)", size =14)
    ax.set_title ("Signal in Central pixel", size=11)
    #plt.legend(loc='upper left', fancybox=True, ncol=1, numpoints=1, prop = prop)
    plt.xlim([0, np.max(range(shapes_darks[0])) +1  ])
    #plt.ylim([0.69,0.76])
    #plt.ylim([-7e-6,7e-6])
    #plt.ylim([-3e-2,3e-2])
    #plt.tight_layout()
    dict['1']=(stacked_surrounding_pixels, stacked_surrounding_pixels_err)
    dict['2']=(stacked_central_pixel, stacked_central_pixel_err/np.sqrt(stacked_central_pixel_numbers)) 
    return dict 


###### Sum all the 9/25 pixels 

def plot_all_pixels (fig, ramps_dict, fmt, label):
    print "Entering 'plot_all_pixels'"
    surrounding_pixels=[]
    surrounding_pixels_err=[]
    if stamp_string == 'five':
        final = 26
    elif stamp_string == 'three':
        final = 10
    else:
        print "error!!!!!"
        sys.exit()

    for i in range(1,final):
        (j,i)=get_pair_index[i]
        ramps_vector=ramps_dict[(j,i)]
        ramps_vector=np.array(ramps_vector)
        #stacked=np.mean(ramps_vector, axis=0)

        stacked, stacked_err, stacked_numbers=[],[],[]
        #for k in range(len(ramps_vector[0])):
        for k in range(ramps_vector.shape[1]):
            if len(ramps_vector) == 0:
                stacked.append(0)
                stacked_err.append(0)
                stacked_numbers.append(1)
                continue

            if k == -1:
                stacked.append(0)
                stacked_err.append(0)
                stacked_numbers.append(1)
            else:
                temp_mean, temp_sigma, indices = sigma_clip.sigma_clip(ramps_vector[:,k],niter=6, nsig=sigma_cut, get_indices=True, verbose=True)
                if indices is None:
                    n=1
                    stacked.append(0.0)
                    stacked_err.append(0.0)
                else:
                    stacked.append(temp_mean)
                    stacked_err.append(temp_sigma)
                    n=len(indices)
                stacked_numbers.append(n)

        stacked=np.array(stacked)
        stacked_err=np.array(stacked_err)
        stacked_numbers=np.array(stacked_numbers)

        surrounding_pixels.append(stacked)
        surrounding_pixels_err.append (stacked_err/np.sqrt(stacked_numbers))

    surrounding_pixels=np.array(surrounding_pixels)
    surrounding_pixels_err=np.array(surrounding_pixels_err)

    stacked_surrounding_pixels=np.sum(surrounding_pixels, axis=0)  # sum of the means in the 8 pixels
    stacked_surrounding_pixels_err=np.sqrt ( np.sum(surrounding_pixels_err**2, axis=0)  ) ## add in cuadrature


    samples=(range(1, len(stacked_surrounding_pixels)+1))

    #fig=plt.figure()
    ax=fig.add_subplot(211)
    #ax.fill_between(range(-100,nsamples+100), -2e-3, 2e-3, facecolor='gray', alpha=0.3)

    plt.errorbar (samples,stacked_surrounding_pixels, yerr=stacked_surrounding_pixels_err, fmt=fmt, markersize=9, label=label)
    ax.set_ylabel(r"$f_{N}$", size =14)
    #ax.set_xlabel("Sample number (time)", size =8)
    ax.set_title ("Sum of all pixels", size=11)
    plt.legend(loc='upper left', fancybox=True, ncol=1, numpoints=1, prop = prop)
    plt.xlim([0, np.max(range(shapes_darks[0])) +1  ])
    #plt.ylim([0.9,1.1])

    #plt.ylim([-0.8e-7,0.8e-7])

    ax=fig.add_subplot(212)
    #ax.fill_between(range(-100,nsamples+100), -2e-3, 2e-3, facecolor='gray', alpha=0.3)

    plt.errorbar (samples,stacked_surrounding_pixels, yerr=stacked_surrounding_pixels_err, fmt=fmt, markersize=9, label=label)
    ax.set_ylabel(r"$f_{N}$", size =14)
    #ax.set_xlabel("Sample number (time)", size =8)
    ax.set_title ("Sum of all pixels (zoomed in)", size=11)
    plt.legend(loc='upper left', fancybox=True, ncol=1, numpoints=1, prop = prop)
    plt.xlim([0, np.max(range(shapes_darks[0])) +1  ])
    #plt.ylim([0.96,1.01])
    #plt.tight_layout()
    print "Exiting 'plot_all_pixels'"


"""
fig=plt.figure()
plot_surrounding_pixels (fig, ramps_dict_all,'k-o', label="All")
#plot_surrounding_pixels (fig, ramps_dict_fbin1,'r-o', label='')#, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut1, cut2))
#plot_surrounding_pixels (fig, ramps_dict_fbin2,'b-o', label='')#, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut2,cut3))
#plot_surrounding_pixels (fig, ramps_dict_fbin3,'g-o', label='')#, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut3,cut4))
plt.suptitle("Chaz's metric")
pp.savefig()
"""

"""
fig=plt.figure()
plot_all_pixels (fig, ramps_dict_all,'k-o', label="All")
#plot_all_pixels (fig, ramps_dict_fbin1,'r-o', label='')#, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut1, cut2))
#plot_all_pixels (fig, ramps_dict_fbin2,'b-o', label='')#, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut2,cut3))
#plot_all_pixels (fig, ramps_dict_fbin3,'g-o', label='')#, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut3,cut4))
plt.suptitle("Chaz's metric")
pp.savefig()
"""


print " BEFORE THE END OF THE CODE!!!!!! "

fig=plt.figure()
dict=plot_surrounding_pixels (fig, ramps_dict_all_jay,'k-o', label=" ")   # returns dictionary with data to plot:  1=(x,y,yerr)
#plot_surrounding_pixels (fig, ramps_dict_fbin1,'r-o', label='')#, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut1, cut2))
#plot_surrounding_pixels (fig, ramps_dict_fbin2,'b-o', label='')#, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut2,cut3))
#plot_surrounding_pixels (fig, ramps_dict_fbin3,'g-o', label='')#, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut3,cut4))
#plt.suptitle("Jay's metric")
plt.tight_layout()
pp.savefig()

print dict

f=open(out_dir+'jay_metric_surrounding.dat', 'w')
for key in dict:
    line="%s "%key
    a=np.concatenate((dict[key][0], dict[key][1]))
    for b in a:
        line+="%g "%b
    line+="\n"
    f.write(line)
f.close()


print "IN BETWEEN PLOTTING"

fig=plt.figure()
plot_all_pixels (fig, ramps_dict_all_jay,'k-o', label="All")
#plot_all_pixels (fig, ramps_dict_fbin1,'r-o', label='')#, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut1, cut2))
#plot_all_pixels (fig, ramps_dict_fbin2,'b-o', label='')#, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut2,cut3))
#plot_all_pixels (fig, ramps_dict_fbin3,'g-o', label='')#, label="%g e$^-$ $\leq$ $Z$ $<$ %g e$^-$"%(cut3,cut4))
#plt.suptitle("Jay's metric")
pp.savefig()
print "END OF PROGRAM."

pp.close()
