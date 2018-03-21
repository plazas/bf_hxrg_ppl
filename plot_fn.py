#!/usr/bin/python

import numpy as np
import math
import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf  import PdfPages
import matplotlib.font_manager as fm
from scipy import optimize

import sys

## Input: reads five ascii files: central pixel and four neighbors
## Each file has 4 columns: linear flux spots | % residual spots |  linear flux flats | % residual flats


sc=True   #sigma clipping?

### DATA
pp=PdfPages("fn.pdf")
print "Output PDF: fn.pdf "
#### PLOTS
#### Do the plotting here
#plt.minorticks_on()
#plt.tight_layout()

### We do not have matplotlib 1.1, with the 'style' package. Modify the matplotlibrc file parameters instead
import matplotlib as mpl

mpl.rc('lines', linewidth=1, color='black', linestyle='-')
mpl.rc('font', family='serif',weight='normal', size=9.0 )
mpl.rc('text',  color='black', usetex=False)
mpl.rc('axes',  edgecolor='black', linewidth=1, grid=False, titlesize=10, labelsize=10, labelweight='normal',labelcolor='black')
mpl.rc ('legend', numpoints=1, fontsize='8', shadow=False, frameon=False)


mpl.rc('lines', linewidth=1, color='black', linestyle='-')
mpl.rc('font', family='serif',weight='normal', size=9.0 )
mpl.rc('text',  color='black', usetex=False)
mpl.rc('axes',  edgecolor='black', linewidth=1, grid=False, titlesize=10, labelsize=10, labelweight='normal',labelcolor='black')
#mpl.rc('axes.formatter', limits=[-4,4])
#mpl.rc('axes.formatter', limits=[-4,4])
mpl.rcParams['xtick.major.size']=5
mpl.rcParams['xtick.minor.size']=4
mpl.rcParams['xtick.major.pad']=8
mpl.rcParams['xtick.minor.pad']=8
mpl.rcParams['xtick.labelsize']= '7.0'
mpl.rcParams['xtick.minor.width']= 1.0
mpl.rcParams['xtick.major.width']= 1.0
mpl.rcParams['ytick.major.size']=5
mpl.rcParams['ytick.minor.size']=4
mpl.rcParams['ytick.major.pad']=8
mpl.rcParams['ytick.minor.pad']=8
mpl.rcParams['ytick.labelsize']= '7.0'
mpl.rcParams['ytick.minor.width']= 1.0
mpl.rcParams['ytick.major.width']= 1.0
mpl.rc ('legend', numpoints=1, fontsize='8', shadow=False, frameon=False)

## Plot parameters
#plt.subplots_adjust(hspace=0.01, wspace=0.01)
prop = fm.FontProperties(size=6.0)
marker_size=6.0
alpha=0.7
loc_label = "upper right"
visible_x, visible_y = True, True
inset=False


factor=1.0 


root="/Users/amalagon/NL_plots/ASCII_FILES_TO_PLOT/"

#directory_with_files="H_FILTER_CENTER_PPL/"


def read_3x3_into_dict (file="", flux=False, nframes=4):
    f=open(file)
    if flux:
        dict={"1":[[],[],[]], "2":[[],[],[]], "3":[[],[],[]], "4":[[],[],[]], \
        "5":[[],[],[]], "6":[[],[],[]], "7":[[],[],[]], "8":[[],[],[]], "9":[[],[],[]]}      #### y (f_N), yerr, x (mean signal in electrons)
    else:
        dict={"1":[[],[]], "2":[[],[]], "3":[[],[]], "4":[[],[]], \
        "5":[[],[]], "6":[[],[]], "7":[[],[]], "8":[[],[]], "9":[[],[]]}

    for line in f:
        s=line.split(' ')
        #print s[0]
        for i in range(1,nframes+1):
            dict[s[0]][0].append(float(s[i]))
            dict[s[0]][1].append(float(s[i+nframes]))
            if flux: dict[s[0]][2].append(float(s[i+2*nframes]))
    f.close()

    for key in dict:
        dict[key][0]= np.array(dict[key][0])
        dict[key][1]= np.array(dict[key][1])
        if flux: dict[key][2]= np.array(dict[key][2])

    return dict



def read_surrounding_into_dict (file="", nframes=4):
    f=open(file)
    dict={"1":[[],[]], "2":[[],[]]}
    for line in f:
        s=line.split(' ')
        #print s[0]
        for i in range(1,nframes+1):
            dict[s[0]][0].append(float(s[i]))
            dict[s[0]][1].append(float(s[i+nframes]))
    f.close()

    for key in dict:
        dict[key][0]= np.array(dict[key][0])
        dict[key][1]= np.array(dict[key][1])
    
    return dict





def fitfunc (x, m):
    x=np.array(x)
    return m*x
pinit=[0.1]

def linear_fit (x,y,y_err):
    pfinal, covar=optimize.curve_fit(fitfunc,x, y, p0=pinit, sigma=y_err,  maxfev=100000)
    return pfinal[0], np.sqrt(covar[0])

### Read files
import ConfigParser
Config = ConfigParser.ConfigParser()
Config.read("config_plot_fn.ini")

Config.get('params', 'OutDirRoot')


# PARAMS
root="/Users/amalagon/NL_plots/ASCII_FILES_TO_PLOT/"
temp_dir="MAR14_Y_BAND_CUBIC/"  #"FEB26_QUAD_DARK_SUBTRACTED_ZERO_FRAME_YES/"
nframes=4
nframes_corner=6




ppl_center_dict=read_3x3_into_dict (root+"H_FILTER_CENTER_PPL/" + temp_dir + "jay_metric.dat", flux=True, nframes=nframes)
#ppl_center_dict=read_3x3_into_dict (root+"CENTER_PPL/" + "jay_metric.dat", flux=True, nframes=nframes)


#ppl_corner_dict=read_3x3_into_dict (root+"CORNER_PPL/jay_metric_corner.dat" )
sim_center_dict=read_3x3_into_dict (root+"CENTER_SIM_NL_CORRECTED/jay_metric_oct21.dat" )
sim_center_nc_dict=read_3x3_into_dict (root+"CENTER_SIM_NL_NOT_CORRECTED/jay_metric_sim_nl_not_corr.dat" )
sim_nothing=read_3x3_into_dict (root+"SIMS_NOTHING/jay_metric.dat")
flats_120v100=read_3x3_into_dict (root+"FLATS_120_VS_FLATS_100/jay_metric_oct21.dat")
#sim_center_bf=read_3x3_into_dict(root+"CENTER_SIM_BF_90RAMPS_V7/jay_metric.dat")



#print ppl_center_dict
#print ppl_corner_dict

#samples=range(1, len(ppl_center_dict['1'][0]) + 1)

samples= range(1, len(ppl_center_dict['2'][0]) + 1)   # THIS IS ACTUALLY MEAN SIGNAL PER FRAME ( (sample + 1 + sample) / 2) in ELECTRONS
print ppl_center_dict['2'][2], ppl_center_dict['5'][2]
#sys.exit()

simulations_flux={'1':[402.957, 677.249, 955.99, 1240.43], '2': [1825.97 ,3078.59 ,4363.04, 5672.92], '3':[402.915, 677.282, 957.607, 1241.39], \
'4': [1827.36, 3079.58, 4361.75, 5672.8], '5':[28285.2, 46977, 65505.4, 83830] , '6': [1832.58 ,3088.91, 4374.23 ,5687.13], \
'7': [402.496, 676.444, 954.201, 1238.62], '8': [1826.36, 3079.46, 4365.73, 5677.67],  '9': [403.374 ,677.778, 956.593, 1240.58]}


prop = fm.FontProperties(size=4.5)


### SIMULATIONS 

size_fn=7.5

###### PAGE 1

fig=plt.figure()
for key in ppl_center_dict:
    ax=fig.add_subplot(3,3,int(key))

    plt.errorbar (simulations_flux[key],sim_center_nc_dict[key][0]/factor, yerr=sim_center_nc_dict[key][1], fmt='b-o', markersize=5, label='sim+noise+NL')
    plt.errorbar (simulations_flux[key],sim_center_dict[key][0]/factor, yerr=sim_center_dict[key][1], fmt='g--o', markersize=5, label='sim+noise+NL, NL corr.')
    plt.errorbar (simulations_flux[key],sim_nothing[key][0]/factor, yerr=sim_nothing[key][1], fmt='r-o', markersize=5, label='sim+no noise')
    
    

    if key == '5':
        
        
        m, m_err=linear_fit (ppl_center_dict[key][2] ,ppl_center_dict[key][0], None)
        print "m, m_err",  m, m_err
        
        plt.ylim([-0.1, 0.1])
        plt.yticks([-0.1, -0.05, 0.0, 0.05, 0.1], fontsize=size_fn)
        plt.xticks([20000,55000,90000], fontsize=size_fn)    
        
        if inset == True:
            a = plt.axes([0.52, 0.52, 0.08, 0.08])
            plt.errorbar (samples,sim_center_nc_dict[key][0], yerr=sim_center_nc_dict[key][1], fmt='b-o', markersize=3, label='')
            #plt.errorbar (samples,flats_120v100[key][0], yerr=flats_120v100[key][1], fmt='m-o', markersize=3, label='flat1 vs flat2')
            plt.errorbar (samples,sim_center_dict[key][0], yerr=sim_center_dict[key][1], fmt='g--o', markersize=3, label='')
            plt.errorbar (samples,sim_nothing[key][0], yerr=sim_nothing[key][1], fmt='r-o', markersize=3, label='')
            #plt.errorbar (samples,ppl_center_dict[key][0], yerr=ppl_center_dict[key][1], fmt='k-o', markersize=3, label='PPL data')
            
            plt.ylim([-0.1, 0.1])
            #plt.xlim([0,5])
            #plt.xticks([])
            plt.yticks([-0.1, -0.05, 0.0, 0.05, 0.1], fontsize=size_fn)
            plt.xticks([20000,60000,90000], fontsize=size_fn)

        #plt.ylim([-0.1,0.1])
    elif key in ['2','4','6','8']:
        print " "
        plt.ylim([-0.0018,0.0018])
        plt.yticks([-0.002, -0.001, 0.0, 0.001, 0.002], fontsize=size_fn)
        plt.xticks([1500,4000, 6000], fontsize=size_fn)
        
    else:
        plt.ylim([-0.0004,0.0004])
        plt.yticks([-0.0004, -0.0002, 0.0, 0.0002, 0.0004], fontsize=size_fn)
        plt.xticks([350,850,1350], fontsize=size_fn)
        
    if key in ['1','4','7']:
        ax.set_ylabel(r"$f_{N}$", size =14)
    else:
        ax.set_ylabel("", size =14)
    if key in ['7','8','9']:
        ax.set_xlabel("Signal (e$^-$)", size =9.5)
        #plt.xlim([0,5])
    else:
        ax.set_xlabel("", size =14)
        #plt.xlim([0,5])





    if key == '2':
        plt.legend(loc='upper right', fancybox=True, ncol=1, numpoints=1, prop = prop, handlelength=3)
        #plt.xlim([20000, 90000])
        #plt.xticks(range(20000,90000,5), fontsize=10)

    #if key in ['2','4','6','8']:
        #plt.xlim([1000, 6000])


    plt.tight_layout()
pp.savefig()



factor=1.0
##### PLOT ONLY THE PPL DATA

#### PAGE 2

fig=plt.figure()
for key in ppl_center_dict:
    ax=fig.add_subplot(3,3,int(key))
    plt.errorbar (ppl_center_dict[key][2] , ppl_center_dict[key][0], yerr=ppl_center_dict[key][1], fmt='k-o', markersize=5, label='PPL data')

    if key == '5':
        
        m, m_err=linear_fit (samples,ppl_center_dict[key][0], None)
        print m, m_err
        
        #plt.ylim([-0.023, 0.023])
        #plt.xlim([5e3,8e4])

        #plt.yticks([-0.025, -0.0125, 0.0, 0.0125, 0.025], fontsize=size_fn)
        #plt.xticks([2e4,5e4,7.6e4], fontsize=size_fn)    

        #plt.yticks([-0.005, -0.0025, 0.0, 0.0025, 0.005], fontsize=size_fn)
        #plt.xticks([2e4,4e4,6e4], fontsize=size_fn)         

        
    #plt.ylim([-0.1,0.1])
    elif key in ['2','4','6','8']:
        print " "
        #plt.ylim([-0.006,0.006])
        #plt.xlim([500,2300])

        #plt.yticks([-0.005, -0.0025, 0.0, 0.0025, 0.005], fontsize=size_fn)
        #plt.xticks([2000,5000,8000], fontsize=size_fn)

        #plt.yticks([-0.005, -0.0025, 0.0, 0.0025, 0.005], fontsize=size_fn)
        #plt.xticks([10000,17500,25000], fontsize=size_fn)
        
        
    else:
        print " "
        #plt.yticks([-0.0005, -0.00025, 0.0, 0.00025, 0.0005], fontsize=size_fn)
        #plt.xticks([500,1500,2400], fontsize=size_fn)
       
        #plt.yticks([-0.0025, -0.00125, 0.0, 0.00125, 0.0025], fontsize=size_fn)
        #plt.xticks([5000,10000,15000], fontsize=size_fn)
        
        
    if key in ['1','4','7']:
        ax.set_ylabel(r"$f_{N}$", size =14)
    else:
        ax.set_ylabel("", size =14)
    if key in ['7','8','9']:
        ax.set_xlabel("Signal (e$^-$)", size =9.5)
        #plt.xlim([0,5])
    else:
        ax.set_xlabel("", size =14)
        #plt.xlim([0,5])
    
    
    #if key in ['1', '2', '3', '4', '5', '6']:
    #    ax.set_xticklabels([float(x) for x in ax.get_xticks()], visible=True)
    
    #if key in ['2', '3', '6', '8', '9']:
    #    ax.set_yticklabels([float(y) for y in ax.get_yticks()], visible=True)
    
    
    #if key in ['7', '8', '9']:
    #    ax.set_xticklabels([ x for x in ax.get_xticks()], visible=True)
    
    
    if key == '1':
        plt.legend(loc='upper left', fancybox=True, ncol=1, numpoints=1, prop = prop, handlelength=3)

    #if key == '5':
    #    plt.xlim([15000, 80000])


    #if key in ['2','4','6','8']:
    #    plt.xlim([1000, 8000])

    #for axi in ax.flat:
    #ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    #ax.yaxis.set_major_locator(plt.MaxNLocator(5))

    #ax.grid(b=True)
plt.tight_layout()
pp.savefig()

#pp.close()
#sys.exit()






##### Plot relative sizes

########## PAGE 3


#data_ppl_center=np.genfromtxt(root+"CENTER_PPL/jay_relative_size_oct25.dat")
data_ppl_center=np.genfromtxt(root+"H_FILTER_CENTER_PPL/"+temp_dir+"jay_relative_size.dat")

#data_sim_nl= np.genfromtxt(root+"CENTER_SIM_NL_CORRECTED/relative_size_oct21.dat")
#data_sim_nl_nc= np.genfromtxt(root+"CENTER_SIM_NL_NOT_CORRECTED/relative_size_sim_NL_not_corr.dat")
#data_sim_nothing=np.genfromtxt(root+"SIMS_NOTHING/relative_size.dat")
#data_flats_120v100=np.genfromtxt(root+"FLATS_120_VS_FLATS_100/relative_size_oct21.dat")
#data_sim_bf=np.genfromtxt(root+"CENTER_SIM_BF_90RAMPS_V7/jay_relative_size.dat")



data_vec=[data_ppl_center] #,data_sim_nl_nc, data_sim_nl, data_sim_nothing] #, data_sim_bf]


c=['k','b', 'g','r']#,'m','c' ]  #ppl, nl corrected, nl non corrected, noiseless sims
label_vec=['PPL data'] #, 'sim+noise+NL', 'sim+noise+NL, NL corr.', 'sim+no noise'] #, 'flat fields', 'sim+BF']

prop = fm.FontProperties(size=10)


fig=plt.figure()
for i,data in enumerate(data_vec):
    x=data[:,2]
    y=data[:,0]
    yerr=data[:,1]
    print c[i], y
    size_final_mean=100*(y - y[0])/y[0]
    ax=fig.add_subplot(111)
    plt.errorbar (x, size_final_mean, yerr=yerr, markersize=8, fmt=c[i]+'-o', label=label_vec[i])
    plt.legend(loc='upper left', fancybox=True, ncol=1, numpoints=1, prop = prop, handlelength=3)
    #plt.xlim([0,5])
    ax.set_xlabel("Average of mean signal in central pixel (e$^-$)", size =12)
    ax.set_ylabel("Relative change in size (%)", size=12)
    #ax.grid(b=True)

ax.set_xticklabels([float(x) for x in ax.get_xticks()], visible=True, size=9)
ax.set_yticklabels([float(y) for y in ax.get_yticks()], visible=True, size=9)
ax.xaxis.set_major_locator(plt.MaxNLocator(6))
pp.savefig()


###### Plot B histogram
######  PAGE 4

#from matplotlib.ticker import NullFormatter
#b = 2*(c1)*(c2)/fc
# new_B = (m/fc)*(NORM/(val0*delta_t/1000))   # delta_s in miliseconds            
# fc,b, c1, c2, m, m_err, m/fc, new_B

#data_B=np.genfromtxt (root+"CENTER_PPL/jay_B.dat")
data_B=np.genfromtxt (root+"H_FILTER_CENTER_PPL/"+temp_dir+"jay_B.dat")
#data_B=np.genfromtxt(root+"jay_B_FINAL_DEC4.dat")


fc=data_B[:,0]
b1=data_B[:,1]
c1=data_B[:,2]
c2=data_B[:,3]
m=data_B[:,4]
me=data_B[:,5]
b2=data_B[:,7]   #### CHANGED to new B, not 

#c2_flats=np.genfromtxt (root+"H_FILTER_CENTER_PPL/jay_c2_flats.dat")
#c2_spots=np.genfromtxt (root+"H_FILTER_CENTER_PPL/jay_c2_spots.dat")

c2_flats=np.genfromtxt (root+"CENTER_PPL/jay_c2_flats.dat")
c2_spots=np.genfromtxt (root+"CENTER_PPL/jay_c2_spots.dat")




vec=[ (fc, 'Fc', 'green'), (b1,"B, method 1", 'blue'), (c1, "c1", 'yellow'), (c2, "c2", 'magenta'), (m, "m", 'cyan'), (b2,"B (1/e$^-$)", 'red'), (c2_flats,"$C_2$ flats", 'red'), (c2_spots,"$C_2$ spots", 'blue')  ]


#print "b: ", np.min(b), np.max(b), np.median(b)
#print "fc: ", np.min(fc), np.max(fc), np.median(fc)

b=b2

import esutil as eu

if sc:

  sigma_cut=3
  mean_b, sigma_b, indices = eu.stat.sigma_clip(b,niter=10, nsig=sigma_cut, get_indices=True, verbose=True)

  b=b[indices]
  fc=fc[indices]

  print "mean b, std b, std err of mean, min, max after 3-sigma clipping: ", mean_b, sigma_b, sigma_b/np.sqrt(len(b)), np.min(b), np.max(b)

  mean_fc, sigma_fc, indices = eu.stat.sigma_clip(fc,niter=10, nsig=sigma_cut, get_indices=True, verbose=True)

  b=b[indices]
  fc=fc[indices]

  print "mean fc, std fc, std err of mean after 3-sigma clipping: ", mean_fc, sigma_fc, sigma_fc/np.sqrt(len(fc))

  print "new b: ", np.mean (b), np.std(b), np.std(b)/np.sqrt(len(b))

  print "Min and max b: ", np.min(b), np.max(b)



label_string="Mean: %g \n Std. Error of Mean: %g" %(np.mean (b), np.std(b)/np.sqrt(len(b)))



####### Bin

bin=eu.stat.Binner( fc, b )
#b.dohist(nperbin=10)
bin.dohist(nbin=8)
bin.calc_stats()
x=np.array(bin['xmean'])
y=np.array(bin['ymean'])
y_err=np.array(bin['yerr'])


prop = fm.FontProperties(size=9)


fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot(fc,b ,'r.', alpha=0.7)
ax.errorbar( x, y, yerr=y_err, ecolor = 'b', fmt='b.-', markersize=9, label=label_string)
plt.legend(loc='upper right', fancybox=True, ncol=1, numpoints=1, prop = prop, handlelength=3)
#plt.ylim([-6e-6, 6e-6])
#plt.xlim([1000,7000])
ax.set_xlabel('F_C')
ax.set_ylabel('B, method 2')
pp.savefig()
#pp.close()
print "hello"


"""
fig=plt.figure()
for data in data_vec:
    if data[1] == 'center': frame=5
    if data[1] == 'n1': frame=2
    if data[1] == 'n2': frame=4
    if data[1] == 'n3': frame=6
    if data[1] == 'n4': frame=8
    ax=fig.add_subplot(3,3,frame)
    x1,y1,x2,y2 = data[0][:,0], data[0][:,1], data[0][:,2],data[0][:,3]
    print np.nanmedian(y1), np.nanmean(y1)
    print np.nanmedian(y2), np.nanmean(y2)
    print " "
    plt.scatter (x1,y1, c='r', alpha=0.2)
    plt.scatter (x2,y2, c='b', alpha=0.2)
    plt.ylabel("NL residual", size=7.5)
    plt.xlabel("Signal (e)", size=7.5)
    plt.tight_layout()
#plt.ylim([-5,5])

pp.savefig()
pp.close()
"""

"""
for v in vec:
    mean_v, sigma_v, indices = eu.stat.sigma_clip(v[0],niter=10, nsig=sigma_cut, get_indices=True, verbose=True)
    #mean_v, sigma_v= np.mean (v[0]), np.std(v[0])
    x=v[0][indices]
    #x=v[0]
    fig=fig=plt.figure()
    ax=fig.add_subplot(111)
    n, bins, patches_out = ax.hist(x, 50, normed=False, facecolor=v[2], alpha=0.9, label="Number of spots: %.4g \n Mean: %.4g \n Std. Dev.: %.4g \n Std. Err.: %.4g" %(len(x),mean_v, sigma_v, sigma_v/np.sqrt(len(x))))
    ax.set_title('%s'%v[1], size=13)
    ax.legend(loc='upper left' , fancybox=True, ncol=1, numpoints=1, prop = prop)
    ax.tick_params(axis='both', which='major', labelsize=9)
    #if v[2] == 'red':
    #    ax.set_title('Only read and shot noise')
    pp.savefig()
"""

##### Plot B2 histogram from data and simulations together

########### PAGE 5

mpl.rc('axes.formatter', limits=[-4,4])
mpl.rc('axes.formatter', limits=[-4,4])
mpl.rcParams['xtick.labelsize']= '10.0'


data_B_sim=np.genfromtxt(root+"CENTER_SIM_BF_90RAMPS_V7/"+"jay_B.dat")
b2_sim=data_B_sim[:,7]   # new B, 7, not 6

if sc:

    mean_b2_sim, sigma_b2_sim, indices = eu.stat.sigma_clip(b2_sim,niter=10, nsig=sigma_cut, get_indices=True, verbose=True)
    b2_sim=b2_sim[indices]

    mean_b2, sigma_b2, indices = eu.stat.sigma_clip(b2,niter=10, nsig=sigma_cut, get_indices=True, verbose=True)
    b2=b2[indices]
else:
    mean_b2_sim, sigma_b2_sim = np.mean(b2_sim), np.std(b2_sim)
    mean_b2, sigma_b2 = np.mean(b2), np.std(b2)



fig=fig=plt.figure()
ax=fig.add_subplot(111)
n, bins, patches_out = ax.hist(b2, 50, normed=False, facecolor='red', alpha=0.9, label="PPL data: \n Number of spots: %.4g \n Mean: %.4g \n Std. Dev.: %.4g \n Std. Err.: %.4g" %(len(b2),mean_b2, sigma_b2, sigma_b2/np.sqrt(len(b2))))
ax.legend(loc='upper left' , fancybox=True, ncol=1, numpoints=1, prop = prop)

n, bins, patches_out = ax.hist(b2_sim, 50, normed=False, facecolor='green', alpha=0.9, label="Simulations: \n Number of spots: %.4g \n Mean: %.4g \n Std. Dev.: %.4g \n Std. Err.: %.4g" %(len(b2_sim),mean_b2_sim, sigma_b2_sim, sigma_b2_sim/np.sqrt(len(b2_sim))))
ax.legend(loc='upper left' , fancybox=True, ncol=1, numpoints=1, prop = prop)

ax.tick_params(axis='both', which='major', labelsize=9)
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False) 
#ax.set_title("B (1/e$^{-}$)")
ax.set_xlabel("B (1/e$^{-}$)")

#pp.savefig(bbox_inches="tight")
pp.savefig()


mpl.rc('axes.formatter', limits=[-6,6])
mpl.rc('axes.formatter', limits=[-6,6])
mpl.rcParams['xtick.labelsize']= '10.0'



###### PLOT Residuals of model for each central spot in flats, for central pixel 

#median_flux_flats=np.genfromtxt(root+"H_FILTER_CENTER_PPL/"+"jay_median_flux_flats_center_pixel.dat")
#median_flux_spots=np.genfromtxt(root+"CENTER_PPL/"+"jay_median_flux_spots_center_pixel.dat")

#if not len(median_flux_flats) == len(median_flux_spots):
#    print "Median flux vectors do not have same length"
#    sys.exit()


#median_flux_flats/=1e4


#data_residuals_flats=np.genfromtxt (root+"H_FILTER_CENTER_PPL/"+"jay_residual_center_pixel_flat.dat")   #_oct25.dat")
#data_residuals=np.genfromtxt (root+"CENTER_PPL/"+"jay_residual_center_pixel_flat.dat")
#data_residuals_spots=np.genfromtxt (root+"CENTER_PPL/"+"jay_residual_center_pixel_spots.dat") #_oct25.dat")


#frame_number=range(1, len(data_residuals_flats[0,:])+1)
#print frame_number

#mean_y = np.median(data_residuals_flats, axis=0)

#fig=fig=plt.figure()
#ax=fig.add_subplot(211)
#for y in data_residuals_flats:
#    ax.errorbar( median_flux_flats, y, yerr=None, ecolor = 'b', fmt='b.', markersize=6, label='', alpha=0.3)
#ax.errorbar( frame_number, mean_y, yerr=None, ecolor = 'r', fmt='r.-', markersize=9, label='', alpha=1.0)
#ax.set_xlabel('Frame Number (time)')
#ax.set_ylabel('Residual (%): (data - model) / model')
#ax.grid(b=True)
#plt.xlim([0,6])

#ax=fig.add_subplot(111)
#for y in data_residuals_flats:
    #print median_flux_flats, y, len(median_flux_flats), len(y)
#    ax.errorbar( median_flux_flats, y, yerr=None, ecolor = 'g', fmt='g.', markersize=7, label='', alpha=0.05)
#ax.errorbar( median_flux_flats, mean_y, yerr=None, ecolor = 'r', fmt='r.-', markersize=11, label='', alpha=1.0)
#print median_flux_flats, mean_y
#ax.set_xlabel('Average signal (e$^-$)')
#ax.set_ylabel(r"Residual (%)")
#ax.grid(b=True)
#plt.ylim([-0.25, 0.25])
#plt.xlim([0,6])
#ax.set_yticklabels([float(y) for y in ax.get_yticks()], visible=True, size=9)
#pp.savefig()




print "Before Page 6 "

################# Plot NL relative residuals for all the pixels 
####### PAGE 6 

median_flux_flats_dict={'1': [], '2': [], '3': [], \
                            '4':[], '5': [], '6': [], \
                            '7':[], '8': [], '9': []}



mean_y_dict={'1': [], '2': [], '3': [], \
                            '4':[], '5': [], '6': [], \
                            '7':[], '8': [], '9': []}


fig=fig=plt.figure()
for i in range(1,10):
    median_flux_flats=np.genfromtxt(root+"H_FILTER_CENTER_PPL/"+temp_dir+"jay_median_flux_flats_pixel_%g.dat"%i)
    data_residuals_flats=np.genfromtxt (root+"H_FILTER_CENTER_PPL/"+temp_dir+"jay_residual_pixel_%g_flat.dat" %i)

    #median_flux_flats=np.genfromtxt(root+"CENTER_PPL/"+"jay_median_flux_flats_pixel_%g.dat"%i)
    #data_residuals_flats=np.genfromtxt (root+"CENTER_PPL/"+"jay_residual_pixel_%g_flat.dat" %i)
    



    mean_y = np.median(data_residuals_flats, axis=0)

    ## Save vectors to use below
    median_flux_flats_dict["%s"%i] = median_flux_flats
    mean_y_dict["%s"%i] = mean_y

    
    ax=fig.add_subplot(3,3,i)
    for y in data_residuals_flats:
        ax.errorbar( median_flux_flats, y, yerr=None, ecolor = 'g', fmt='g.', markersize=7, label='', alpha=0.1)
    ax.errorbar( median_flux_flats, mean_y, yerr=None, ecolor = 'r', fmt='r.-', markersize=11, label='', alpha=1.0)    
    if i in [7,8,9]: ax.set_xlabel('Average signal \n (e$^-$)', size=9.5)
    if i in [1,4,7]: ax.set_ylabel(r"Residual (%)", size=9)
    plt.ylim([-0.25, 0.25])
    #ax.set_yticklabels([float(y) for y in ax.get_yticks()], visible=True, size=7)
    plt.yticks([-0.25, 0.0, 0.25], fontsize=size_fn)
    plt.xticks([1.0e4, 4.0e4, 7.0e4, 1.0e5], fontsize=size_fn)
plt.tight_layout()
pp.savefig()



print "Before page 6, part 2"

####### PAGE 6, part 2: Absolute residuals (not relative). Note that I use (almost) the same variables as above. 

median_flux_flats_dict_absolute={'1': [], '2': [], '3': [], \
                            '4':[], '5': [], '6': [], \
                            '7':[], '8': [], '9': []}


mean_y_dict_absolute={'1': [], '2': [], '3': [], \
                            '4':[], '5': [], '6': [], \
                            '7':[], '8': [], '9': []}


fig=fig=plt.figure()
for i in range(1,10):
    median_flux_flats_absolute=np.genfromtxt(root+"H_FILTER_CENTER_PPL/"+temp_dir+"jay_median_flux_flats_pixel_%g.dat"%i)
    data_residuals_flats_absolute=np.genfromtxt (root+"H_FILTER_CENTER_PPL/"+temp_dir+"jay_residual_absolute_pixel_%g_flat.dat" %i)

    #median_flux_flats=np.genfromtxt(root+"CENTER_PPL/"+"jay_median_flux_flats_pixel_%g.dat"%i)
    #data_residuals_flats=np.genfromtxt (root+"CENTER_PPL/"+"jay_residual_pixel_%g_flat.dat" %i)
    



    mean_y_absolute = np.median(data_residuals_flats_absolute, axis=0)

    ## Save vectors to use below
    median_flux_flats_dict_absolute["%s"%i] = median_flux_flats_absolute
    mean_y_dict_absolute["%s"%i] = mean_y_absolute

    
    ax=fig.add_subplot(3,3,i)
    for y in data_residuals_flats_absolute:
        ax.errorbar( median_flux_flats_absolute, y, yerr=None, ecolor = 'b', fmt='b.', markersize=7, label='', alpha=0.1)
    ax.errorbar( median_flux_flats_absolute, mean_y_absolute, yerr=None, ecolor = 'r', fmt='r.-', markersize=11, label='', alpha=1.0)    
    if i in [7,8,9]: ax.set_xlabel('Average signal \n (e$^-$)', size=9.5)
    if i in [1,4,7]: ax.set_ylabel(r"Absolute Residual", size=9)
    plt.ylim([-20, 20])
    #ax.set_yticklabels([float(y) for y in ax.get_yticks()], visible=True, size=7)
    #plt.yticks([-0.25, 0.0, 0.25], fontsize=size_fn)
    plt.xticks([1.0e4, 4.0e4, 7.0e4, 1.0e5], fontsize=size_fn)
plt.tight_layout()
pp.savefig()








######## APPLES TO APPLES COMPARISON OF ERROR INDUCED BY NL ON SIGNAL
### PAGE 8 




#delta_t=3.0 #seconds
#simulated_mean_spot_fluxes={'1': 139.028, '2': 521.738, '3':163.384, \
#                            '4':531.722, '5':5047.73, '6':520.802, \
#                            '7': 147.19, '8': 520.826, '9': 151.987}

## H filter 
#delta_t=3.350
#simulated_mean_spot_fluxes={'1': 845.621, '2': 1105.43, '3':847.195, \
#                            '4':1207.35, '5':2772.44, '6':1192.33, \
#                            '7': 844.803, '8': 1109.37, '9': 845.586}


simulated_mean_spot_fluxes_data=np.genfromtxt (root+"H_FILTER_CENTER_PPL/"+temp_dir+"jay_median_spot_fluxes.dat")
simulated_mean_spot_fluxes={'1': simulated_mean_spot_fluxes_data[0], '2': simulated_mean_spot_fluxes_data[1], '3':simulated_mean_spot_fluxes_data[2], \
                            '4':simulated_mean_spot_fluxes_data[3], '5':simulated_mean_spot_fluxes_data[4], '6':simulated_mean_spot_fluxes_data[5], \
                            '7': simulated_mean_spot_fluxes_data[6], '8': simulated_mean_spot_fluxes_data[7], '9': simulated_mean_spot_fluxes_data[8]}
#delta_t=3.350528
#delta_t=3.0
delta_t=0.837632

    
    
## Dictionary of fluxes for average ramps ine ach pixel of spot:
prop = fm.FontProperties(size=7)


#ppl_center_dict_NO_NORM=read_3x3_into_dict (root+"PPL_DATA_NO_NORM/jay_metric.dat", flux=True)
ppl_center_dict_NO_NORM=read_3x3_into_dict (root+"H_FILTER_CENTER_PPL/"+temp_dir+"jay_metric_no_NORM.dat", flux=True, nframes=nframes)




  

error_NO_NORM={'1': [], '2': [], '3': [], \
                            '4':[], '5': [], '6': [], \
                            '7':[], '8': [], '9': []}


res={'1': [], '2': [], '3': [], \
                            '4':[], '5': [], '6': [], \
                            '7':[], '8': [], '9': []}


residual_functions_dict={'1': [], '2': [], '3': [], \
                            '4':[], '5': [], '6': [], \
                            '7':[], '8': [], '9': []}
    

### 1. Interpolate residual funcion

print "      "
from scipy.interpolate import interp1d
print mean_y_dict, "hola"

#sys.exit(1)

#mean_y=np.array([0.08,0.09,0.005,0.05,0.01])/100
#mean_y=np.array([0.08,0.07,0.06,0.05,0.04])


#Generate the residual function objects by interpolating the mean residual curves per pixel (read data from Fig. 6)

for key in residual_functions_dict:
    median_flux_flats=median_flux_flats_dict[key]
    mean_y=mean_y_dict[key]
    
    median_flux_flats[0]=0.0   # Define zero residual at zero signal 
    mean_y[0]=0.0

    #residual=interp1d(median_flux_flats, np.abs(mean_y)/100) #,bounds_error=False, fill_value=np.mean(np.abs(mean_y))/100)  ## Divide by 100 because mean_y is in %
    residual=interp1d(median_flux_flats, np.abs(mean_y)/100)
    #print median_flux_flats, np.mean(np.abs(mean_y))
    residual_functions_dict[key]=residual

    print "key, median_flux_flats, np.abs(mean_y)/100: ", key, median_flux_flats, np.abs(mean_y)/100




    
def get_error_no_norm_ramp (pixel_number_string, end_sample=5):
    simulated_ramp=np.array(range(1,end_sample+1))*simulated_mean_spot_fluxes[pixel_number_string]*delta_t
    print "key, simulated ramp: ", key, simulated_ramp
    error=simulated_ramp*residual_functions_dict[pixel_number_string](simulated_ramp)
    print "error: ", error
    vec=[]
    for i in range(len(error)-1):
        print "pixel_number_string: ", pixel_number_string, i, (error[i+1]-error[i])/ delta_t
        vec.append((error[i+1]-error[i])/ delta_t)

    vec2=[]
    first=vec[0]
    for x in vec:
        vec2.append(x-first)
    
    return np.array(vec2)


def get_residual_func (pixel_number_string, end_sample=5):
    simulated_ramp=np.array(range(1,end_sample+1))*simulated_mean_spot_fluxes[pixel_number_string]*delta_t
    return residual_functions_dict[pixel_number_string](simulated_ramp)


fig=plt.figure()
for key in ppl_center_dict:
    ax=fig.add_subplot(3,3,int(key))
    #plt.errorbar (ppl_center_dict[key][2] , np.abs(1000*ppl_center_dict_NO_NORM[key][0]), yerr=ppl_center_dict[key][1], fmt='k-o', markersize=5, label='PPL data')
    plt.errorbar (ppl_center_dict[key][2] , (1000*ppl_center_dict_NO_NORM[key][0]), yerr=ppl_center_dict[key][1], fmt='k-o', markersize=5, label='PPL data')


    error_NO_NORM[key]=get_error_no_norm_ramp(key, end_sample=nframes+1)
    print "error_NO_NORM[key: ", error_NO_NORM[key]
    #plt.errorbar (ppl_center_dict[key][2] , np.abs(error_NO_NORM[key]), yerr=None, fmt='r-o', markersize=5, label='NL error')
    plt.errorbar (ppl_center_dict[key][2] ,   (error_NO_NORM[key]), yerr=None, fmt='r-o', markersize=5, label='NL error')
    
    #if key == '5':
        
    #    m, m_err=linear_fit (samples,1000*ppl_center_dict_NO_NORM[key][0]/factor, None)
    #    print m, m_err
        
    #else:
    #    print " "
    #    plt.xlim([500,2300])

    if key in ['1','4','7']:
        ax.set_ylabel(r"|$F_{k}-F_{1}|$ (e$^-$/s)", size =9)
    else:
        ax.set_ylabel("", size =14)
    if key in ['7','8','9']:
        ax.set_xlabel("Signal (e$^-$)", size =9.5)
    else:
        ax.set_xlabel("", size =14)
    
    
    
    
    if key == '5':
        plt.legend(loc='upper left', fancybox=True, ncol=1, numpoints=1, prop = prop, handlelength=3)

    if key == '5':
        #plt.xlim([15000, 80000])
        #plt.yticks([0.0, 50, 100, 150], fontsize=size_fn)
        #plt.xticks([1.5e4,4.5e4,7.5e4], fontsize=size_fn)
        print " "
    elif key in ['1','3','7','9']:
        #plt.yticks([0.0, 1.0, 2.0, 3.0,4.0], fontsize=size_fn)
        #plt.xticks([500,1500,2500], fontsize=size_fn)
        print " "
    else:
        print " "
        #plt.xticks([2e3,5e3,8e3], fontsize=size_fn)
        #plt.yticks([0.0, 10.0, 20.0, 30.0], fontsize=size_fn)

        
plt.tight_layout()
pp.savefig()



#### Plot only the residual function

#### PAGE 8

fig=plt.figure()
for key in ppl_center_dict:
    ax=fig.add_subplot(3,3,int(key))
    res[key]=get_residual_func(key, end_sample=nframes+1)
    print "res[key]: ", res[key]
    plt.plot ( (res[key]), 'b-o', markersize=5, label='Residual function')
    
    if key == '5':
        
        m, m_err=linear_fit (samples,1000*ppl_center_dict_NO_NORM[key][0]/factor, None)
        print m, m_err
        
    else:
        print " "
        #plt.xlim([500,2300])

    if key in ['1','4','7']:
        ax.set_ylabel(r"Interpolated NL residual", size =8)
    else:
        ax.set_ylabel("", size =14)
    if key in ['7','8','9']:
        ax.set_xlabel("Frame number", size =9)
    else:
        ax.set_xlabel("", size =14)
    
    
    
    
    if key == '5':
        plt.legend(loc='upper left', fancybox=True, ncol=1, numpoints=1, prop = prop, handlelength=3)

    if key == '5':
        print " "
        #plt.xlim([15000, 80000])


    if key in ['2','4','6','8']:
        print " "
        #plt.xlim([1000, 8000])

    ax.set_yticklabels([float(y) for y in ax.get_yticks()], visible=True, size=7)
        

plt.tight_layout()
pp.savefig()



##### Plot charge conservation

### PAGE 9

#surr=read_surrounding_into_dict(root+"CENTER_PPL/"+"jay_metric_surrounding.dat")
surr=read_surrounding_into_dict(root+"H_FILTER_CENTER_PPL/"+temp_dir+"jay_metric_surrounding.dat", nframes=nframes)

fig=plt.figure()
samples=range(1, nframes+1)

ax=fig.add_subplot (1,1,1)
ax.errorbar( samples, surr['1'][0]/factor, yerr=surr['1'][1], ecolor = 'b', fmt='b.-', markersize=12, label='', alpha=1.0)
ax.errorbar( samples, surr['2'][0]/factor, yerr=surr['1'][1], ecolor = 'r', fmt='r*-', markersize=12, label='', alpha=1.0)
ax.set_xlabel("Frame difference number", size =12)
ax.set_ylabel(r"$f_N$", size =15)
#plt.ylim([-0.023,0.023])
#plt.xticks([1, 2, 3, 4], fontsize=11)
#plt.yticks([-0.02,-0.01,0.0, 0.01, 0.02], fontsize=11)

#ax.set_xticklabels(['', 1, '', 2, '', 3, '', 4, ''], visible=True, size=8)
#ax.set_yticklabels([float(x) for x in ax.get_yticks()], visible=True, size=8)
#plt.ylim([-0.023,0.023])
plt.tight_layout()
pp.savefig()





mpl.rc('axes.formatter', limits=[-5,5])
######### PLOT OF CORNERS

#### PAGE 10 

#data_corner1=read_3x3_into_dict (root+"CORNER_PPL/SECOND_RUN/REGION1_xc_lt_0_yc_lt_0/"+"jay_metric.dat", flux=True)   #5,6,8,9
#data_corner2=read_3x3_into_dict (root+"CORNER_PPL/SECOND_RUN/REGION2_xc_gt_0_yc_lt_0/"+"jay_metric.dat", flux=True)   #2,3,5,6
#data_corner3=read_3x3_into_dict (root+"CORNER_PPL/SECOND_RUN/REGION3_xc_lt_0_yc_gt_0/"+"jay_metric.dat", flux=True)   #4,5,7,8
#data_corner4=read_3x3_into_dict (root+"CORNER_PPL/SECOND_RUN/REGION4_xc_gt_0_yc_gt_0/"+"jay_metric.dat", flux=True)   #1,2,4,5


data_corner1=read_3x3_into_dict (root+"H_FILTER_CENTER_PPL/FEB15_FILMSTRIP_V2_DARK_SUB_YES_ZERO_FRAME_YES_REGION1/"+"jay_metric.dat", flux=True, nframes=nframes_corner)
data_corner2=read_3x3_into_dict (root+"H_FILTER_CENTER_PPL/FEB15_FILMSTRIP_V2_DARK_SUB_YES_ZERO_FRAME_YES_REGION2/"+"jay_metric.dat", flux=True, nframes=nframes_corner)
data_corner3=read_3x3_into_dict (root+"H_FILTER_CENTER_PPL/FEB15_FILMSTRIP_V2_DARK_SUB_YES_ZERO_FRAME_YES_REGION3/"+"jay_metric.dat", flux=True, nframes=nframes_corner)
data_corner4=read_3x3_into_dict (root+"H_FILTER_CENTER_PPL/FEB15_FILMSTRIP_V2_DARK_SUB_YES_ZERO_FRAME_YES_REGION4/"+"jay_metric.dat", flux=True, nframes=nframes_corner)



mean_corner={"1":[[],[],[]], "2":[[],[],[]], "3":[[],[],[]] , "4":[[],[],[]]}

mean_corner["1"][0]=(data_corner1["5"][0]+data_corner1["6"][0]+data_corner1["8"][0]+data_corner1["9"][0] )*0.25
mean_corner["1"][1]=np.sqrt ( data_corner1["5"][1]**2+data_corner1["6"][1]**2+data_corner1["8"][1]**2+data_corner1["9"][1]**2)
mean_corner["1"][2]=(data_corner1["5"][2]+data_corner1["6"][2]+data_corner1["8"][2]+data_corner1["9"][2] )*0.25


mean_corner["2"][0]=(data_corner2["2"][0]+data_corner2["3"][0]+data_corner2["5"][0]+data_corner2["6"][0] )*0.25
mean_corner["2"][1]=np.sqrt ( data_corner2["2"][1]**2+data_corner2["3"][1]**2+data_corner2["5"][1]**2+data_corner2["6"][1]**2)
mean_corner["2"][2]=(data_corner2["2"][2]+data_corner2["3"][2]+data_corner2["5"][2]+data_corner2["6"][2] )*0.25


mean_corner["3"][0]=(data_corner3["4"][0]+data_corner3["5"][0]+data_corner3["7"][0]+data_corner3["8"][0] )*0.25
mean_corner["3"][1]=np.sqrt ( data_corner3["4"][1]**2+data_corner3["5"][1]**2+data_corner3["7"][1]**2+data_corner3["8"][1]**2)
mean_corner["3"][2]=(data_corner3["4"][2]+data_corner3["5"][2]+data_corner3["7"][2]+data_corner3["8"][2] )*0.25


mean_corner["4"][0]=(data_corner4["1"][0]+data_corner4["2"][0]+data_corner4["4"][0]+data_corner4["5"][0] )*0.25
mean_corner["4"][1]=np.sqrt ( data_corner4["1"][1]**2+data_corner4["2"][1]**2+data_corner4["4"][1]**2+data_corner4["5"][1]**2)
mean_corner["4"][2]=(data_corner4["1"][2]+data_corner4["2"][2]+data_corner4["4"][2]+data_corner4["5"][2] )*0.25



fig=plt.figure()
for key in mean_corner:
    ax=fig.add_subplot (2,2,int(key))
    ax.errorbar( mean_corner[key][2], mean_corner[key][0]/factor, yerr=mean_corner[key][1], ecolor = 'k', fmt='k.-', markersize=6, label='', alpha=1.0)
    #plt.ylim([-0.006,0.006])
    #if key in ['1', '2']:
    #    ax.set_xticklabels([int(x) for x in ax.get_xticks()], visible=False)
    #if key in ['2', '4']:
    #    ax.set_yticklabels([int(y) for y in ax.get_yticks()], visible=False)

    if key in ['3', '4']:
        ax.set_xlabel("Signal (e$^-$)", size =10)

    if key in ['1', '3']:
        ax.set_ylabel(r"$f_N$", size =13)

    #plt.xlim([5000,35000])
    #ax.set_xticklabels([int(x) for x in ax.get_xticks()], visible=True, size=7.5)
    #ax.set_yticklabels([float(y) for y in ax.get_yticks()], visible=True, size=7.5)
    #plt.yticks([-0.005, -0.0025, 0.0, 0.0025, 0.005], fontsize=9.5)
    #plt.xticks([5e3, 10e3, 15e3, 20e3, 25e3, 30e3], fontsize=9.5)
    
plt.tight_layout()
plt.suptitle ("")
pp.savefig()



##### Plot "A" NORM for both flats and spots
### PAGE 11 

"""

data_NORM_flats=np.genfromtxt(root+"H_FILTER_CENTER_PPL/"+"jay_NORM_flats.dat")
data_NORM_spots=np.genfromtxt(root+"H_FILTER_CENTER_PPL/"+"jay_NORM_spots.dat")

#m_flats,  indices = eu.stat.sigma_clip(b2_sim,niter=10, nsig=sigma_cut, get_indices=True, verbose=True)
#b2_sim=b2_sim[indices]

#mean_b2, sigma_b2, indices = eu.stat.sigma_clip(b2,niter=10, nsig=sigma_cut, get_indices=True, verbose=True)
#b2=b2[indices]


fig=fig=plt.figure()
ax=fig.add_subplot(211)
n, bins, patches_out = ax.hist(data_NORM_flats, 50, normed=False, facecolor='red', alpha=0.9) #, label="PPL data \n Number of spots: %.4g \n Mean: %.4g \n Std. Dev.: %.4g \n Std. Err.: %.4g" %(len(b2),mean_b2, sigma_b2, sigma_b2/np.sqrt(len(b2))))
#ax.legend(loc='upper left' , fancybox=True, ncol=1, numpoints=1, prop = prop)

ax=fig.add_subplot(212)
n, bins, patches_out = ax.hist(data_NORM_spots, 50, normed=False, facecolor='green', alpha=0.9)
#n, bins, patches_out = ax.hist(b2_sim, 50, normed=False, facecolor='green', alpha=0.9, label="Simulations \n Number of spots: %.4g \n Mean: %.4g \n Std. Dev.: %.4g \n Std. Err.: %.4g" %(len(b2_sim),mean_b2_sim, sigma_b2_sim, sigma_b2_sim/np.sqrt(len(b2_sim))))
#ax.legend(loc='upper left' , fancybox=True, ncol=1, numpoints=1, prop = prop)

ax.tick_params(axis='both', which='major', labelsize=9)

plt.suptitle("F* (norm of metric)")

pp.savefig()
"""




pp.close()







