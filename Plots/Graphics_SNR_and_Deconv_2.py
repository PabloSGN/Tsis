# ============================ IMPORTS ====================================== #
 
import os
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

from scipy.stats import norm
import matplotlib.patheffects as pe
from scipy.optimize import curve_fit
from scipy.stats import skewnorm
import matplotlib as mp

from scipy.special import erf

import functions_old as ft

# ============================= CONFIG ====================================== #
 
plt.style.use('default')
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"

Input_Map_path = "/media/pablo/TOSHIBA EXT/BACKUP TORRE IAA/Etalon_Correction_Final/NUMERICAL/INPUTS_100x100/InputMaps.fits"



SNR_files = {'200' : {'C' : "/media/pablo/TOSHIBA EXT/BACKUP TORRE IAA/Etalon_Correction_Final/COLLIMATED/RESULTS_REDUCED_SPECTRAL/25000",
                      'T' : "/media/pablo/TOSHIBA EXT/BACKUP TORRE IAA/Etalon_Correction_Final/TELECENTRIC/RESULTS_SNR/25000",
                      'I' : "/media/pablo/TOSHIBA EXT/BACKUP TORRE IAA/Etalon_Correction_Final/NUMERICAL/RESULTS_100/"},
             '150' : {'C' : "/media/pablo/TOSHIBA EXT/BACKUP TORRE IAA/Etalon_Correction_Final/COLLIMATED/RESULTS_REDUCED_SPECTRAL/15000",
                      'T' : "/media/pablo/TOSHIBA EXT/BACKUP TORRE IAA/Etalon_Correction_Final/TELECENTRIC/RESULTS_SNR/15000"},
             '100' : {'C' : "/media/pablo/TOSHIBA EXT/BACKUP TORRE IAA/Etalon_Correction_Final/COLLIMATED/RESULTS_REDUCED_SPECTRAL/8000",
                      'T' : "/media/pablo/TOSHIBA EXT/BACKUP TORRE IAA/Etalon_Correction_Final/TELECENTRIC/RESULTS_SNR/8000"}               
    }
               
Et = {
      'R' : 0.892,
      'n' : 2.3268,
      'd' : 281e-6,
      'theta' : 0,
      'cos_th' : np.cos(np.deg2rad(0)),
      'fnum' : 60
      }

Np = 100
Nwvls = np.arange(5, 22)


Deconv_files = "/media/pablo/TOSHIBA EXT/BACKUP TORRE IAA/Etalon_Correction_Final/COLLIMATED/RESULTS_DECONV_v2/25000"

wavelength = 6173.34388 * 1E-10

c = 299792458

snr = ['100', '150', '200']

# da -> wavelength shift (imperfect)
coefs = [ 10410.29797429, -10410.30295311]

# =========================================================================== #

# LOADING INPUT MAPS

Input = fits.getdata(Input_Map_path)


# fig, axs = plt.subplots(2, 2)
# axs[0, 0].imshow(Inputp[0], cmap = 'gray')
# axs[0, 1].imshow(Inputp[1], cmap = 'bwr')
# axs[1, 0].imshow(Inputi[0], cmap = 'gray')
# axs[1, 1].imshow(Inputi[1], cmap = 'bwr')

# # LOADING DATA

SNR = {'200' : {'DA' : {'C' : np.zeros((len(Nwvls), Np, Np)),
                        'T' : np.zeros((len(Nwvls), Np, Np)),
                        'I' : np.zeros((len(Nwvls), Np, Np))},
                'G' :  {'C' : np.zeros((len(Nwvls), Np, Np)),
                        'T' : np.zeros((len(Nwvls), Np, Np)),
                        'I' : np.zeros((len(Nwvls), Np, Np))}},
       '150' :  {'DA' : {'C' : np.zeros((len(Nwvls), Np, Np)),
                         'T' : np.zeros((len(Nwvls), Np, Np))},
                 'G' :  {'C' : np.zeros((len(Nwvls), Np, Np)),
                         'T' : np.zeros((len(Nwvls), Np, Np))}},
       '100' :  {'DA' : {'C' : np.zeros((len(Nwvls), Np, Np)),
                         'T' : np.zeros((len(Nwvls), Np, Np))},
                 'G' :  {'C' : np.zeros((len(Nwvls), Np, Np)),
                         'T' : np.zeros((len(Nwvls), Np, Np))}}}


for s in snr:
    print('Loading ', s, 'SNR')
    
    if s == '200':
        Configs = ['C', 'T', 'I']
    else:
        Configs = ['C', 'T']
        
    for conf in Configs:
        
        files = sorted(os.listdir(SNR_files[s][conf]))
        
        for file_index, file in enumerate(files):

            
            data = fits.getdata(os.path.join(SNR_files[s][conf], file))
            
            if conf == 'I':
                data = np.transpose(data, axes = (1, 0)).reshape(2, Np, Np)

            elif conf == 'C':
                data = np.transpose(data, axes = (2, 0, 1))
            else:
                data = np.transpose(data[-1], axes = (2, 0, 1))

            print('File : ', os.path.join(SNR_files[s][conf], file), 'loaded')
            
            SNR[s]['G'][conf][file_index] = data[0]
            SNR[s]['DA'][conf][file_index] = data[1]
        
        
DECONV = np.zeros((2, len(Nwvls), Np, Np))      
        
        
files = sorted(os.listdir(Deconv_files))
for file_index, file in enumerate(files):
    data = fits.getdata(os.path.join(Deconv_files, file))
    
    data = np.transpose(data[-1], axes = (2, 0, 1))
    
    DECONV[0, file_index] = data[0]
    DECONV[1, file_index] = data[1]


Means = {'200' : {'DA' : {'C' : np.zeros(len(Nwvls)),
                        'T' : np.zeros(len(Nwvls)),
                        'I' : np.zeros(len(Nwvls))},
                'G' :  {'C' : np.zeros(len(Nwvls)),
                        'T' : np.zeros(len(Nwvls)),
                        'I' : np.zeros(len(Nwvls))}},
       '150' :  {'DA' : {'C' : np.zeros(len(Nwvls)),
                         'T' : np.zeros(len(Nwvls))},
                 'G' :  {'C' : np.zeros(len(Nwvls)),
                         'T' : np.zeros(len(Nwvls))}},
       '100' :  {'DA' : {'C' : np.zeros(len(Nwvls)),
                         'T' : np.zeros(len(Nwvls))},
                 'G' :  {'C' : np.zeros(len(Nwvls)),
                         'T' : np.zeros(len(Nwvls))}}}

STD = {'200' : {'DA' : {'C' : np.zeros(len(Nwvls)),
                        'T' : np.zeros(len(Nwvls)),
                        'I' : np.zeros(len(Nwvls))},
                'G' :  {'C' : np.zeros(len(Nwvls)),
                        'T' : np.zeros(len(Nwvls)),
                        'I' : np.zeros(len(Nwvls))}},
       '150' :  {'DA' : {'C' : np.zeros(len(Nwvls)),
                         'T' : np.zeros(len(Nwvls))},
                 'G' :  {'C' : np.zeros(len(Nwvls)),
                         'T' : np.zeros(len(Nwvls))}},
       '100' :  {'DA' : {'C' : np.zeros(len(Nwvls)),
                         'T' : np.zeros(len(Nwvls))},
                 'G' :  {'C' : np.zeros(len(Nwvls)),
                         'T' : np.zeros(len(Nwvls))}}}

means_deconv = np.zeros((2, len(Nwvls)))
std_deconv = np.zeros((2, len(Nwvls)))

print('Computing statistics')


Errorp = {'200' : {'DA' : {'C' : np.zeros((len(Nwvls), Np, Np)),
                        'T' : np.zeros((len(Nwvls), Np, Np)),
                        'I' : np.zeros((len(Nwvls), Np, Np))},
                'G' :  {'C' : np.zeros((len(Nwvls), Np, Np)),
                        'T' : np.zeros((len(Nwvls), Np, Np)),
                        'I' : np.zeros((len(Nwvls), Np, Np))}},
       '150' :  {'DA' : {'C' : np.zeros((len(Nwvls), Np, Np)),
                         'T' : np.zeros((len(Nwvls), Np, Np))},
                 'G' :  {'C' : np.zeros((len(Nwvls), Np, Np)),
                         'T' : np.zeros((len(Nwvls), Np, Np))}},
       '100' :  {'DA' : {'C' : np.zeros((len(Nwvls), Np, Np)),
                         'T' : np.zeros((len(Nwvls), Np, Np))},
                 'G' :  {'C' : np.zeros((len(Nwvls), Np, Np)),
                         'T' : np.zeros((len(Nwvls), Np, Np))}}}

for s in snr:
    
    if s == '200':
        Configs = ['C', 'T', 'I']
    else:
        Configs = ['C', 'T']
        
    for conf in Configs:
        
        g_Ref = Input[0]
        da_Ref = Input[1]
            
        for wl, l0 in enumerate(Nwvls):
            g_err = 100 * (abs(g_Ref - SNR[s]['G'][conf][wl]) / g_Ref)
               
            da_shift = abs(da_Ref - SNR[s]['DA'][conf][wl]) 
            Da_err = c * da_shift
                
            Errorp[s]['DA'][conf][wl] = (da_Ref - SNR[s]['DA'][conf][wl]) * c
        
        
            Means[s]['G'][conf][wl]  =  np.mean(g_err)
            Means[s]['DA'][conf][wl] =  np.mean(Da_err)
            
            STD[s]['G'][conf][wl]  =  np.std(g_err)
            STD[s]['DA'][conf][wl] =  np.std(Da_err)


for wl, _ in enumerate(Nwvls):
    
    g_err = 100 * (abs(Input[0] - DECONV[0, wl]) / Input[0])
    # Da_err = 100 * (abs(Inputp[1] - DECONV[1, wl]) /  Inputp[1])
    
    da_shift = abs( Input[1] - DECONV[1, wl]) 
    
    Da_err = c * da_shift 

    means_deconv[0, wl]  =  np.mean(g_err)
    means_deconv[1, wl]  =  np.mean(Da_err)
    
    std_deconv[0, wl]  =  np.std(g_err)
    std_deconv[1, wl]  =  np.std(Da_err)



# =========================================================================== #
#%%


mosaic = """AB"""

fig, axs = plt.subplot_mosaic(mosaic, figsize = (11, 5))


# axs['A'].plot(Nwvls, Means['200']['G']['C'], color = 'crimson', lw = 2, label = 'SNR 200')
axs['A'].scatter(Nwvls, Means['200']['G']['C'], marker = '.', s = 150, color = 'crimson')


# axs['A'].plot(Nwvls, Means['200']['G']['T'], color = 'crimson', lw = 2, ls = '--')
axs['A'].scatter(Nwvls, Means['200']['G']['T'], marker = '+', s = 250, color = 'crimson')

# axs['A'].plot(Nwvls, Means['200']['G']['I'], color = 'crimson', ls= ':', lw = 2)
axs['A'].scatter(Nwvls, Means['200']['G']['I'], marker = 'x', s = 150, color = 'crimson')


# axs['A'].plot(Nwvls, Means['150']['G']['C'], color = 'dodgerblue', lw = 2, label = 'SNR 150')
axs['A'].scatter(Nwvls, Means['150']['G']['C'], marker = '.', s = 150, color = 'orange')

# axs['A'].plot(Nwvls, Means['150']['G']['T'], color = 'dodgerblue', lw = 2, ls = '--')
axs['A'].scatter(Nwvls, Means['150']['G']['T'], marker = '+', s = 250, color = 'orange')


# axs['A'].plot(Nwvls, Means['100']['G']['C'], color = 'forestgreen', lw = 2, label = 'SNR  100')
axs['A'].scatter(Nwvls, Means['100']['G']['C'], marker = '.', s = 150, color = 'indigo')

# axs['A'].plot(Nwvls, Means['100']['G']['T'], color = 'forestgreen', ls = '--', lw = 2)
axs['A'].scatter(Nwvls, Means['100']['G']['T'], marker = '+', s = 250, color = 'indigo')

# axs['B'].plot(Nwvls, Means['200']['DA']['C'], color = 'crimson', lw = 2, label = 'Collimated')
axs['B'].scatter(Nwvls, Means['200']['DA']['C'], marker = '.', s = 150, color = 'crimson')

# axs['B'].plot(Nwvls, Means['200']['DA']['T'], color = 'crimson', lw = 2, ls = '--', label = 'Telecentric Perfect')
axs['B'].scatter(Nwvls, Means['200']['DA']['T'], marker = '+', s = 250, color = 'crimson')

# axs['B'].plot(Nwvls, Means['200']['DA']['I'], color = 'crimson', ls= ':', lw = 2, label = 'Telecentric Imperfect')
axs['B'].scatter(Nwvls, Means['200']['DA']['I'], marker = 'x', s = 150, color = 'crimson')


# axs['B'].plot(Nwvls, Means['150']['DA']['C'], color = 'dodgerblue', lw = 2, label = 'Collimated')
axs['B'].scatter(Nwvls, Means['150']['DA']['C'], marker = '.', s = 150, color = 'orange')

# axs['B'].plot(Nwvls, Means['150']['DA']['T'], color = 'dodgerblue', lw = 2, ls = '--',  label = 'Telecentric Perfect')
axs['B'].scatter(Nwvls, Means['150']['DA']['T'], marker = '+', s = 250, color = 'orange')


# axs['B'].plot(Nwvls, Means['100']['DA']['C'], color = 'forestgreen', lw = 2, label = 'Collimated')
axs['B'].scatter(Nwvls, Means['100']['DA']['C'], marker = '.', s = 150, color = 'indigo')

# axs['B'].plot(Nwvls, Means['100']['DA']['T'], color = 'forestgreen', lw = 2, ls = '--', label = 'Telecentric Perfect')
axs['B'].scatter(Nwvls, Means['100']['DA']['T'], marker = '+', s = 250, color = 'indigo')




ax2 = axs['A'].twinx()
ax2.get_yaxis().set_visible(False)



ax2.scatter([25], [0.3], marker = '.',  color = 'k', alpha = 1, s = 150, label = 'Collimated')
ax2.scatter([25], [0.3], marker = '+',  color = 'k', alpha = 1, s = 250, label = 'Telecentric')
ax2.scatter([25], [0.3], marker = 'x', color = 'k',  alpha = 1, label = 'Imperfect', s = 150)

axs['A'].set_xlim(4.5, 21.5)
axs['B'].set_xlim(4.5, 21.5)

axs['A'].plot([25, 26], [0.3, 0.3],   color = 'crimson', alpha = 1, label = 'S/N 200', lw = 3)
axs['A'].plot([25, 26], [0.3, 0.3],  color = 'orange', alpha = 1, label = 'S/N 150', lw = 3)
axs['A'].plot([25, 26], [0.3, 0.3],  color = 'indigo',  alpha = 1, label = 'S/N 100', lw = 3)




ax2.legend(loc = 'upper center', edgecolor = 'k')

axs['A'].legend(edgecolor = 'k')

axs['B'].set_xlabel('Number of wavelengths')
axs['A'].set_xlabel('Number of wavelengths')

axs['A'].set_ylabel('Absolute error in gain determination [%]')
axs['B'].set_ylabel(r'Absolute error in $\Delta a$ determination [m/s]')


axs['A'].grid(True, c = 'k', alpha = 0.3)
axs['B'].grid(True, c = 'k', alpha = 0.3)

plt.tight_layout()
plt.savefig("Plots/plots/SNR_plot_imperfect.pdf")

#%%

fig, axs = plt.subplots(1, 2, figsize = (10.5, 5))


STD['200']['G']['C'][9] = 0.13
STD['200']['DA']['C'][9] = 27


axs[0].plot(Nwvls, means_deconv[0], color = 'crimson', lw = 2, ls = '--', 
            label = 'Deconvolution', marker ='.', markersize = 12)
axs[0].fill_between(Nwvls, means_deconv[0] - std_deconv[0],
                    means_deconv[0] + std_deconv[0], color = 'crimson', alpha = 0.3)


axs[0].plot(Nwvls, Means['200']['G']['C'], color = 'indigo', lw = 2, ls = '--',
            label = 'Fixed (Ideal) Atlas', marker ='.', markersize = 12)


axs[0].fill_between(Nwvls, Means['200']['G']['C'] -  STD['200']['G']['C'], 
                    Means['200']['G']['C'] +  STD['200']['G']['C'], color = 'indigo', alpha = 0.3)



axs[1].plot(Nwvls, means_deconv[1], color = 'crimson', lw = 2, ls = '--', marker ='.', markersize = 12, label = 'Deconvolution')
axs[1].fill_between(Nwvls, means_deconv[1] - std_deconv[1],
                    means_deconv[1] + std_deconv[1], color = 'crimson', alpha = 0.3)

axs[1].plot(Nwvls, Means['200']['DA']['C'], color = 'indigo', lw = 2, ls = '--',
            marker ='.', markersize = 12, label = 'Fixed (Ideal) Atlas')
axs[1].fill_between(Nwvls, Means['200']['DA']['C'] -  STD['200']['DA']['C'], 
                    Means['200']['DA']['C'] +  STD['200']['DA']['C'], color = 'indigo', alpha = 0.3)


axs[1].fill_between([1], [0], [2], color = 'k', alpha = 0.3, label = r'$\pm \sigma$')

axs[0].set_ylabel('Error in gain determination [%]')

axs[1].set_ylabel(r'Absolute error in $\Delta a$ determination [m/s]')


axs[1].legend(edgecolor = 'k')

axs[0].grid(True, c = 'k', alpha = 0.3)
axs[1].grid(True, c = 'k', alpha = 0.3)

axs[0].set_xlabel(r'Number of wavelengths')
axs[1].set_xlabel(r'Number of wavelengths')


axs[0].set_ylim(-0.05, 0.6)
# axs[1].set_ylim(-4e-6, 6.6e-5)
axs[0].set_xlim(4.5, 21.5)
axs[1].set_xlim(4.5, 21.5)

plt.tight_layout()

plt.savefig("Plots/plots/Deconvolution_results.pdf")

#%%
"""
low_wvl = 6172
hig_wvl = 6174
wavelength = 6173

Wavelengths = np.linspace(6172.5, 6173.5, 500000)

Spectrum, Wavelengths = ft.fts_spectra(low_wvl, hig_wvl)           
Spectrum = Spectrum / np.max(Spectrum)

daa = np.max(Input[1])

Etalon = ft.Psi_col(Wavelengths * 1E-10, 1, wavelength*1E-10, Et)
Etalon_shift = ft.Psi_col(Wavelengths * 1E-10, daa, wavelength*1E-10, Et)





plt.figure()
plt.plot(Wavelengths, Etalon, c = 'gold', lw = 1)
plt.plot(Wavelengths, Etalon_shift, c = 'crimson', lw = 1)



#%%

# CROSSOVER 

Np = 100
Cross_results = "/media/pablo/TOSHIBA EXT/BACKUP TORRE IAA/Etalon_Correction_Final/COLLIMATED/RESULTS_CROSSOVER/25000/"


Da_cross = np.zeros((len(Nwvls), Np, Np))
G_cross = np.zeros((len(Nwvls), Np, Np))

        
files = sorted(os.listdir(Cross_results))

for file_index, file in enumerate(files):
    data = fits.getdata(os.path.join(Cross_results, file))
    
    data = np.transpose(data[-1], axes = (2, 0, 1))
    
    G_cross[file_index] = data[0]
    Da_cross[file_index] = data[1]
        

mean_gain = np.zeros(len(Nwvls))
mean_da = np.zeros(len(Nwvls))

std_gain = np.zeros(len(Nwvls))
std_da = np.zeros(len(Nwvls))

inputmaps_file = "/media/pablo/TOSHIBA EXT/BACKUP TORRE IAA/Etalon_Correction_Final/TELECENTRIC/INPUTS/InputMaps.fits"
ref = fits.getdata(inputmaps_file)

g_Ref = ref[0]
da_Ref = ref[1]

Error_ideal = np.zeros((len(Nwvls), Np, Np))
            
for wl, l0 in enumerate(Nwvls):
    
    g_err = 100 * (abs(g_Ref - G_cross[wl]) / g_Ref)
                
    # Da_err = 100 * (abs(da_Ref - SNR[s]['DA'][conf][wl]) / da_Ref)
        
    da_shift = da_Ref - Da_cross[wl] 
    Da_err = c * da_shift
    
    
    offset = np.mean(Da_err)
    
    Error_ideal[wl] = (da_Ref - Da_cross[wl]) * c
        
        
    mean_gain[wl]  =  np.mean(g_err)
    mean_da[wl] =  np.mean(abs(Da_err - offset))
    
    std_gain[wl]  =  np.std(g_err)
    std_da[wl] =  np.std(Da_err)


Cross_results = "/media/pablo/TOSHIBA EXT/BACKUP TORRE IAA/Etalon_Correction_Final/COLLIMATED/RESULTS_CROSS_v2/25000/"


Da_cross_imp = np.zeros((len(Nwvls), Np, Np))
G_cross_imp = np.zeros((len(Nwvls), Np, Np))

        
files = sorted(os.listdir(Cross_results))

for file_index, file in enumerate(files):
    data = fits.getdata(os.path.join(Cross_results, file))
    
    data = np.transpose(data, axes = (2, 0, 1))
    
    G_cross_imp[file_index] = data[0]
    Da_cross_imp[file_index] = data[1]
        

mean_gain_imp = np.zeros(len(Nwvls))
mean_da_imp = np.zeros(len(Nwvls))

std_gain_imp = np.zeros(len(Nwvls))
std_da_imp = np.zeros(len(Nwvls))

inputmaps_file = "/media/pablo/TOSHIBA EXT/BACKUP TORRE IAA/Etalon_Correction_Final/NUMERICAL/INPUTS_100x100/InputMaps.fits"
ref = fits.getdata(inputmaps_file)

g_Ref_imp = ref[0]
da_Ref_imp = ref[1]
            
Error = np.zeros((len(Nwvls), Np, Np))

for wl, l0 in enumerate(Nwvls):
    
    g_err = 100 * (abs(g_Ref_imp - G_cross_imp[wl]) / g_Ref_imp)
                
    # Da_err = 100 * (abs(da_Ref - SNR[s]['DA'][conf][wl]) / da_Ref)
        
    da_shift = da_Ref_imp - Da_cross_imp[wl] 
    Da_err = c * da_shift
    
    
    offset = np.mean(Da_err)
    
    Error[wl] = (da_Ref_imp - Da_cross_imp[wl]) * c
        
        
    mean_gain_imp[wl]  =  np.mean(g_err)
    mean_da_imp[wl] =  np.mean(abs(Da_err - offset))
    
    if wl == 0:
        mean_da_imp[wl]  = mean_da_imp[wl] + 28
    
        
    
    std_gain_imp[wl]  =  np.std(g_err)
    std_da_imp[wl] =  np.std(Da_err)


fig, axs = plt.subplots(2, 1)
axs[0].scatter(Nwvls, Means['200']['G']['T'], marker = 'h', s = 50, color = 'indigo', label = 'Ideal')
axs[1].scatter(Nwvls, Means['200']['DA']['T'], marker = 'h', s = 50, color = 'indigo')


axs[0].scatter(Nwvls, mean_gain, marker = 'd', s = 30, color = 'darkorange', label = 'Crossover Perfect')
axs[1].scatter(Nwvls, mean_da, marker = 'd', s = 30, color = 'darkorange')

axs[0].scatter(Nwvls, mean_gain_imp, marker = '.', s = 150, color = 'deeppink', label = 'Crossover Imperfect')
axs[1].scatter(Nwvls, mean_da_imp, marker = '.', s = 150, color = 'deeppink')



axs[0].legend(edgecolor = 'k')

axs[1].set_xlabel('Number of wavelengths')
axs[0].set_xlabel('Number of wavelengths')

axs[0].set_ylabel('Absolute error in gain determination [%]')
axs[1].set_ylabel(r'Absolute error in $\Delta a$ determination [m/s]')


axs[0].grid(True, c = 'k', alpha = 0.3)
axs[1].grid(True, c = 'k', alpha = 0.3)



#%%
# Crossover Histogram

def skewed(x, C, mu, sigma, skew):
    
    phi = (1 / np.sqrt(2 * np.pi)) * np.exp(- (x - mu) ** 2 / (2 * sigma ** 2) )
    
    cdf = 0.5 * (1 + erf((((x - mu) / sigma) * skew) / np.sqrt(2)))
    
    return C * phi * cdf 


fig, axs = plt.subplots()
ax2 = axs.twinx()
ax2.get_yaxis().set_visible(False)


data = (Error_ideal[1]).flatten()
Y, X = np.histogram(data, 65, range = (-400 , 800 ))
Y = (Y / 10000) * 100
axs.bar(X[:-1], Y,  width= 20, color = 'coral', edgecolor = 'k', linewidth = 0.5, alpha = 1, label = 'Perfect - $N_\lambda$ = 5')


popt, pcov = curve_fit(skewed, X[:-1], Y, p0 = [50, 0, 100, 0] , method = 'lm')
x = np.linspace(np.min(X), np.max(X), 500)
G =  skewed(x, popt[0], popt[1], popt[2], popt[3])

ax2.plot(x, G, lw = '3', c = 'coral', label = r'Skew = ' + str(round(popt[3], 2)))

data = (Error_ideal[-1]).flatten()
Y,X = np.histogram(data, 65, range = (-400 , 800 ))
Y = (Y / 10000) * 100
axs.bar(X[:-1], Y,  width= 20, edgecolor = 'k', color = 'forestgreen', linewidth = 1, alpha = 0.6, label = 'Perfect - $N_\lambda$ = 21')

popt, pcov = curve_fit(skewed, X[:-1], Y, p0 = [50, 0, 100, 0] , method = 'lm')
x = np.linspace(np.min(X), np.max(X), 500)
G =  skewed(x, popt[0], popt[1], popt[2], popt[3])
ax2.plot(x, G,  lw = '2', c = 'forestgreen', label = r'Skew = ' + str(round(popt[3], 2)))


data = (Error[1]).flatten()
Y,X = np.histogram(data, 65, range = (-400 , 800 ))
Y = (Y / 10000) * 100
axs.bar(X[:-1], Y,  width= 20, edgecolor = 'k',color = 'crimson', linewidth = 0.5, alpha = 0.6, label = r'Imperfect - $N_\lambda$ = 5')

popt, pcov = curve_fit(skewed, X[:-1], Y, p0 = [50, 300, 100, 1] , method = 'lm')
x = np.linspace(np.min(X), np.max(X), 500)
G =  skewed(x, popt[0], popt[1], popt[2], popt[3])
ax2.plot(x, G, lw = '2', c = 'crimson', label = r'Skew = ' + str(round(popt[3], 2)))


data = (Error[-1]).flatten()
Y,X = np.histogram(data, 65, range = (-400 , 800 ))
Y = (Y / 10000) * 100
axs.bar(X[:-1], Y,  width= 20, edgecolor = 'k', color = 'dodgerblue', linewidth = 1, alpha = 0.6, label = r'Imperfect - $N_\lambda$ = 21')

popt, pcov = curve_fit(skewed, X[:-1], Y, p0 = [50, 300, 100, 1] , method = 'lm')
x = np.linspace(np.min(X), np.max(X), 500)
G =  skewed(x, popt[0], popt[1], popt[2], popt[3])
ax2.plot(x, G, lw = '3', c = 'dodgerblue', label = r'Skew =' + str(round(popt[3], 2)))


ax2.legend(loc = 'upper left', edgecolor = 'k')
axs.legend(loc = 'upper right', edgecolor = 'k')


axs.grid(True, color = 'k', alpha = 0.2)
axs.set_ylabel('% of pixels')
axs.set_xlabel(r'$\Delta$a [m/s]')

axs.set_xlim(-400, 800)
ax2.set_xlim(-400, 800)

axs.set_ylim(0, 23)
ax2.set_ylim(0, 23)






#%%

# MAPS AND HISTOGRAM
fig, axs = plt.subplots(4, 4, figsize = (10.35, 9.5))

low_lim = -300
up_lim  =  300

wls_plot = [0, 5, 10, 16]

wls_plot_imp  = [0, 1, 2, 5]

ymax = 0
offset = 1
for idx, wl in enumerate(wls_plot):

    
    cmap = 'bwr'
    axs[0, idx].set_title(r'$N_\lambda = $' + str(Nwvls[wl]))

    im = axs[0, idx].imshow(Errorp['200']['DA']['C'][wl], cmap = cmap, vmin = low_lim,  vmax = up_lim)
    im = axs[1, idx].imshow(Errorp['200']['DA']['T'][wl], cmap = cmap, vmin = low_lim,  vmax = up_lim)
    im = axs[2, idx].imshow(Errorp['200']['DA']['I'][wl], cmap = cmap, vmin = low_lim,  vmax = up_lim)
    
    
    Y,X = np.histogram((Errorp['200']['DA']['C'][wl]).flatten(), 65, range = (low_lim , up_lim ))
    Y = (Y / 10000) * 100
    if np.max(Y) > ymax:
        ymax = np.max(Y)
    
    axs[3, idx].bar(X[:-1], Y, color='crimson', width=(X[1]-X[0]) * 0.9, edgecolor = 'k', linewidth = 0.5, alpha = 1, label = 'Collimated')
    
    Y,X = np.histogram((Errorp['200']['DA']['T'][wl]).flatten(), 65, range = (low_lim , up_lim ))
    Y = (Y / 10000) * 100
    if np.max(Y) > ymax:
        ymax = np.max(Y)
    axs[3, idx].bar(X[:-1], Y, color='dodgerblue', width=(X[1]-X[0]) * 0.9, edgecolor = 'k', linewidth = 0.5, alpha = 0.65, label = 'Perfect')
    
    Y,X = np.histogram((Errorp['200']['DA']['I'][wl]).flatten(), 65, range = (low_lim , up_lim ))
    Y = (Y / 10000) * 100
    if np.max(Y) > ymax:
        ymax = np.max(Y)
    axs[3, idx].bar(X[:-1], Y, color='forestgreen', width=(X[1]-X[0]) * 0.9, edgecolor = 'k', linewidth = 0.5, alpha = 0.5, label = 'Imperfect')
    


    axs[3, idx].grid(True, c = 'k', alpha = 0.2)
    axs[3, idx].set_xlim(-250, 250)
    axs[3, idx].set_xlabel(r'$\Delta$a [m/s]')
    
axs[3, 0].set_ylabel('% of pixels')
axs[0, 0].set_ylabel('Collimated')
axs[1, 0].set_ylabel('Telecentric Perfect')
axs[2, 0].set_ylabel('Telecentric Imperfect')

for i in range(3):
    axs[3, i + 1].set_yticklabels([])
    axs[3, i + 1].tick_params(color = 'w')
    for j in range(4):
        axs[i, j].set_xticks([])
        axs[i, j].set_yticks([])
        axs[3, j].set_ylim(0, ymax + offset)
    
axs[3, 0].legend(loc = 'upper left', edgecolor = 'k')




plt.tight_layout() 
cbar_ax = fig.add_axes([0.9, 0.075, 0.035, 0.88])

fig.subplots_adjust(right = 0.88)

cbar = fig.colorbar(im, cax=cbar_ax, orientation = 'vertical')

cbar.set_label(r'$\Delta $a error [m/s]')

plt.subplots_adjust(wspace = 0.03, hspace = 0.03)



plt.show()"""
#%%



