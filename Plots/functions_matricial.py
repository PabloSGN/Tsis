# ============================ IMPORTS ====================================== #
 
import os
import numpy as np
import matplotlib.pyplot as plt

import scipy.signal as sgn

from scipy.integrate import simps
 

pi = np.pi
# ============================================================================ #

def fts_spectra(wli,wlf):

    """'Kitt Peak FTS-Spectral-Atlas'"""
    # print('Enter end wavelength (3290 - 12508 A)')
    #path = '/Users/dorozco/Dropbox (IdAdA)/Python/'
    file = 'Plots/fts.npz'
    # np.savez('fts.npz', fts=fts, fts_w=fts_w)
    data = np.load(file)
    fts = data['fts']
    fts_w = data['fts_w']

    indx = np.where((fts_w > int(wli)) & (fts_w < int(wlf)))

    return fts[indx], fts_w[indx]

def Finesse(R):
    return 4 * R / (1 - R) ** 2

def a_tun(Et, l0):
    
    m  = round(2 * Et['n'] * Et['d'] / l0)                    # Order of the resonance peak
    dh = (m * l0 - 2 * Et['n'] * Et['d']) / (2 * Et['n'])     # Variation of thickness to tune again to wvl0
    thick = Et['d'] + dh

    return thick * Et['n'] * Et['cos_th']
    
def Psi_single(wvls, Da, l0, Et):
    
    a = a_tun(Et, l0)
    
    da_lambd = Da / wvls
    
    return 1 / (1 + Finesse(Et['R']) * np.sin(2 * pi * a * da_lambd) ** 2)

def Psi(wvls, Da, l0, Et):
    
    """
    Wvls -> Nwavelengths x 1
    Da   -> Ncoords x 1
    """
    
    a = a_tun(Et, l0)
    
    da_lambd = Da @ (1 / wvls).T
    
    return 1 / (1 + Finesse(Et['R']) * np.sin(2 * pi * a * da_lambd) ** 2)


def dPsi(wvls, Da, l0, Et):
    
    """
    Wvls -> Nwavelengths x 1
    Da   -> Ncoords x 1
    """
    
    a = a_tun(Et, l0)
    
    da_lambd = Da @ (1 / wvls).T
    
    a = a_tun(Et, l0)
    R = Et['R']
    
    num = 8 * pi * a * R * (-1 + R) ** 2 * np.sin(4 * pi * a * da_lambd)
    den = wvls[:, None].T * ((1 + R ** 2) - 2 * R * np.cos(4 * pi * a * da_lambd)) ** 2
    
    return - num / den[0]


def Profile(scan_wls, etalon_wls, g, da, phi_real, Et, Ncont):
    
    N_coords = len(da)
    N_wvls = len(scan_wls)
    Prof = np.zeros((N_coords, N_wvls))
    
    et_wvls_vector = etalon_wls[:, np.newaxis]
    Da_vector = da[:, np.newaxis]  
    
    for i, wvl in enumerate(scan_wls):
        
        Etalon = Psi(et_wvls_vector * 1E-10, Da_vector, wvl * 1E-10, Et)
        Etalon_Object = phi_real * Etalon
        Integral = simps(Etalon_Object, x = etalon_wls, axis = -1)
        
        Prof[:, i] = Integral
        
    Prof = Prof / Prof[:, Ncont, None]   
    
    return g[:, None] * Prof

def Profile_1p(scan_wls, etalon_wls, g, da, phi_real, Et, Ncont):
    
    N_wvls = len(scan_wls)
    Prof = np.zeros(N_wvls)
    
    for i, wvl in enumerate(scan_wls):
        
        Etalon = Psi_single(etalon_wls * 1E-10, da, wvl * 1E-10, Et)
        Etalon_Object = phi_real * Etalon
        Integral = simps(Etalon_Object, x = etalon_wls)
        Prof[i] = Integral
        
    Prof = Prof / Prof[Ncont]   
    
    return g * Prof
    

def function(phi_obs, scan_wls, etalon_wls, g, da, phi_real, Et, Ncont):
    
    Scan = Profile(scan_wls, etalon_wls, g, da, phi_real, Et, Ncont)
        
    return phi_obs - Scan
   
    
def df_a(scan_wls, etalon_wls, g, da, phi_real, Et, Ncont):
    
    N_coords = len(da)
    N_wvls = len(scan_wls)
    DF = np.zeros((N_coords, N_wvls))

    et_wvls_vector = etalon_wls[:, np.newaxis]
    Da_vector = da[:, np.newaxis]  
    
    # Normalization and its derivative
    Etalon  = Psi(et_wvls_vector * 1E-10, Da_vector, scan_wls[Ncont] * 1E-10, Et)
    Norm = simps(Etalon * phi_real, x = etalon_wls, axis = -1)
    
    dEtalon = dPsi(et_wvls_vector * 1E-10, Da_vector, scan_wls[Ncont] * 1E-10, Et)
    dNorm = simps(dEtalon * phi_real, x = etalon_wls, axis = -1)
    
    for i, wvl in enumerate(scan_wls):
        Etalon  =  Psi(et_wvls_vector * 1E-10, Da_vector, wvl * 1E-10, Et)
        dEtalon = dPsi(et_wvls_vector * 1E-10, Da_vector, wvl * 1E-10, Et)
        
        Etalon_Object  = phi_real * Etalon
        dEtalon_Object = phi_real * dEtalon
        
        Integral  = simps(Etalon_Object,  x = etalon_wls, axis = -1)
        dIntegral = simps(dEtalon_Object, x = etalon_wls, axis = -1)
        
        DF[:, i] =  (dIntegral * Norm - Integral * dNorm) / Norm ** 2
        
    return - g[:, None] * DF
        
    
def df_g(scan_wls, etalon_wls, da, phi_real, Et, Ncont):
    
    N_coords = len(da)
    N_wvls = len(scan_wls)
    DF = np.zeros((N_coords, N_wvls))
    
    et_wvls_vector = etalon_wls[:, np.newaxis]
    Da_vector = da[:, np.newaxis]  
    
    for i, wvl in enumerate(scan_wls):
        
        Etalon = Psi(et_wvls_vector * 1E-10, Da_vector, wvl * 1E-10, Et)
        Etalon_Object = phi_real * Etalon
        Integral = simps(Etalon_Object, x = etalon_wls, axis = -1)
        
        DF[:, i] = Integral
        
    DF = DF / DF[:, -1, None]   
    
    return - DF
    
def Jacobian(scan_wls, etalon_wls, g, da, phi_real, Et, Ncont):
    
    N_coords = len(da)
    N_wvls = len(scan_wls)
    
    J = np.zeros((N_coords, N_wvls, 2))
    
    J[:, :, 0] =  df_g(scan_wls, etalon_wls, da, phi_real, Et, Ncont)
    J[:, :, 1] =  df_a(scan_wls, etalon_wls, g, da, phi_real, Et, Ncont)
    
    return J 

def invert(H):

    try:
        # u, s, v = np.linalg.svd(H)
        # Ainv = (u * s[:, None, :]) @ v
        Ainv = np.linalg.inv(H)
    except:
        
        Ainv = np.zeros((3, 3))
    
    return Ainv

def Newton(phi_obs, scan_wls, etalon_wls, phi_real, x0, Et, Ncont, iters, mode = 'normal'):
    
    x = x0
    
    if mode == 'full':
        Output = np.zeros((iters, np.shape(x)[0], 2))
    
    for i in range(iters):
        # print('iter', i)
        J = Jacobian(scan_wls, etalon_wls, x[:, 0], x[:, 1], phi_real, Et, Ncont)
        F = function(phi_obs, scan_wls, etalon_wls, x[:, 0], x[:, 1], phi_real, Et, Ncont)
        
        F = F[:, :, np.newaxis]
        Jt = np.transpose(J, axes = (0, 2, 1))

        H = Jt @ J 

        HI = invert(H)

        JF = Jt @ F

        Change = HI @ JF
        Change = np.squeeze(Change)
                
        x = x - Change
        
        if mode == 'full':
            Output[i] = x
    
    if mode == 'full':
        return Output
    else:
        return x

#%%

