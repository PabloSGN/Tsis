# ============================ IMPORTS ====================================== #
 
import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import simps

import telecentric_functions as tl


# ============================= CONFIG ====================================== #
 
pi = np.pi

# =========================================================================== #

def fts_spectra(wli,wlf):

    """'Kitt Peak FTS-Spectral-Atlas'"""
    # print('Enter end wavelength (3290 - 12508 A)')
    #path = '/Users/dorozco/Dropbox (IdAdA)/Python/'
    file = 'fts.npz'
    # np.savez('fts.npz', fts=fts, fts_w=fts_w)
    data = np.load(file)
    fts = data['fts']
    fts_w = data['fts_w']

    indx = np.where((fts_w > int(wli)) & (fts_w < int(wlf)))

    return fts[indx], fts_w[indx]

def a_tuning(Et, l0):
    
    wvl = l0+l0/(16*Et['n']**2*Et['fnum']**2) #Peak wavelength
    m  = round(2 * Et['n'] * Et['d'] * np.cos(Et['theta']) / l0)                               #Order of the resonance peak
    dh=(m*wvl-2*Et['n']*Et['d'])/(2*Et['n']) #Variation of thickness to tune again to wvl0    #Variation of thickness to tune again to wvl0
    thick = Et['d'] + dh 
    return thick * Et['n'] 
    
def Psi_single(wvl, Da, l0, Et):

    R = Et['R']
    
    F = (4 * R) / ((1 - R) ** 2)
        
    nh = a_tuning(Et, l0)
    
    a = (2 * pi * nh) * Da / wvl # Creating Ncoord x Nwavelengths array
        
    f = 2 * Et['fnum']
    b = 1/(np.sqrt(2) * Et['n'] * f)
    
    alpha1 = tl.alpha1(a, b, F)
    alpha2 = tl.alpha2(a, b, F)
    gamma1 = tl.gamma1(a, F)
    gamma2 = tl.gamma2(a, b, F)
    gamma3 = tl.gamma3(F)
    gamma4 = tl.gamma4(a, b, F)
    gamma5 = tl.gamma5(a, F)
        
    Re = tl.RealE(alpha1, gamma1, gamma2)
    Im = tl.imE(R, alpha2, gamma3, gamma4, gamma5)
    
    return Re ** 2 + Im ** 2


def Psi(wvl, Da, l0, Et):
    
    """
    Wvls -> Nwavelengths x 1
    Da   -> Ncoords x 1
    """

    R = Et['R']
    
    F = (4 * R) / ((1 - R) ** 2)
        
    nh = a_tuning(Et, l0)
    
    a = (2 * pi * nh) * Da @ (1 / wvl).T # Creating Ncoord x Nwavelengths array
        
    f = 2 * Et['fnum']
    b = 1/(np.sqrt(2) * Et['n'] * f)
    
    alpha1 = tl.alpha1(a, b, F)
    alpha2 = tl.alpha2(a, b, F)
    gamma1 = tl.gamma1(a, F)
    gamma2 = tl.gamma2(a, b, F)
    gamma3 = tl.gamma3(F)
    gamma4 = tl.gamma4(a, b, F)
    gamma5 = tl.gamma5(a, F)
        
    Re = tl.RealE(alpha1, gamma1, gamma2)
    Im = tl.imE(R, alpha2, gamma3, gamma4, gamma5)
    
    return Re ** 2 + Im ** 2

def dPsi(wvl, Da, l0, Et):

    """
    Wvls -> Nwavelengths x 1
    Da   -> Ncoords x 1
    """    

    R = Et['R']
    
    F = (4 * R) / ((1 - R) ** 2)
        
    nh = a_tuning(Et, l0)
    
    a = (2 * pi * nh) * Da @ (1 / wvl).T
        
    f = 2 * Et['fnum']
    b = 1/(np.sqrt(2) * Et['n'] * f)
    
    alpha1 = tl.alpha1(a, b, F)
    alpha2 = tl.alpha2(a, b, F)
    gamma1 = tl.gamma1(a, F)
    gamma2 = tl.gamma2(a, b, F)
    gamma3 = tl.gamma3(F)
    gamma4 = tl.gamma4(a, b, F)
    gamma5 = tl.gamma5(a, F)
            
    Re = tl.RealE(alpha1, gamma1, gamma2)
    Im = tl.imE(R, alpha2, gamma3, gamma4, gamma5)
    
    der_Re = tl.der_RealE(a, b, F, alpha1, gamma1, gamma2)
    der_Im = tl.der_ImE(R, a, b, F, alpha2, gamma3, gamma4, gamma5)
    
    da = 2 * pi * nh * (1 / wvl[:, None]).T
    
    return (2 * (Re * der_Re + Im * der_Im) * da)[0]


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
        
    Prof = Prof / Prof[:, -1, None]   
    
    return g[:, None] * Prof

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

