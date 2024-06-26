#%%
# ============================ IMPORTS ====================================== #
 
import time
import numpy as np

# Integrals
from scipy.integrate import simps, nquad
from scipy.interpolate import RegularGridInterpolator

import general_func as gf

# Own Libs
import etalon_funtions as etf

# ============================= CONFIG ====================================== #


pi = np.pi

#Import array with spectral transmission for different incidence angles
namefits='transmission2D'
data=gf.read_fits(f"Plots/{namefits}.fits")

# Wavelength vector
wavelength = 617.334388e-9
low_wvl = wavelength-0.15e-9
high_wvl = wavelength+0.15e-9
Nwvls=data.shape[0]
wvl_vector= np.linspace(low_wvl, high_wvl, Nwvls)

#Chief ray angle vector
Ntheta=data.shape[1]
theta_max=1 #Maximum deviation of chief ray
theta_vector=np.linspace(0,theta_max,Ntheta)

#Interpolation
interp = RegularGridInterpolator((wvl_vector, theta_vector), data,\
bounds_error=False, fill_value=np.min(data), method='linear')

# ============================================================================ #


def fts_spectra(wli : float, wlf : float):
    """
    Function to generate the object and xaxis for integrals

    Parameters
    ----------
    wli : float
        Lower limit of xaxis.
    wlf : float
        Upper limit of xaxis.

    Returns
    -------
    Object
        array containing the object.
    Xaxis
        array containing the xaxis.

    """
    """'Kitt Peak FTS-Spectral-Atlas'"""

    file = 'Plots/fts.npz'
    data = np.load(file)
    fts = data['fts']
    fts_w = data['fts_w']

    indx = np.where((fts_w > int(wli)) & (fts_w < int(wlf)))

    return fts[indx], fts_w[indx]


# --------------------------------------------------------------------------- #

# MAIN FUNCTION.
def da_tuning(Et, l0):
    
    wvl = l0+l0/(16*Et['n']**2*Et['fnum']**2) #Peak wavelength
    m  = round(2 * Et['n'] * Et['d']  / l0)                               #Order of the resonance peak
    dh=(m*wvl-2*Et['n']*Et['d'])/(2*Et['n']) #Variation of thickness to tune again to wvl0    #Variation of thickness to tune again to wvl0
    thick = Et['d'] + dh 
    return (thick * Et['n']) 

def Psi(wvls, l0, Da, da_ref, Et, theta):
   
    #da_ref = da_tuning(Et, wavelength * 1E-10)
    da_tun = da_tuning(Et, l0)
    da_mod = (da_tun / da_ref) * Da

    wvls = wvls * np.abs(2 - da_mod) # To shift correctly the profile

    prof = interp((wvls, theta))

    return prof


# --------------------------------------------------------------------------- #

def dPsi(wvls, l0, Da, da_ref, Et, theta, h = 5e-8):

    f_plus  = Psi(wvls, l0,  Da + h, da_ref, Et, theta = theta)
    f_minus = Psi(wvls, l0,  Da - h, da_ref, Et, theta = theta)
        
    return (f_plus - f_minus) / (2 * h)  

def dPsi_theta(wvls, l0, Da, Et, da_ref, theta, h = 1e-7):

    f_plus  = Psi(wvls, l0,  Da, da_ref, Et, theta = theta + h)
    f_minus = Psi(wvls, l0,  Da, da_ref, Et, theta = theta - h)
        
    return (f_plus - f_minus) / (2 * h)  

# --------------------------------------------------------------------------- #

def Profile_1p(scan_wls, etalon_wls, Profiles, g, phi_real, Ncont):
    """
    Compute integrals of Object x Psi and normalization.
    """
    
    N_wvls = len(scan_wls)
    Prof = np.zeros(N_wvls)
        
    for i, _ in enumerate(scan_wls):
        
        Etalon = Profiles[i]
        
        Integral = simps(phi_real * Etalon, x = etalon_wls)
        
        Prof[i] = Integral
                
    Prof = Prof / Prof[Ncont]   

    return g * Prof


def function_1p(phi_obs, scan_wls, etalon_wls, Profiles, g, phi_real, Ncont):
    
    """
    Compute function.
    """
    
    Scan = Profile_1p(scan_wls, etalon_wls, Profiles, g, phi_real, Ncont)
        
    return phi_obs - Scan


def df_a_1p(scan_wls, etalon_wls, Profiles, der_Profiles, g, phi_real, Ncont):
    
    """
    Derivative of the function with respect to a
    """
 
    N_wvls = len(scan_wls)
    DF = np.zeros(N_wvls)
        
    # Normalization and its derivative
    Etalon  = Profiles[Ncont]
    Norm = simps(Etalon * phi_real, x = etalon_wls, axis = -1)
    
    dEtalon = der_Profiles[Ncont]
    dNorm = simps(dEtalon * phi_real, x = etalon_wls, axis = -1)

    for i, _ in enumerate(scan_wls):

        Etalon  = Profiles[i]
        dEtalon = der_Profiles[i]
        
        Etalon_Object  = phi_real * Etalon
        dEtalon_Object = phi_real * dEtalon
        
        Integral  = simps(Etalon_Object,  x = etalon_wls, axis = -1)
        dIntegral = simps(dEtalon_Object, x = etalon_wls, axis = -1)
        
        DF[i] =  (dIntegral * Norm - Integral * dNorm) / Norm ** 2
        
    return - g * DF

# --------------------------------------------------------------------------- #

def df_g_1p(scan_wls, etalon_wls, Profiles, phi_real, Ncont):
    
    """
    Derivative of the function with respect to g
    """

    N_wvls = len(scan_wls)
    DF = np.zeros(N_wvls)
    
    for i, _ in enumerate(scan_wls):
                
        Etalon_Object = phi_real * Profiles[i]
        Integral = simps(Etalon_Object, x = etalon_wls, axis = -1)
        DF[i] = Integral
        
    DF = DF / DF[Ncont]
    
    return - DF

# --------------------------------------------------------------------------- #

def Jacobian_1p(scan_wls, etalon_wls, Profiles, der_Profiles, g, phi_real, Ncont):
    
    """
    Compute Jacobian (derivatives of function with respect of both params)
    """
    
    N_wvls = len(scan_wls)
    
    J = np.zeros((N_wvls, 2))
    
    J[:, 0] =  df_g_1p(scan_wls, etalon_wls, Profiles, phi_real, Ncont)
    J[:, 1] =  df_a_1p(scan_wls, etalon_wls, Profiles, der_Profiles, g, phi_real, Ncont)
    
    return J
 
def invert(H):
    
    """
    Invert matrix
    """

    try:
        Ainv = np.linalg.inv(H)
    except:
        print('Unable to invert')
        Ainv = np.zeros((2, 2))
    
    return Ainv

# --------------------------------------------------------------------------- #

def Newton_1p(phi_obs, scan_wls, etalon_wls, phi_real, x0, Ncont, 
              iters, Et, mode = "normal", fit = "Da"):
    
    # Some info of the input params
    Nwvls = len(scan_wls)
    N_etalon_wls = len(etalon_wls) 
    x = np.zeros(2)

    x[0] = x0[0]

    if fit == "Da":
        x[1] = x0[1]
    elif fit == "Th":
        x[1] = x0[2] 
    
    da = x0[1]
    th = x0[2]

    if mode == 'full':
        Output = np.zeros((iters, np.shape(x)[0], 2))

    # For loop with iterations for the newtons's method    
    for it in range(iters):

        if it > 0:
            if fit == "Da":
                da = x[1]
            elif fit == "Th":
                if 0 <= x[1] < 1:
                    th = x[1]
                elif x[1] < 0:
                    th = 0.05
                    x[1] = th
                    print('Underestimated')
                else:
                    th = th + 0.2
                    x[1] = th
                    print('Overestimated')

        #print(f"Th : {th}")

        # Compute the profiles in the first place to avoid repetitions.

        # Start by computing etalon profiles. 
        tic = time.time()
        Counter = 0
        Profiles = np.zeros((Nwvls, N_etalon_wls))

        dr = da_tuning(Et, wavelength)

        for i_wvls, wvl in enumerate(scan_wls):
            if fit == "Da":
                Profiles[i_wvls] = Psi(wvls = etalon_wls * 1E-10, l0 = wvl * 1E-10, 
                                    Da = da, da_ref=dr, Et = Et, theta = th) 
            elif fit == "Th":
                Profiles[i_wvls] = Psi(wvls = etalon_wls * 1E-10, l0 = wvl * 1E-10, 
                                    Da = da, da_ref=dr, Et = Et, theta = th) 

            Counter += 1

        # And derivatives. 
        der_Profiles = np.zeros((Nwvls, N_etalon_wls))
        for i_wvls, wvl in enumerate(scan_wls):

            if fit == "Da":
                der_Profiles[i_wvls] = dPsi(wvls = etalon_wls * 1E-10, l0 = wvl * 1E-10, Da = da, 
                                            da_ref = dr, Et = Et, theta = th) 
            elif fit == "Th":
                der_Profiles[i_wvls] = dPsi_theta(wvls = etalon_wls * 1E-10, l0 = wvl * 1E-10, Da = da, 
                                                  da_ref = dr, Et = Et, theta = th) 
            else:
                raise Exception('Please enter a valid mode for "fit" (Da or Th)')

            Counter += 1

        tac = time.time()
        
        # Jacobian
        J = Jacobian_1p(scan_wls, etalon_wls, Profiles, der_Profiles, x[0], phi_real, Ncont)
        
        # Function value
        F = function_1p(phi_obs, scan_wls, etalon_wls, Profiles, x[0], phi_real, Ncont)
        
        # Application of Newtons method to derive next step of x

        Jt = np.transpose(J, axes = (1, 0))
        H = Jt @ J 
        HI = invert(H)
        JF = Jt @ F

        Change = HI @ JF
        Change = np.squeeze(Change)         

        #print(f"SUM:{np.sum(F)}")     
        
        #print(f"Step:  {x[1]} - {Change[1]} = {x[1] - Change[1]}")
        x = x - Change

        if mode == 'full':
            Output[it] = x
    
    if mode == 'full':
        return Output
    else:
        return x



