# ============================ IMPORTS ====================================== #
 
import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import simps, nquad
 
import telecentric_functions as tl
import etalon_funtions as etf

# ============================= CONFIG ====================================== #
 
plt.style.use('dark_background')
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

# ------------------------- ETALON FUNCTIONS -------------------------------- #

def Finesse(R):
    return 4 * R / (1 - R) ** 2

def a_tun_collimated(Et, l0):
    
    m  = round(2 * Et['n'] * Et['d'] / l0)                    # Order of the resonance peak
    dh = (m * l0 - 2 * Et['n'] * Et['d']) / (2 * Et['n'])     # Variation of thickness to tune again to wvl0
    thick = Et['d'] + dh

    return thick * Et['n'] * Et['cos_th']

def a_tuning(Et, l0):
    
    wvl = l0+l0/(16*Et['n']**2*Et['fnum']**2) #Peak wavelength
    m  = round(2 * Et['n'] * Et['d'] * np.cos(Et['theta']) / l0)                               #Order of the resonance peak
    dh=(m*wvl-2*Et['n']*Et['d'])/(2*Et['n']) #Variation of thickness to tune again to wvl0    #Variation of thickness to tune again to wvl0
    thick = Et['d'] + dh 
    return thick * Et['n'] 


def Psi_col(wvl, Da, l0, Et):
    
    a = a_tun_collimated(Et, l0)
    
    return 1 / (1 + Finesse(Et['R']) * np.sin(2 * pi * a * Da / wvl) ** 2)


def dPsi_R(wvl, Da, l0, Et, R):
    
    a = a_tun_collimated(Et, l0)
    A = 2 * pi * a * Da / wvl
    
    numerator = 4 * (R ** 2 - 1) * np.sin(A) 
    denominator = (4 * R * np.sin(A) + (R - 1) ** 2) ** 2

    return numerator / denominator

def dPsi_c(wvl, Da, l0, Et):
    
    a = a_tun_collimated(Et, l0)
    R = Et['R']
    
    num = 8 * pi * a * R * (-1 + R) ** 2 * np.sin(4 * pi * a * Da / wvl)
    den = wvl * ((1 + R ** 2) - 2 * R * np.cos(4 * pi * a * Da / wvl)) ** 2
    
    return - num / den

def Psi_tel(wvl, Da, l0, Et):
    
    R = Et['R']
    
    F = (4 * R) / ((1 - R) ** 2)
        
    nh = a_tuning(Et, l0)
    
    a = (2 * pi * nh) * Da  * (1 / wvl) # Creating Ncoord x Nwavelengths array
        
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


def Numerical_TelecentricEtalon(wvls, l0, Da, Et):

    Nl = int(len(wvls))
    Inum = np.zeros(Nl)

    tau  = etf.transm(Et['R'])
    xi   = np.sin(np.deg2rad(Et['angle1'])) * Et['fnum'] 
    eta  = np.sin(np.deg2rad(Et['angle2'])) * Et['fnum'] 
    xi0  = xi
    eta0 = eta
    theta3 = 0
    f = Et['fnum'] * 2 
    lims = [[0, 1], [0, 2 * np.pi]]  # Integral limits
    accur = {'limit': 50, 'epsabs':0,'epsrel': 1.49e-8}  # Optional arguments
    j = -1

    F= etf.F(Et['R'])

    for wvli in wvls:
        
        j += 1
        k = 2 * np.pi / wvli
        
        params = (tau, Et['R'], F, xi, eta, xi0, eta0, k, f, wvli, theta3, Et['n'], Et['fnum'], Et['d'], Da, l0)
        H11tr = nquad(etf.H11pintr_, lims, args=params, opts=accur)
        H11ti = nquad(etf.H11pinti_, lims, args=params, opts=accur)
        H11t = (H11tr[0] + 1j * H11ti[0]) / (np.pi * 1 ** 2)
        Inum[j] = np.real(H11t * np.conj(H11t))

    return Inum
# --------------------------- INTEGRAL FUNCTIONS ---------------------------- #

def Int_col(phi_real, central_wvl, etalon_wls, Da, Et):
    
    I = simps(phi_real * Psi_col(etalon_wls * 1E-10, Da, central_wvl * 1E-10, Et), 
              x = etalon_wls)
    
    return I

def Prof_col(phi_real, scan_wls, etalon_wls, g, Da, Et, Ncont):
    
    Scan = [Int_col(phi_real, wl, etalon_wls, Da, Et) for wl in scan_wls] 
    Scan = np.array(Scan)
    
    Norm = Int_col(phi_real, scan_wls[Ncont], etalon_wls, Da, Et)
    
    Scan /= Norm
    
    return g * Scan

def fun_col(phi_obs, scan_wls, etalon_wls, g, da, phi_real, Et, Ncont):
    
    Scan = Prof_col(phi_real, scan_wls, etalon_wls, g, da, Et, Ncont)
        
    return phi_obs - Scan



def df_g_col(phi_obs, scan_wls, etalon_wls, g, da, phi_real, Et, Ncont):
    
    Scan = [Int_col(phi_real, wl, etalon_wls, da, Et) for wl in scan_wls] 
    Scan = np.array(Scan)
    
    Norm = Int_col(phi_real, scan_wls[Ncont], etalon_wls, da, Et)
    
    return - (Scan / Norm)

def df_a_col(phi_obs, scan_wls, etalon_wls, g, da, phi_real, Et, Ncont):
    
    Int_ls = np.array([Int_col(phi_real, wl, etalon_wls, da, Et) for wl in scan_wls])
    Norm = Int_col(phi_real, scan_wls[Ncont], etalon_wls, da, Et)
    
    Int_ls_der = np.array([simps(phi_real * dPsi_c(etalon_wls * 1E-10, da, wl * 1E-10, Et), 
                  x = etalon_wls) for wl in scan_wls])
    
    Norm_der = simps(phi_real * dPsi_c(etalon_wls * 1E-10, da, scan_wls[Ncont] * 1E-10, Et), 
                  x = etalon_wls)
    
    return - g * (Int_ls_der * Norm  - Int_ls * Norm_der) / Norm ** 2

# ------------------------------- NEWTON ------------------------------------ #

def svdsolve(A):
    u, s, v = np.linalg.svd(A)
    
    Ainv = np.dot(v.transpose(), np.dot(np.diag(s**-1), u.transpose()))
    
    return Ainv

def Jacobian_col(phi_obs, scan_wls, etalon_wls, phi_real, g, Da, Et, Ncont):
    
    Jac = np.zeros((len(scan_wls), 2))
    
    Jac[:, 0] = df_g_col(phi_obs, scan_wls, etalon_wls, g, Da, phi_real, Et, Ncont)
    Jac[:, 1] = df_a_col(phi_obs, scan_wls, etalon_wls, g, Da, phi_real, Et, Ncont)
    
    return Jac

def Newton(phi_obs, scan_wls, etalon_wls, phi_real, g, Da, Et, Ncont, iters, 
           full = False):

    params = [g, Da]
    
    if full:
        full_params = np.zeros((iters + 1, 2))
        full_params[0] = params
        
    for i in range(iters):
    
        Jac = Jacobian_col(phi_obs, scan_wls, etalon_wls, phi_real, params[0], params[1], Et, Ncont)
        fun = fun_col(phi_obs, scan_wls, etalon_wls, params[0], params[1], phi_real, Et, Ncont)
    
        H = np.matmul(np.transpose(Jac), Jac)
        HI = svdsolve(H)
    
        Jf = np.matmul(np.transpose(Jac), fun)
    
        Change = np.matmul(HI, Jf)
    
        params -= Change
        
        if full:
            full_params[i + 1] = params

        print('\nJacobian :\n ', Jac)
        print('\nFunction :\n', fun)
        print('\nJac * Function\n', Jf)
        print('\nHI\n', HI)
        print('\nChange\n', Change)
        
    if full:
        return full_params
    else:
        return params
















