# ============================ IMPORTS ====================================== #
 
import os
import numpy as np
import matplotlib.pyplot as plt
 
# ============================= CONFIG ====================================== #
 
plt.style.use('dark_background')
 
# ============================= AUXILIAR ==================================== #

def sec(x):
    return 1 / np.cos(x)

def alpha1(a,b,F):
    return 2*a*b**2*np.sqrt(F)

def der_alpha1(b,F):
    #Derivative respect to 'a'
    return 2 * b ** 2 * np.sqrt(F)

def alpha2(a,b,F):
    alpha1=2*a*b**2*np.sqrt(F)
    return 2*alpha1*np.sqrt(F+1)

def der_alpha2(a,b,F):
    #Derivative respect to 'a'
    return 4*b**2*np.sqrt(F*(F+1))

def gamma1(a, F):
    return np.sqrt(F) * np.sin(a)

def der_gamma1(a,F):
    return np.sqrt(F)*np.cos(a)

def gamma2(a,b,F):
    return np.sqrt(F) * np.sin(a * (1- b ** 2))

def der_gamma2(a,b,F):
    #Derivative respect to 'a'
    return np.sqrt(F) * (1-b) * np.cos(a * (1- b ** 2))

def gamma3(F):
    return np.sqrt(F / (F+1))

def gamma4(a,b,F):
    return np.tan(a / 2 * (1-b ** 2)) / np.sqrt(F + 1)
 
def der_gamma4(a,b,F):
    #Derivative respect to 'a'
    return (1 - b) / (2 * np.sqrt(F+1)) *sec(a / 2 * (1-b ** 2)) ** 2

def gamma5(a,F):
    return np.tan(a/2)/np.sqrt(F+1)

def der_gamma5(a,F):
    return 1 / (2 * np.sqrt(F+1)) * sec(a/2) **2


# ============================= ETALON FS =================================== #


def RealE(alpha1, gamma1, gamma2):
    
    Ere = (2 / alpha1) * (np.arctan(gamma1) - np.arctan(gamma2))
    
    return Ere

def imE(R, alpha2,gamma3,gamma4,gamma5):

    Eim=(2 / alpha2) * (1 + R) / (1 - R) * \
    (np.log(((1+gamma3) ** 2 + gamma4 ** 2) / ((1-gamma3) ** 2 + gamma4 ** 2)) - \
     np.log(((1+gamma3) ** 2 + gamma5 ** 2) / ((1-gamma3) ** 2 + gamma5 ** 2)))
    
    return Eim
    
def der_RealE(a, b, F, alpha1, gamma1, gamma2):

    A = 1 / a * (np.arctan(gamma2) - np.arctan(gamma1))
    
    B = der_gamma1(a, F) / (1 + gamma1 ** 2)
    
    C = der_gamma2(a, b, F) / (1 + gamma2 ** 2)
    
    return 2 / (alpha1) * (A + B - C)    
    
def der_ImE(R,a,b,F,alpha2,gamma3,gamma4,gamma5):
 
    #Derivatives of alpha and gammas
    alpha2prima = der_alpha2(a,b,F)
    gamma4prima = der_gamma4(a,b,F)
    gamma5prima = der_gamma5(a,F)
    
    #Denominators of derivative of ratio
    denom1=((1-gamma3)**2+gamma4**2)*((1+gamma3)**2+gamma4**2)
    denom2=((1-gamma3)**2+gamma5**2)*((1+gamma3)**2+gamma5**2)

    der1=-alpha2prima/alpha2**2\
    *np.log(((1+gamma3)**2+gamma4**2)/((1-gamma3)**2+gamma4**2))\
    -(1/alpha2)*8*gamma3*gamma4*gamma4prima/denom1

    der2=-alpha2prima/alpha2**2\
    *np.log(((1+gamma3)**2+gamma5**2)/((1-gamma3)**2+gamma5**2))\
    -(1/alpha2)*8*gamma3*gamma5*gamma5prima/denom2
    
    return 2 * (1 + R) / (1-R) *(der1-der2)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

