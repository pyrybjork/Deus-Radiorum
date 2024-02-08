#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: sergeykoldobskiy
"""
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

def num_integral(F,x):
    '''use if len(F) = len(x)'''
    F1 = np.array(F[:-1])
    F2 = np.array(F[1:])
    x1 = np.array(x[:-1])
    x2 = np.array(x[1:])
    INT=(F2*x2-F1*x1)/((np.log(F2/F1)/np.log(x2/x1)+1))
    return INT

def rigidity(Tn,p):
    'rigidity (GV) from kinetic energy per nucleon (GeV)'
    return p.A/p.Z * np.sqrt(Tn * (Tn + 2 * p.M/p.A)) 

class Proton:
    n = 'Proton'
    M = 0.938272
    Z = 1
    A = 1
    
class Alpha:
    n = 'Alpha'
    M =  3.7273
    Z = 2
    A = 4


def Jmod(Tn,phi,p,LIS):
    '''
    LIS modulated using FF approach [(m^2 s sr GeV)^-1]
    phi should be in GV
    '''
    if LIS=='Bos20':
        Jlis=Jlis_Bos20   
    Phi = phi * p.Z / p.A
    return Jlis(Tn+Phi,p) * Tn * (Tn+ 2*p.M/p.A)/(Tn+Phi)/(Tn+Phi+2*p.M/p.A)
    
def Jlis_Bos20(Tn,pa):
    '''
    originally in (m^2 s sr GV)^-1
    '''
    def L (x):
        return 1/(1+np.exp(-x))
    def G (x): 
        return np.exp(-np.power(x,2))
    
    if pa.Z == 1 or pa.n =='PAlpha':
        pa = Proton
        R = rigidity(Tn,pa)
        Rtilda = np.log(R)
        a = 3.2181e+2
        b = 1.1729e+4 
        c = 2.9232e+3
        d = 1.0621e+4
        f = 1.3867e+3
        g = 1.0399e+4
        h = 6.7166e-1
        i = 1.0528e+4
        l = 2.8240e+3
        m = 1.7430e-3
        n = 1.4742e+4
        o = 2.6617e+3
        p = 5.2830e-2
        q = 1.7160e+2
        r = 1.5000e-1
        s = 1.9222e+4
        t = 9.4790e-1
        u = 8.4900e-1
        J = np.where(R <= 2.5,
            a*R*R + b*(L(R)*L(R))+(c+d*Rtilda)*G(R/L(R)) - f - g*Rtilda - h * L(R) * np.power(i, G(R/L(R)) ),
            np.power(R,-2.7)*(- l - m*R + n*L(R) + o*G(p*R) + (q*Rtilda - np.power(s,(L(R)))*r)*np.cos(t+u*Rtilda)) 
                )
        J = J  * (Tn*pa.A+pa.M) / rigidity(Tn,pa) 
    elif pa.n == 'Alpha':
        R = rigidity(Tn,Alpha)
        Rtilda = np.log(R)
        a = 2.5869E+02
        b = 8.8956E+01
        c = 3.6399E+02
        d = 1.6337E+03
        f = 1.9527E+00
        g = 1.0469E+02
        h = 3.1229E+02
        i = 4.4168E+02
        l = 1.6856E+02
        m = 4.4660E+03
        n = 5.9618E+03
        o = 3.1158E+09
        p = 2.0719E+08
        q = 8.5045E+04
        r = 2.9055E+00
        s = 2.7152E+09
        t = 6.5850E+08
        u = 8.5045E+04 
        v = 1.6836E+02
        z = 5.5377E-05
        J = np.where(R <= 2.5,
            a - b * R  + c * np.sqrt(np.sin(R))+(d*R*G((f*R)**2) - g * R)*np.sin(R)-(h+i*R)*G(f*R),
            np.power(R,-2.0)*(l + m/R - n/R**2 + o/(p+q*R) + r*np.tanh(s/(t+u*R)) - v*z**(3/R))
                )
        
    J = J  * (Tn*pa.A+pa.M) / rigidity(Tn,pa) 
    return J 


hpa_to_gcm2 = 0.980665 
phi_9Sep2021_V2023 = 425.9 /1000
phi_Aug2021_V2023 = 432.6 / 1000
energies = np.logspace(-3,3,101)
J_Proton = Jmod(energies,phi_Aug2021_V2023,Proton,'Bos20')
HEPD_Proton_CR_2247 = np.array( [[0.05, 0.0587, 0.05435, 1037.17, 236.235],
                        [0.0587, 0.069, 0.06385, 1246.27, 276.624],
                        [0.069, 0.081, 0.075, 1363.73, 199.333],
                        [0.081, 0.0952, 0.0881, 1607.23, 324.497],
                        [0.0952, 0.112, 0.1036, 1726.01, 261.938],
                        [0.112, 0.131, 0.1215, 2022.32, 248.735],
                        [0.131, 0.154, 0.1425, 2208.48, 309.394],
                        [0.154, 0.181, 0.1675, 2427.18, 314.978],
                        [0.181, 0.213, 0.197, 2733.85, 436.517],
                        [0.213, 0.25, 0.2315, 2902.58, 688.096]] )

YF_ioniz_P_U6_10 = pd.read_csv('I_proton_U6_U10.txt',sep='\t',index_col='depth')
YF_ioniz_A_U6_10 = pd.read_csv('I_alpha_U6_U10.txt',sep='\t',index_col='depth')
    
def YF_interpol(depth,YF_file,E=None):
    '''
    Returns: energy, YF

    '''
    if depth > YF_file.index.max():
        YF_depth = YF_file.loc[YF_file.index.max()].values
    elif depth < YF_file.index.min():
        YF_depth = YF_file.loc[YF_file.index.min()].values
    else:
        min_ind = np.argpartition(np.abs(YF_file.index-depth),1)[:2]
        YF_1 = YF_file.iloc[min_ind[0]]
        YF_2 = YF_file.iloc[min_ind[1]]
        depth_1 = YF_file.index[min_ind[0]]
        depth_2 = YF_file.index[min_ind[1]]
        # breakpoint()
        diff_1 = (np.abs(depth_1 - depth))
        diff_2 = (np.abs(depth_2 - depth))
        if ~np.isfinite(diff_1):
            diff_1 = 0
        if ~np.isfinite(diff_2):
            diff_2 = 0
            
        diff_sum = np.abs(diff_1)+np.abs(diff_2)
        YF_depth = 10** (np.log10(YF_2) * np.abs(diff_1)/diff_sum + np.log10(YF_1) * np.abs(diff_2)/diff_sum )
    if E is not None:
        E_inteprol = interp1d(10+np.log10(YF_file.columns.astype('float')),np.log10(YF_depth),bounds_error=False,fill_value=np.nan)
        YF_depth_E = np.power(10,E_inteprol(10+np.log10(E)))
        return E, YF_depth_E
    else:
        return YF_file.columns.values.astype('float'), YF_depth.values

class C:
    def __init__(self, shape=(0,), dtype=float):
        """First item of shape is ingnored, the rest defines the shape"""
        self.shape = shape
        self.data = np.zeros((100,*shape[1:]),dtype=dtype)
        self.capacity = 100
        self.size = 0

    def update(self, x):
        if self.size == self.capacity:
            self.capacity *= 4
            newdata = np.zeros((self.capacity,*self.data.shape[1:]))
            newdata[:self.size] = self.data
            self.data = newdata

        self.data[self.size] = x
        self.size += 1

    def finalize(self):
        return self.data[:self.size]    



def calc_model_ionization_dose(depth,YF_file_p,YF_file_a,LIS_model='LIS_VP15', yf_units='depth',output='ioniz'):
    from _common.F_phys import part_heavy_AMS
    E = np.logspace(-2,3,51)
    E_p = E
    E_a = E
    
    if LIS_model == 'LIS_Bos20':
        J_p = Jmod(E_p,phi_9Sep2021_V2023,Proton,'Bos20') * 1e-4 #m2 to cm2
        J_a = Jmod(E_a,phi_9Sep2021_V2023,Alpha,'Bos20') * 1e-4 #m2 to cm2

    heavy_AMS = part_heavy_AMS(rigidity(E_a, Alpha))
    J_heavy = J_a * heavy_AMS #*4 is included in heavy_AMS

    ioniz_calc = C()
    
    
    E_p, YF_p = YF_interpol(depth,YF_file_p,E_p)
    E_a, YF_a = YF_interpol(depth,YF_file_a,E_a)
    P_c = 0.6 #GV
    ioniz_calc.update(np.nansum(num_integral(J_p*YF_p,E_p)[rigidity(E_p[1:],Proton)>P_c]) +\
                      np.nansum(num_integral(J_heavy*YF_a,E_a)[rigidity(E_a[1:],Alpha)>P_c]) \
                      ) 
    result = ioniz_calc.finalize()
    #here ionization is obtained in (ion par) / g / s
    #now let's transfer it to J/kg/s
    # 1 ion part is 35 eV
    # 1 eV = 1.6e-19 J
    print ('Atmospheric depth is '+str(depth)+' g/cm^2')
    print ('Ionization in units of [J/kg/s]:  '+ str(result * 35 * 1.6e-19 / 1e-3) )
    print ('Ionization in units of [Ion pairs / g / s]:  '+str(result))
    if output =='ioniz':
        result = result * 35 * 1.6e-19 / 1e-3    
        return result 
    elif output =='dose':
        return result
    
LIS_model = 'LIS_Bos20'

depth = 1000 #g/cm2
model_U6_10 = calc_model_ionization_dose(depth,YF_ioniz_P_U6_10,YF_ioniz_A_U6_10,LIS_model)


