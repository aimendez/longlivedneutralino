######################################################################
'''
LLP Flux Integral Ec (2.4) in [arXiv:1910.12839v2]
@Authors:
    - Pablo Candia (pcandia@uc.cl)
    - Giovanna Cottin ( @uai.cl)
    - Andres Mendez (aimendez@uc.cl)
    - Victor Munoz (vmmunoz2@uc.cl)
'''
######################################################################

import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline
from scipy import integrate
import os 
import time
import multiprocessing as mp
import sys 
from tqdm import tqdm
from config import BM 
cwd = os.getcwd()

#-------------------- BENCHMARK -------------------------#

BENCHMARK = sys.argv[1] #'BM1'
MESON = BM[BENCHMARK]['MESON'][0]
MODEL = BM[BENCHMARK]['MODEL']

#--------------------- CONSTANTS -------------------------#
# Earth Radius
Rt=6371

# Conversion factor GeV<->cm-1
CM_to_GeVM1=5.063e+13
GeVM1_to_cm=(5.06e13)**-1

# Masses 
mDp = 1.869 #GeV
mKp = 0.4937 #GeV
me = 0.511e-3 #GeV 

# decay length (ctau) in rest frame
ctaurest_Dp = (1.04e-12)*(1.52e+24)
ctaurest_Kp = (1.24e-08)*(1.52e+24)
ctaurest_Pip= (2.6e-08)*(1.52e+24) #GeV-1

#------------------ MESON DEPENDENT COLS ------------------#
if MESON == 'D':
    m_meson = mDp
    ctaurest = ctaurest_Dp
    col_names = ['H', 'L','C','E','Dp','Dm','Dsp','Dsm','R','X']
    col_flux1 = 'Dp'
    col_flux2 = 'Dm'
    E0 = 3

elif MESON == 'K':
    m_meson = mKp 
    ctaurest = ctaurest_Kp
    col_names = ['H', 'L', 'C', 'E', 'Pip', 'Pim', 'Kp', 'Km', 'R', 'X'] 
    col_flux1 = 'Kp'
    col_flux2 = 'Km'
    E0 = 0

#-------------------- READ DATA --------------------------#
df = pd.read_csv(cwd+f"/datafiles/MesonFlux/{MESON}_flux_{MODEL}.txt", header=None, delim_whitespace=True, names = col_names)


#------------------ USEFUL FORMULAS FROM  [arXiv:1910.12839v2] (victor paper) EC (2.4) - ----------------------#
def sqrMom_P(E, m): # GeV
    ''' Relativistic three-momentum p_i formula
        p_i^2 = E^2 - m^2
    '''
    if E<m:
        raise ValueError('Error ')
    else: 
         return np.sqrt( (E**2 - m**2) )
        
def beta(E,m):
    if E < m:
        raise ValueError('Error ')
    else:
        return np.sqrt(E**2-m**2)/E
    
def gamma(E,m):
    '''Relativistic factor gamma
    '''
    return E/m

def dec_length_Dp(E): # GeV-1
    '''First term in integral Ec(2.3)
    '''
    return gamma(E,m_meson)*beta(E,m_meson)*ctaurest
    
def lambdaFunc(mp ,md1, md2): # no dim
    '''Aux lambda function Ec(2.6)
    '''
    return  np.sqrt(1 + (md1**4)/(mp**4) + (md2**4)/(mp**4) - (2*md1**2)/(mp**2) - (2*md2**2)/(mp**2) - (2*md2**2*md1**2)/(mp**4) ) 

def twobody_dist(Ep,mp,md1,md2): # GeV M1
    '''Fraction of decaying parents that produce an A Ec(2.5)
        dn(E,Ep)/dE
    '''
    return 1/(sqrMom_P(Ep,mp)*lambdaFunc(mp,md1,md2)) 

#----------------------- FLUX INTEGRAL -----------------------------#
def flux_E( cos_i, x_i ):
    ''' Meson Flux as a function of the energy FLUX=FLUX(E) for a given set of values (cos(theta)_i, X_i).
        Second term in integral (2.3)
        dPHI_p(Ep)/dEp
    '''
    df_tmp = df[(df.C==cos_i) & (df.X==x_i)]
    return interpolate.InterpolatedUnivariateSpline(x=df_tmp.E, y=df_tmp[col_flux1], k=1), interpolate.InterpolatedUnivariateSpline(x=df_tmp.E, y=df_tmp[col_flux2], k=1)

    
def lim_integral(mM,mN,ml,EN):
    '''Integal limits given by Ec(2.13) 
    '''
    lim_integral = [( 1/(2 *mN**2) )*(EN *(-ml**2+mM**2+mN**2)-np.sqrt( (EN-mN)*(EN+mN)*(-ml-mM+mN)*(ml-mM+mN)*(-ml+mM+mN)*(ml+mM+mN)) ),  
                    ( 1/(2 *mN**2) )*(EN *(-ml**2+mM**2+mN**2)+np.sqrt( (EN-mN)*(EN+mN)*(-ml-mM+mN)*(ml-mM+mN)*(-ml+mM+mN)*(ml+mM+mN)) ) ]
    return lim_integral

def dL_Dp(E): # cm
    '''Conversion factor
    ''' 
    return dec_length_Dp(E)*GeVM1_to_cm


def Int_LogEp(LogEp, ma, x, ic, ix, rho):
    '''Function to integrate in Ec(2.4)  
    '''
    flux_D = flux_array[ ix + 50*ic ]
    return   np.exp(LogEp) * (flux_D[0](np.exp(LogEp)) + flux_D[1](np.exp(LogEp))) *twobody_dist(np.exp(LogEp),m_meson,me,ma)* (rho*dL_Dp(np.exp(LogEp)))**-1


#----------------- INTEGRAL COMPUTATION ---------------------------#
def compute_I(MA, cos):
    ''' Computation of integral (2.3) for a given mass MA and direction cos(theta).
    '''
    ic = int(list(df.C.unique()).index(cos))
    df_tmp = pd.DataFrame()
    X_ = sorted( df[df.C == cos].X.unique() )
    for ix, X in enumerate(tqdm(X_)):
        rho = df[(df.C==cos) & (df.X==X)].R.values[0]
        L = df[(df.C==cos) & (df.X==X)].L.values[0]
        H = df[(df.C==cos) & (df.X==X)].H.values[0]
        for E in np.array(df['E'].unique()[E0:25]):
            lim = lim_integral(m_meson,me,MA,E)
            I = integrate.quad(Int_LogEp, np.log(lim[0]), np.log(lim[1]), args=(MA,X, ic, ix, rho), epsabs= 0.00, epsrel= 1e-2, limit=300 )
            df_tmp = df_tmp.append(pd.DataFrame({'E':[E],'X': [X], 'L':[L], 'H':[H], 'cosine':[cos], 'm':[MA],'rho':[rho],'flux':[I[0]]}))
    return df_tmp 

flux_array = np.array([flux_E(cos, x) for cos in df.C.unique() for x in sorted(df[df.C== cos].X.unique())  ])


#=====================================================================================================================================#
#=====================================================================================================================================#
#=====================================================================================================================================#


if __name__ == '__main__':
    
    start = time.time()
    #-------------- BENCHMARK ---------------------#
    ma_vec = BM[BENCHMARK]['MA']
    args = [ (ma,cos) for ma in ma_vec for cos in df.C.unique()]

    #------------- MULTIPROCESSING ----------------#
    p = mp.Pool(processes = mp.cpu_count()-1)
    print('BENCHMARK SELECTED: ', BENCHMARK)
    print('PPROCESSES: ', mp.cpu_count()-1)
    print('CALCULATING INTEGRAL FLUX ...')
    result = pd.concat(p.starmap(compute_I, args))
    p.close()
    p.join()

    #----------- SAVE RESULTS AS CSV --------------#
    result = result.sort_values(by=['E','m', 'cosine', 'X']).reset_index(drop=True)
    result = result[['flux', 'E', 'X', 'cosine', 'm', 'L', 'rho', 'H']] 
    result.to_csv(cwd+f'/datafiles/NeutralinoFlux/NeutralinoFlux_{BENCHMARK}_{MODEL}.csv')
    print(f'/datafiles/NeutralinoFlux/NeutralinoFlux_{BENCHMARK}_{MODEL}.csv \t SAVED')
   
   #----------------------------------------------#
    end = time.time()
    print("COMPLETE")
    print('TOTAL TIME: ' + str(round(end-start,4)) + '(s)')
    print()
    print('########################################################################################')
    print()


#=====================================================================================================================================#
#=====================================================================================================================================#
#=====================================================================================================================================#

