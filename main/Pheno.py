######################################################################
'''
RpV Pheno. BR and ctau in [arXiv:1910.12839v2]
@Authors
    - Pablo Candia (pcandia@uc.cl)
    - Giovanna Cottin ( @uai.cl)
    - Andres Mendez (aimendez@uc.cl)
    - Victor Munoz (vmmunoz2@uc.cl)
'''
######################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import math
import pandas as pd
import sys 
from math import sqrt
from math import pi as PI
from config import BM 
import time


#-------------------- BENCHMARK -------------------------#
BENCHMARK = sys.argv[1] #'BM1'

#--------------------- CONSTANTS -------------------------#
# Fundamental physical constants
SW2 = 0.231 #sin^2 theta_w
TW = sqrt(SW2)/sqrt(1 - SW2) # tan theta_w
g1 = 0.356 # U(1) gauge coupling constant
g2 = 0.650 # SU(2) gauge coupling constant
hbar = 1.05*10**-34
JtoGeV = 6.24*10**9 # 1 J = 6.24*10^9 GeV
Inv_GeV_to_m = 1.97e-16 # 1/GeV = 1.97e-16 m

# Trilinear bino couplings (eqs. (24) to (28) in [1511.07436])

gslep = g2*TW/sqrt(2)
gsneut = g2*TW/sqrt(2)
gsuL = -g2*TW/(3*sqrt(2))
gsdL = 5*g2*TW/(3*sqrt(2))
gsdR = -2*g2*TW/(3*sqrt(2))

# Masses in GeV

## Mesons

mdPlus = 1.87
mdstarPlus = 2.01
mdsPlus = 1.97
mdsstarPlus = 2.11
mkstarPlus = 0.892
mkPlus = 0.494
mk0l = 0.498
mk0s = 0.498
mkstar = 0.896
mpi0 = 0.135
mpiPlus = 0.140
meta = 0.548
metap = 0.957
mrho = 0.776
mrhoPlus = 0.775
momega = 0.783
mphi = 1.02

## Leptons

mtau = 1.78
mmu = 0.106
me = 0.000511

## Quarks

mup = 2.3e-3
mdown = 4.8e-3
mcharm = 1.28
mstrange = 95e-3

# Meson decay constants (from [1511.07436] and [0612278])

FdPlus = 0.205
FdstarPlus = FdPlus
FdsPlus = 0.259
FdsstarPlus = FdsPlus 
Fkstar = 0.230
FkstarPlus = 0.230 # Same as Fkstar?
FkPlus = 0.156
Fk0l = FkPlus/sqrt(2)
Fk0s = FkPlus/sqrt(2)
Fpi0 = 0.130
FpiPlus = Fpi0
Frho = 0.22
FrhoPlus = Frho # I couldn't find this one
Fomega = Frho # I couldn't find this one
Fphi = 0.23
Feta = -0.142
Fetap = 0.038

# Lifetimes (lt) in seconds

dPluslt = 1.04*10**-12
kstarPluslt = 7.35*10**-20
pi0lt = 8.52*10**-17
piPluslt = 2.6*10**-8
kPluslt = 1.24*10**-8
k0slt = 0.895*10**-10
k0llt = 5.12*10**-8
dsPluslt = 504*10**-15

# Widths in GeV

dPluswidth = JtoGeV*hbar/dPluslt
kstarPluswidth = JtoGeV*hbar/kstarPluslt
dstarPluswidth = 83.4*10**-6
rhowidth = 0.149
rhoPluswidth = rhowidth
pi0width = JtoGeV*hbar/pi0lt
piPluswidth = JtoGeV*hbar/piPluslt
omegawidth = 8.49*10**-3
phiwidth = 4.26*10**-3
kPluswidth = JtoGeV*hbar/kPluslt
k0swidth = JtoGeV*hbar/k0slt
k0lwidth = JtoGeV*hbar/k0llt
kstarwidth = 47.4*10**-3
etawidth = 1.31*10**-6
etapwidth = 0.188*10**-3
dsPluswidth = JtoGeV*hbar/dsPluslt
dsstarPluswidth = 1.9*10**-3

# Mass, decay constant and width dictionaries (Note that for pseudoscalar mesons we use F := f_S_M as defined in eq. (35) of [1511.07436])

mesons = {"K0L": {'F': Fk0l*mk0l**2/(mstrange + mdown), 'm': mk0l, 'width': k0lwidth},
          "K0S": {'F': Fk0s*mk0s**2/(mstrange + mdown), 'm': mk0s, 'width': k0swidth},
          "rho": {'F': Frho, 'm': mrho, 'width': rhowidth},
          "rho+": {'F': FrhoPlus, 'm': mrhoPlus, 'width': rhoPluswidth},
          "omega": {'F': Fomega, 'm': momega, 'width': omegawidth},
          "phi": {'F': Fphi, 'm': mphi, 'width': phiwidth},
          "eta": {'F': Feta*meta**2/(mstrange + mstrange), 'm': meta, 'width': etawidth},
          "etap": {'F': Fetap*meta**2/(mstrange + mstrange), 'm': metap, 'width': etapwidth},
          "Pi0": {'F': Fpi0*mpi0**2/(mdown + mdown), 'm': mpi0, 'width': pi0width},
          "Pi+": {'F': FpiPlus*mpiPlus**2/(mup + mdown), 'm': mpiPlus, 'width': piPluswidth},
          "K*": {'F': Fkstar, 'm': mkstar, 'width': kstarwidth},
          "D+": {'F': FdPlus*mdPlus**2/(mcharm + mdown), 'm': mdPlus, 'width': dPluswidth},
          "D*+": {'F': FdstarPlus, 'm': mdstarPlus, 'width': dstarPluswidth},
          "Ds+": {'F': FdsPlus*mdsPlus**2/(mcharm + mstrange), 'm': mdsPlus, 'width': dsPluswidth},
          "Ds*+": {'F': FdsstarPlus, 'm': mdsstarPlus, 'width': dsstarPluswidth},
          "K+": {'F': FkPlus*mkPlus**2/(mup + mstrange), 'm': mkPlus,'width': kPluswidth},
          "K*+": {'F': FkstarPlus, 'm': mkstarPlus, 'width': kstarPluswidth}}

# Dictionary with RpV parameters and the neutralino decay final states associated to them. It is understood that the neutral mesons are produced together with a neutrino.

lambda_dict = {"lam111": {'charged': ["Pi+","rho+"], 'neutral': ["Pi0","eta","etap","rho","omega"], 'lepton': 'e'},
               "lam112": {'charged': ["K+","K*+"], 'neutral': ["K0L","K0S","K*"], 'lepton': 'e'},
               "lam121": {'charged': ["D+","D*+"], 'neutral': ["K0L","K0S","K*"], 'lepton': 'e'},
               "lam122": {'charged': ["Ds+","Ds*+"], 'neutral': ["eta","etap","phi"], 'lepton': 'e'},
               "lam211": {'charged': ["Pi+","rho+"], 'neutral': ["Pi0","eta","etap","rho","omega"], 'lepton': 'mu'},
               "lam212": {'charged': ["K+","K*+"], 'neutral': ["K0L","K0S","K*"], 'lepton': 'mu'},
               "lam221": {'charged': ["D+","D*+"], 'neutral': ["K0L","K0S","K*"], 'lepton': 'mu'},
               "lam222": {'charged': ["Ds+","Ds*+"], 'neutral': ["eta","etap","phi"], 'lepton': 'mu'}}

pseudoscalar_meson_list = ['Pi0','Pi+','eta','etap','K+','K0L','K0S','D+','Ds+']

vector_meson_list = ['rho','rho+','omega','phi','K*','K*+','D*+','Ds*+']

lepton_mass = {'e': me,
               'mu': mmu,
               'tau': mtau,
               'nu': 0.0}

# Effective four-fermion couplings (Eq. (30) to (33)). We consider the same values for all soft sfermion masses and family degeneracy.

def G(lepton,meson,param):
    if lepton == 'nu' and meson in pseudoscalar_meson_list:
        return param*(gsneut - 0.5*gsdL - 0.5*gsdR) # G^S_\nu
    elif lepton in ['e','mu','tau'] and meson in pseudoscalar_meson_list:
        return param*(0.5*gsuL + 0.5*gsdR - gslep) # G^S_\ell
    elif lepton == 'nu' and meson in vector_meson_list:
        return param*(0.25*gsdL + 0.25*gsdR) # G^T_\nu
    elif lepton in ['e','mu','tau'] and meson in vector_meson_list:
        return param*(0.25*gsuL + 0.25*gsdR) # G^T_\ell
    else:
        raise ValueError('Wrong meson and/or lepton.')

# Lambda function

def LambdaSQ(x,y,z):
    if x**2 + y**2 + z**2 - 2*x*y - 2*x*z - 2*y*z > 0:
        return sqrt(x**2 + y**2 + z**2 - 2*x*y - 2*x*z - 2*y*z)
    else:
        return 0

# Branching ratio of a meson decaying to neutralino + lepton. It receives a parameter (lam_P) that corresponds to the production RpV parameter. This parameter has to be specified and passed as a string 'prod' in the function that calulates the neutralino decay width. It also receives a meson and lepton that MUST BE CONSISTENT WITH THE CHOICE OF 'lam_P'.

def BrMeson(mchi,meson,msferm,lepton,lam_P):
    M_mass = mesons[meson]['m']
    decay_const = mesons[meson]['F']
    M_width = mesons[meson]['width']
    ml = lepton_mass[lepton]
    param = lam_P/(msferm**2)
    if M_mass > (ml + mchi):
        if meson in pseudoscalar_meson_list:
            partial_width_pseudo = LambdaSQ(M_mass**2,mchi**2,ml**2) * G(lepton,meson,lam_P/msferm**2)**2 * decay_const**2 * (M_mass**2 - mchi**2 - ml**2)/(64*PI*M_mass**3)
            return partial_width_pseudo/(M_width + partial_width_pseudo)
        elif meson in vector_meson_list:
            partial_width_vector = LambdaSQ(M_mass**2,mchi**2,ml**2) * G(lepton,meson,lam_P/msferm**2)**2 * decay_const**2 * (M_mass**2*(M_mass**2 + mchi**2 + ml**2) - 2*(mchi**2 - ml**2)**2)/(3*PI*M_mass**3)
            return partial_width_vector/(M_width + partial_width_vector)
        else:
            raise ValueError('Wrong meson name.')
    else:
        return 0

# Neutralino partial decay width to a meson and a lepton. This function receives a lambda primer trilinear RpV parameter 'trili', the neutralino mass, the meson, the lepton, and the sfermion mass scale.

def Gamma_Chi(mchi,meson,msferm,lepton,trili):
    M_mass = mesons[meson]['m']
    decay_const = mesons[meson]['F']
    ml = lepton_mass[lepton]
    param = trili/(msferm**2)
    if mchi > (ml + M_mass):
        if meson in pseudoscalar_meson_list:
            return LambdaSQ(mchi**2,M_mass**2,ml**2) * G(lepton,meson,trili/msferm**2)**2 * decay_const**2 * (mchi**2 + ml**2 - M_mass**2)/(128*PI*mchi**3)
        elif meson in vector_meson_list:
            return LambdaSQ(mchi**2,M_mass**2,ml**2) * G(lepton,meson,trili/msferm**2)**2 * decay_const**2 * (2*(mchi**2 - ml**2)**2 - M_mass**2*(M_mass**2 + mchi**2 + ml**2))/(2*PI*mchi**3)
        else:
            raise ValueError('Wrong meson name.')
    else:
        return 0

# Neutralino total width. This function receives the neutralino mass, the sfermion mass scale, and the values of the production and decay RpV parameters. It also receives strings lam_P, lam_D that specify what production and decay parameters are activated. This choice also determines what decay channels are available and contribute to the neutralino decay width, as specified in lambda_dict. If one lambda is used for production and decay of the neutralino, then the argument lam_P has to coincide with lam_D.

def width_chi(mchi, msferm, prod, dec, lam_P, lam_D):
    Gamma_neut = 0
    for lam in lambda_dict:
        if prod == lam:
            for mes in lambda_dict[lam]['neutral']:
                Gamma_neut += 2*Gamma_Chi(mchi,mes,msferm,'nu',lam_P)
            for mes in lambda_dict[lam]['charged']:
                Gamma_neut += 2*Gamma_Chi(mchi,mes,msferm,lambda_dict[lam]["lepton"],lam_P)
        if prod != dec:
            if dec == lam:
                for mes in lambda_dict[lam]['neutral']:
                    Gamma_neut += 2*Gamma_Chi(mchi,mes,msferm,'nu',lam_D)
                for mes in lambda_dict[lam]['charged']:
                    Gamma_neut += 2*Gamma_Chi(mchi,mes,msferm,lambda_dict[lam]["lepton"],lam_D)
    return Gamma_neut

# Branching ratios of neutralino decays. This function takes a lepton and a meson and verifies whether these final states are allowed given a choice of the production and decay parameters. If they are, it returns a (generally) non-zero branching ratio.

def Br_chi(mchi, mes, msferm, lepton, prod, dec, lam_P, lam_D):
    Br = 0
    total_width = width_chi(mchi, msferm, prod, dec, lam_P, lam_D)
    if mes in lambda_dict[prod]['charged'] and lepton == lambda_dict[prod]['lepton']:
        Br += Gamma_Chi(mchi,mes,msferm,lepton,lam_P)/total_width
    if mes in lambda_dict[prod]['neutral'] and lepton == 'nu':
        Br += Gamma_Chi(mchi,mes,msferm,lepton,lam_P)/total_width
    if dec != prod:
        if mes in lambda_dict[dec]['charged'] and lepton == lambda_dict[dec]['lepton']:
            Br += Gamma_Chi(mchi,mes,msferm,lepton,lam_D)/total_width
        if mes in lambda_dict[dec]['neutral'] and lepton == 'nu':
            Br += Gamma_Chi(mchi,mes,msferm,lepton,lam_D)/total_width
    return Br
 
# Electron-like branching ratios of neutralino decays (Eq. (4.8) in 1910.12839)

def Br_elike(mchi, msferm, prod, dec, lam_P, lam_D):
    Br = 0    
    total_width = width_chi(mchi, msferm, prod, dec, lam_P, lam_D)
    for lam in ['lam111','lam112','lam121','lam122']:
        if prod == lam:
            for mes in lambda_dict[lam]['neutral']:
                Br += 2*Gamma_Chi(mchi,mes,msferm,'nu',lam_P)/total_width
            for mes in lambda_dict[lam]['charged']:
                Br += 2*Gamma_Chi(mchi,mes,msferm,lambda_dict[lam]["lepton"],lam_P)/total_width
        if prod != dec:
            if dec == lam:
                for mes in lambda_dict[lam]['neutral']:
                    Br += 2*Gamma_Chi(mchi,mes,msferm,'nu',lam_D)/total_width
                for mes in lambda_dict[lam]['charged']:
                    Br += 2*Gamma_Chi(mchi,mes,msferm,lambda_dict[lam]["lepton"],lam_D)/total_width
    for lam in ['lam211','lam212','lam221','lam222']:
        if prod == lam:
            for mes in lambda_dict[lam]['neutral']:
                Br += 2*Gamma_Chi(mchi,mes,msferm,'nu',lam_P)/total_width
        if prod != dec:
            if dec == lam:
                for mes in lambda_dict[lam]['neutral']:
                    Br += 2*Gamma_Chi(mchi,mes,msferm,'nu',lam_D)/total_width
    return Br

# Lifetime in m

def ctau_chi(mchi, msferm, prod, dec, lam_P, lam_D):
    return Inv_GeV_to_m / width_chi(mchi, msferm, prod, dec, lam_P, lam_D)

# Methods to generate list with masses, decay lengths and branching ratios.

# This method generates a list where the production and decay parameters are set to the same value. It generates a list where the mass and the value of the trilinear parameter change. It also works when the production and decay parameters coincide.

def GetPheno(mchi_vec, prod_meson, prod_lepton, lamP_vec, lamD_vec, msferm_t, prod, dec,equal_coupling):
    mt_to_km = 1e-3
    ct_list, BrP_list, BrD_list, Br_list, Br_charged_list, lamP_list, lamD_list, m_list = ([] for i in range(8))
    for mchi in mchi_vec:
        for lamP in lamP_vec:
            for lamD in lamD_vec:
                if (equal_coupling == True) and (lamD != lamP):
                    continue
                Br_charged = 0
                m_list.append(mchi)
                lamP_list.append(lamP)
                lamD_list.append(lamD)
                ct_list.append(mt_to_km*ctau_chi(mchi, msferm_t, prod, dec, lamP, lamD))
                BrD_list.append(Br_elike(mchi, msferm_t, prod, dec, lamP, lamD))
                BrP_list.append(BrMeson(mchi, prod_meson, msferm_t, prod_lepton, lamP))
                Br_list.append(BrMeson(mchi, prod_meson, msferm_t, prod_lepton, lamP)*Br_elike(mchi, msferm_t, prod, dec, lamP, lamD))
                for mes in lambda_dict[dec]['charged']:
                    Br_charged += 2*Br_chi(mchi, mes, msferm_t, lambda_dict[dec]['lepton'], prod, dec, lamP, lamD)
                for mes in lambda_dict[prod]['charged']:
                    Br_charged += 2*Br_chi(mchi, mes, msferm_t, lambda_dict[prod]['lepton'], prod, dec, lamP, lamD)
                Br_charged_list.append(Br_charged)
    return pd.DataFrame.from_dict({'ctau': ct_list, 'm': m_list, prod: lamP_list, dec: lamD_list, 'BrP': BrP_list, 'BrD': BrD_list, 'Br': Br_list, 'Br_charged': Br_charged_list})



#=====================================================================================================================================#
#=====================================================================================================================================#
#=====================================================================================================================================#

if __name__ == '__main__':

    start = time.time()
    #----------------------------------------------#

    meson = BM[BENCHMARK]['MESON']
    ma_vec = BM[BENCHMARK]['MA']
    lam_prod = BM[BENCHMARK]['LAM_PROD']
    lam_prod_vec = BM[BENCHMARK]['LAM_PROD_RANGE']
    lam_dec = BM[BENCHMARK]['LAM_DEC']
    lam_dec_vec = BM[BENCHMARK]['LAM_DEC_RANGE']
    msferm_t =  BM[BENCHMARK]['MSFERM']
    lepton = BM[BENCHMARK]['LEPTON']
    equal_coupling = BM[BENCHMARK]['EQUAL_COUPLING']

    df_ct_eq_lam = GetPheno(ma_vec, meson, lepton, lam_prod_vec, lam_dec_vec, msferm_t, lam_prod, lam_dec, equal_coupling)
    df_ct_eq_lam.to_csv(f'./datafiles/NeutralinoPheno/NeutralinoPheno_{BENCHMARK}.csv', index=None, header=True)
    print(f'./datafiles/NeutralinoPheno/NeutralinoPheno_{BENCHMARK}.csv \t SAVED')

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

