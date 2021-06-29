######################################################################
'''
Number of Events dN_ev/dEdcos in [arXiv:1910.12839v2]
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
import warnings
warnings.filterwarnings("ignore")
cwd = os.getcwd()


#-------------------- BENCHMARK -------------------------#
BENCHMARK = sys.argv[1] #'BM1'
MODEL = BM[BENCHMARK]['MODEL']

#--------------------------------------------------------#


def GetFlux(BENCHMARK):
	df=pd.read_csv( cwd + f'/datafiles/NeutralinoFlux/NeutralinoFlux_{BENCHMARK}_{MODEL}.csv', index_col = 0).reset_index(drop=True)
	return df

def GetPheno(BENCHMARK):
	df_ct = pd.read_csv(cwd + f'/datafiles/NeutralinoPheno/NeutralinoPheno_{BENCHMARK}.csv' )
	return df_ct

def GetAeff(detector = 'SK'):
	df_surface = pd.read_csv( cwd + f"/datafiles/Aeff/{detector}_Aeff.txt", header=None, delim_whitespace=True, names = ['A_top','A_side','cos_theta','ctau','A_total'])
	df_surface = df_surface.sort_values(by=['cos_theta', 'ctau'], ascending=True)
	x = np.array(df_surface['ctau'].unique(), dtype='float64')
	y = np.array(df_surface['cos_theta'].unique(), dtype='float64')
	z = np.array(df_surface['A_total'].values, dtype='float64')
	Z = z.reshape(len(y), len(x))
	return RectBivariateSpline(y, x, Z, kx=1, ky=1)

def Chi_sq_bin(s,b,d):
	return s + b -d + d*np.log(d/(s+b))

def ctau_lab_func(ct_rest,m,E):
    return ct_rest*(E/m)


def ComputeNEvents(BENCHMARK='BM1', detector = 'SK'):
	df = GetFlux(BENCHMARK)
	cos_values = df["cosine"].unique()
	df["m_round"] = df["m"].round(2)
	mass_values = df["m_round"].unique()
	E_vec=df["E"].unique()

	df_ct = GetPheno(BENCHMARK)
	ctau_values = df_ct['ctau']
	Br_values = df_ct['Br']
	df_ct['m_round'] = df_ct["m"].round(2)
	df_ct = df_ct.sort_values(by=['m', BM[BENCHMARK]['LAM_PROD'], BM[BENCHMARK]['LAM_DEC']])

	Aeff = GetAeff(detector)
	eficiencia = 0.75
	fn = eficiencia*2*np.pi*(5326)*24*60*60
	 
	list_tmp = []
	evts_per_cos_g = []
	Evts_g = []
	E_multiGeV=E_vec[0:-9] 
	print('COMPUTING EVENT INTEGRAL ...\n')
	i = 0
	for idx, row in tqdm(df_ct.iterrows(), total=df_ct.shape[0], position=0, leave=True):
		ct_r = row['ctau']
		m_val = row['m_round']
		lambda_prod = row[BM[BENCHMARK]['LAM_PROD']]
		lambda_dec =  row[BM[BENCHMARK]['LAM_DEC']]
		Br_val = row['Br']
		for icos, cos in enumerate(cos_values):
			ALP_surface=[]
			for ie,e in enumerate(E_vec):
				df_tmp = df[ (df.cosine == cos) & (df.m_round == m_val) & (df.E == e)]
				PHI_val= df_tmp.flux.values
				X_val= df_tmp.X.values
				dflux_dXdE = interpolate.interp1d(X_val,PHI_val)
				L_func = interpolate.interp1d(X_val, df[ (df.cosine == cos) ].L.unique())
				ALP_surf = lambda LogX, Ea, ma, ct: np.exp(-(L_func(np.exp(LogX))*ma)*(Ea*ct)**-1)*np.exp(LogX)*float(dflux_dXdE(np.exp(LogX)))
				I=integrate.quad(ALP_surf,np.log(min(X_val)),np.log(max(X_val)),args=(e,m_val,ct_r), epsabs= 1e-10, epsrel= 1e-10,limit=300)
				ALP_surface.append(I[0])

			ALP_flux  = interpolate.interp1d(E_vec,ALP_surface)
			def Integrando_Ea(LogEa,ma,ctaur,cosTh):
				ct_val=ctau_lab_func(ctaur,ma,np.exp(LogEa))
				return np.exp(LogEa)*ALP_flux(np.exp(LogEa))*Aeff(cosTh,ct_val)
        
			I = integrate.quad((Integrando_Ea),np.log(min(E_multiGeV)),np.log(max(E_multiGeV)), args=(m_val,ct_r,cos), epsabs= 1e-10, epsrel= 1e-10,limit=300)
			evts_per_cos_g.append(fn*I[0]*Br_val)
			nevs = fn*I[0]*Br_val
			dict_tmp = { 
						 'nev': nevs ,
			             'E':e,
			             'ctau':ct_r,
			             'm':m_val,
			              BM[BENCHMARK]['LAM_PROD']:lambda_prod,
						  BM[BENCHMARK]['LAM_DEC']:lambda_dec,
			             'Br': Br_val,
			             'cos':cos
			           }
			list_tmp.append(pd.DataFrame(dict_tmp, index=[idx]))
			i+=1
		Evts_g.append(evts_per_cos_g)
	df = pd.concat(list_tmp)
	return df




#=====================================================================================================================================#
#=====================================================================================================================================#
#=====================================================================================================================================#

if __name__ == '__main__':
	start = time.time()
	#---------------------------------------------------------------------#
	detector = 'SK'
	df = ComputeNEvents(BENCHMARK, detector)

	#---------------------------------------------------------------------#
	df.to_csv( cwd + f'/datafiles/EventDistribution/NEVENTS_{BENCHMARK}_{MODEL}.csv')
	print(cwd + f'/datafiles/EventDistribution/NEVENTS_{BENCHMARK}_{MODEL}.csv - SAVED')

	#---------------------------------------------------------------------#
	end = time.time()
	print("COMPLETE")
	print('TOTAL TIME: ' + str(round(end-start,4)) + '(s)')
	print()
	print('########################################################################################')
	print()

	
#=====================================================================================================================================#
#=====================================================================================================================================#
#=====================================================================================================================================#
