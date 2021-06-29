import numpy as np

me = 0.000511
mpi0 = 0.135
mkPlus = 0.494


BM  = {
		"BM1" : {	#BM1 D-MESONS lam121 = lam112 (fig 5a) #CHECKED
					"MESON": 'D+',
					"MA": np.linspace(0.5, 1.0, 40),
					"LAM_PROD": 'lam121',
					"LAM_PROD_RANGE": np.logspace(-1.0, 1.0, 40),
					"LAM_DEC": 'lam112',
					"LAM_DEC_RANGE": np.logspace(-1.0, 1.0, 40),
					"MSFERM": 1000,
					"LEPTON": 'e',
					"EQUAL_COUPLING":True,
					"MODEL": "SYBILL",
				},

		"BM1_AUX" : {	#BM1 AUX D-MESONS - Same as BM1 but bigger range for Ma (fig 2)
					"MESON": 'D+',
					"MA": np.linspace(0.2, 1.0, 20) ,
					"LAM_PROD": 'lam121',
					"LAM_PROD_RANGE": np.logspace(-3.0, 1.0, 20),
					"LAM_DEC": 'lam112',
					"LAM_DEC_RANGE": np.logspace(-3.0, 1.0, 20),
					"MSFERM": 1000,
					"LEPTON": 'e',
					"EQUAL_COUPLING":True,
					"MODEL": "SYBILL",
				},

		"BM2" : {	#BM2 D-MESONS lam121 (fig 5b) # CHECKED
					"MESON": 'D+',
					"MA": np.linspace(0.5, 1.0, 40) ,
					"LAM_PROD": 'lam121',
					"LAM_PROD_RANGE": np.logspace(-1.0, 1.0, 40), #np.logspace(-3,1,40),
					"LAM_DEC": 'lam121',
					"LAM_DEC_RANGE": np.logspace(-1.0, 1.0, 40), # np.logspace(-3,1,40),
					"MSFERM": 1000,
					"LEPTON": 'e',
					"EQUAL_COUPLING":True,
					"MODEL":"SYBILL",
				},

		

		"BM3": {    #BM3 KAONS (fig 2, fig6a)
					"MESON": 'K+',
					"MA": np.linspace(mpi0+me,  mkPlus - 2*me, 15) ,
					"LAM_PROD": 'lam112',
					"LAM_PROD_RANGE": np.logspace(-4.0, 1.0, 20),
					"LAM_DEC": 'lam111',
					"LAM_DEC_RANGE": np.logspace(-4.0, 1.0, 20),
					"MSFERM": 1000,
					"LEPTON": 'e',
					"EQUAL_COUPLING":True,
					"MODEL":"SYBILL",
			},

		
	}





