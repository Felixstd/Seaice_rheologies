import numpy as np 

#-------------------------------------
#------ sea ice characteristics ------
rhoi = 900
Pstar = 5e4
d_mean = 1000
var_state = True

#-------------------------------------
#--------- Ice strength --------------
C = 20

#------------------------------------
#-------- For standard VP -----------
ellip_ratio      = 2
ellip_ratio2     = ellip_ratio**2
ellip_ratio2_inv = 1/ellip_ratio2

#------------------------------------
#-------- For GRAM-mu(I) ------------

mu0       = 0.1  # friction in quasi-static
muinf     = 0.9  # friction in gaseous state
I0        = 1e-3 # critical inertial number 
eII_max   = 1e-4 # maximum shear strain
delta_max = 0.4 # maximum dilatancy
delta_c   = 0
dilat_type = 'const'


#------------------------------------
# Viscous capping
max_zeta = 1e12