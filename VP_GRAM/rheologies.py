import numpy as np
import matplotlib.pyplot as plt

import vp_gram_rheo as rh


# Starting by generating data
data = rh.generate_data(1e-1,
                        [0.5,1], 
                        [0.8, 1], 
                        201)

# Computing strain variables 
rh.compute_strain_var(data=data)


rheos = ['svp']

for rheo_n in rheos:
    rh.elliptical_vp(data = data, rheo_n = rheo_n)

    rh.compute_sigma(data = data, rheo_n = rheo_n)
    
    rh.compute_stress_inv(data = data, rheo_n = rheo_n)
    
    rh.plot_stresses(data = data, rheo_n = rheo_n, arrows=True)




