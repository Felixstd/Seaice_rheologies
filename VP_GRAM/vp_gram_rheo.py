import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import Ellipse
from vp_gram_param import *


plt.style.use('science.mplstyle')
def compute_strain(u, v, dx):
    
    # deformation vectors
    dudx = (u[1:,1:] - u[:-1,1:])/dx
    dvdx = (v[1:,1:] - v[:-1,1:])/dx
    dudy = (u[1:,1:] - u[1:,:-1])/dx
    dvdy = (v[1:,1:] - v[1:,:-1])/dx
    
    e11 = dudx
    e22 = dvdy 
    
    e12 = 0.5*(dudy + dvdx)   
    
    return e11, e22, e12
    

def generate_data(type_field, 
                  eps_def, 
                  eps_h, 
                  eps_A,
                  eps_P,
                  size, 
                  dx = 1000):
    """
    Generating random data to compute stresses for the 
    new rheologies. 

    Args:
        eps_def (_type_): _description_
        eps_h (_type_): _description_
        size (_type_): _description_
    """
    
    #velocity vectors 
    
    if type_field == 'rand':
        u = np.random.random((size,size))*eps_def[0] + eps_def[1]
        v = np.random.random((size,size))*eps_def[0] + eps_def[1]
        
    elif type_field == 'pure_shear':
        
        nx, ny = size, size
        x = np.arange(nx) * dx
        y = np.arange(ny) * dx
        X, Y = np.meshgrid(x, y, indexing='ij') 
        amp = eps_def[0]
        
        u = amp*X
        v = -amp*Y
    
    #ice state vectors
    if var_state:
        h = np.random.random((size-1,size-1))*eps_h[0]+eps_h[1]
        A = np.random.random((size-1,size-1))*eps_A[0]+eps_A[1]
        P_rand = np.random.random((size-1,size-1))*eps_P[0]+eps_P[1]
    else:
        h = np.ones((size-1,size-1))
        A = np.ones((size-1,size-1))
        P_rand = np.ones((size-1,size-1))*Pstar
    

    e11, e22, e12 = compute_strain(u, v, dx)  
    
    
    init_data = {'u':u, 
                 'v':v, 
                 'h':h, 'A':A, 
                 'P_rand':P_rand,
                 'e11':e11, 'e22':e22, 
                 'e12':e12}
    
    return init_data

def compute_strain_var(data):
    
    data['ep'] = data['e11'] + data['e22']
    data['em'] = data['e11'] - data['e22']
    data['e12s'] = data['e12']
    
    data['eTmp'] = np.sqrt(data['em']**2+4*data['e12s']**2)
    
    data['e1'] = 1/2*(data['ep']+data['eTmp'])
    data['e2'] = 1/2*(data['ep']-data['eTmp'])
    
    data['eI']  = data['e1'] + data['e2']
    data['eII'] = data['e1'] - data['e2'] 
    
    
def ice_strength(h, A, Pstar = Pstar, C = C):
    
    return Pstar*h*np.exp(-C*(1-A)) 

def inertial_number(eII, P, h, 
                    rho = rhoi, d = d_mean):
    
    return eII*np.sqrt(rho*h/P)*d 

def friction_coefficient(I, mu0, muinf, I0):
    
    delta_mu = muinf - mu0
    
    return mu0 + delta_mu/(I/I0+1)

def dilatancy(eII, eII_max = eII_max, delta_max = delta_max):
    
    eII_ratio = eII/eII_max
    
    return delta_max*eII_ratio*np.exp(1-eII_ratio)
    

def elliptical_vp(data, rheo_n):
    """
    Function to compute pressure and 
    viscosities for the standard VP rheology

    Args:
        data (_type_): _description_
    """
    
    #renaming variables
    em = data['em']
    ep = data['ep']
    e12 = data['e12']
    eII = data['eII']
    h = data['h']
    A = data['A']
    data[rheo_n] = {}    
    
    #computing delta
    delta = np.sqrt(ep**2 + ellip_ratio2_inv*em**2 + 4*ellip_ratio2_inv*e12**2)+1e-12
    
    #computing ice strength
    strength_vp = ice_strength(h, A)
    
    #computing bulk viscosity
    zeta = strength_vp/(2*delta)
    zeta_cap = np.minimum(zeta, max_zeta)
    
    #computing shear viscosity
    eta = zeta*ellip_ratio2_inv
    
    #saving in dictionnary
    data[rheo_n]['zeta'] = zeta_cap
    data[rheo_n]['eta'] = eta
    data[rheo_n]['P'] = strength_vp 
    
def GRAM_dilat(data, rheo_n, dilat = dilat_type):
    
    
    #renaming variables
    em = data['em']
    ep = data['ep']
    e12 = data['e12']
    eII = data['eII']
    h = data['h']
    A = data['A']
    P_rand = data['P_rand']
    data[rheo_n] = {}    
    
    #ice strength
    P = ice_strength(h, A, eII, 'gram')
    P_H = ice_strength(h, A, eII, 'hibler')
    
    mask_p_rand = (P_rand >= P_H)
    P[mask_p_rand] = P_H[mask_p_rand]-P_rand[mask_p_rand]
    
    
    
    
    # inertial number 
    i_num = inertial_number(eII, P, h)
    
    
    # dilatancy
    if dilat == 'const':
        delta = np.ones_like(h)*delta_c
        mu = np.ones_like(h)*10*np.pi/180
        
        print(delta)
    else:
        delta = dilatancy(eII)
        
        # friction coefficient
        mu = np.arctan(friction_coefficient(i_num, mu0, muinf, I0))
    

    #friction angle
    phi = delta + mu
    
    mask_p = (P == P_H)
    e_kk = np.zeros_like(h)
    e_kk[~mask_p] = 2*eII[~mask_p]*np.tan(delta[~mask_p])
    e_kk[mask_p] = ep[mask_p]
    
    eta = np.minimum(np.sin(phi)*P/(eII+1e-12), max_zeta)
    
    zeta = np.zeros_like(eta)

    #saving in dictionnary
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['P'] = P
    data[rheo_n]['ep'] = e_kk
    
        
    return None
    
def compute_viscs(data, rheo_n):
    
    if rheo_n == 'svp':
        elliptical_vp(data = data, rheo_n = rheo_n)
        
    elif rheo_n == 'muGRAM':
        GRAM_dilat(data = data, rheo_n = rheo_n)
        
    
    
    
    
def compute_sigma(data, rheo_n):
    
    
    zeta = data[rheo_n]['zeta']
    eta  = data[rheo_n]['eta']
    press = data[rheo_n]['P']
    em = data['em']
    ep = data['ep']
    e12 = data['e12']
    e11 = data['e11']
    e22 = data['e22']
    
    if rheo_n == 'svp':
        
        data[rheo_n]['sig11'] = zeta*ep + eta*em - press/2
        data[rheo_n]['sig22'] = zeta*ep - eta*em - press/2
        data[rheo_n]['sig12'] = 2*e12*eta
        
    elif rheo_n == 'muGRAM':
        
        ekk_delta = data[rheo_n]['ep']
        
        data[rheo_n]['sig11'] = 2*eta*e11 - eta*ekk_delta - press
        data[rheo_n]['sig22'] = 2*eta*e22 - eta*ekk_delta - press
        data[rheo_n]['sig12'] = 2*e12*eta
        
def compute_stress_inv(data, rheo_n):
    
    sig11 = data[rheo_n]['sig11'] 
    sig22 = data[rheo_n]['sig22'] 
    sig12 = data[rheo_n]['sig12']
    
    press = data[rheo_n]['P']
    
    sigp = sig11 + sig22 
    sigm = sig11 - sig22 
    
    sigtmp = np.sqrt(sigm**2+4*sig12**2)
    
    #--------------------------------------------------
    # computing principal stresses 
    sig1 = data[rheo_n]['sig1'] = 1/2*(sigp + sigtmp)
    sig2 = data[rheo_n]['sig2'] = 1/2*(sigp - sigtmp)
    
    # normalized
    sig1n = sig1/press
    sig2n = sig2/press
    
    data[rheo_n]['sig1n'] = sig1n
    data[rheo_n]['sig2n'] = sig2n
    
    #--------------------------------------------------
    # computing stress invariants
    data[rheo_n]['sigI'] = 0.5*(sig1+sig2)
    data[rheo_n]['sigII'] = 0.5*(sig1-sig2)

    #normalized
    data[rheo_n]['sigIn'] = 0.5*(sig1n+sig2n)
    data[rheo_n]['sigIIn'] = 0.5*(sig1n-sig2n)
    
    
def plot_stresses(data, rheo_n, arrows = False):
    
    
    sigIn  = data[rheo_n]['sigIn']
    sigIIn = data[rheo_n]['sigIIn']
    sigI  = data[rheo_n]['sigI']
    sigII = data[rheo_n]['sigII']
    sig1   = data[rheo_n]['sig1']
    sig2  = data[rheo_n]['sig2']
    eI     = data['eI']
    eII    = data['eII']
    
    
    #----------------------------------------------------
    # Stress invariants 
    
    fig, ax = plt.subplots(figsize = (8, 6))

    plt.scatter(sigI.flatten(), sigII.flatten(), color = 'royalblue')
    
    if rheo_n == 'svp':
        elli=Ellipse(xy=((-1)/2.,0), width=(1)/ellip_ratio, height=(1), angle=-90,edgecolor='b', fc='None', lw=0.5, color = 'k')
        ax.add_patch(elli)
        plt.ylim(-0.25,0.4)
        plt.xlim(-1.2,0.2)
        
    if arrows :
        qpfac=10

        eu=np.hypot(eI[::qpfac,::qpfac],eII[::qpfac,::qpfac])
        
        ax.quiver(sigI[::qpfac,::qpfac],sigII[::qpfac,::qpfac],
                  eI[::qpfac,::qpfac]/eu,eII[::qpfac,::qpfac]/eu, 
                  scale=10, color='darkgreen')
    
    plt.xlabel(r'$\sigma_\mathrm{I}/P$')
    plt.ylabel(r'$\sigma_\mathrm{II}/P$')
    
    plt.savefig('stress_inv_{}.png'.format(rheo_n), dpi = 500)
    
    fig, ax = plt.subplots(figsize = (8, 6))
    
    plt.scatter(sig2.flatten(), sig1.flatten(), color = 'royalblue',zorder=1)
    plt.grid(zorder = 0)
    
    plt.xlabel(r'$\sigma_2$')
    plt.ylabel(r'$\sigma_1$')
    
    plt.savefig('stress_princ_{}.png'.format(rheo_n), dpi = 500)
