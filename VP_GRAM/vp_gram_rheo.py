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
    data[rheo_n]['P'] = strength_vp/2
    data[rheo_n]['press0'] = strength_vp

def MC_ps(data, rheo_n ):
    
    #renaming variables
    em = data['em']
    ep = data['ep']
    e12 = data['e12']
    eII = data['eII']
    h = data['h']
    A = data['A']
    P_rand = data['P_rand']
    data[rheo_n] = {}    
    
    mu = 0.7
        
    #ice strength
    P_H = ice_strength(h, A)
    
    zeta = np.minimum(P_H/np.fabs(2*ep+1e-20), max_zeta)
    plt.figure()
    plt.pcolormesh(zeta)
    plt.savefig('visc_bulk1.png')
    
    P = 0.5*P_H
    
    plt.figure()
    pc = plt.pcolormesh(P)
    plt.colorbar(pc)
    plt.savefig('P1.png')
    
    eta = np.minimum(mu*(P - zeta*ep)/((eII+1e-20)), max_zeta)
    
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['P'] = P
    data['ep'] = 0
    
def TMC(data, rheo_n):
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
    
    eta_vp = zeta_cap*ellip_ratio2_inv
    eta_MC = mu/(2*eII)*(strength_vp/2-zeta_cap*ep)
    
    eta = np.minimum(eta_vp, eta_MC)
    
    
    data[rheo_n]['zeta'] = zeta_cap
    data[rheo_n]['eta'] = eta
    data[rheo_n]['P'] = strength_vp
    
def MC_TD(data, rheo_n):
    
    em = data['em']
    ep = data['ep']
    e12 = data['e12']
    eII = data['eII']
    h = data['h']
    A = data['A']
    data[rheo_n] = {} 
    
     #constants
    gamma = 2*(1+kt)*np.sqrt((1+kt)/3)
        
    D = np.tan(np.radians(delta_c))
    alpha = 0.95
    
    mu_mac = D+np.tan(mu)
    
    #ice strength
    P = ice_strength(h, A)
    cyc = (2-kt)/3    
   
    
    gam_over_d = gamma/D 
    gam_over_d_2 = gam_over_d**2
    
    #dilatancy ratio
    l = ep/(eII+1e-20)
    
    x = (-(6*(1+kt)-2*gam_over_d_2*l**2) \
            + 2*gam_over_d*l*np.sqrt(gam_over_d_2*l**2 + 3*(1+kt)))/9 + kt
    
    x = np.minimum(x, alpha*kt)
    
    zeta = (x+cyc)*P/(np.copysign(np.maximum(np.fabs(ep), 1e-20),ep)) 
    
    eta = -mu_mac/(2*(eII + 1e-20)*gamma)*(x-kt)*(x+1)**(1/2)*P
    
    press = cyc * P

    ### save in the dictionary
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['P'] = press
    data[rheo_n]['press0'] = P

    
def compute_viscs(data, rheo_n):
    
    if rheo_n == 'svp':
        elliptical_vp(data = data, rheo_n = rheo_n)

    elif rheo_n == 'MC_S':
        MC_ps(data = data, rheo_n = rheo_n)
        
    elif rheo_n == 'TMC':
        TMC(data = data, rheo_n = rheo_n)
        
    elif rheo_n == 'MCTD':
        MC_TD(data = data, rheo_n = rheo_n)
    
    
def compute_sigma(data, rheo_n):
    
    
    zeta = data[rheo_n]['zeta']
    eta  = data[rheo_n]['eta']
    press = data[rheo_n]['P']
    em = data['em']
    ep = data['ep']
    e12 = data['e12']
    
        
    data[rheo_n]['sig11'] = zeta*ep+eta*em-press
    data[rheo_n]['sig22'] = zeta*ep-eta*em-press
    data[rheo_n]['sig12'] = 2*e12*eta
        
def compute_stress_inv(data, rheo_n):
    
    sig11 = data[rheo_n]['sig11'] 
    sig22 = data[rheo_n]['sig22'] 
    sig12 = data[rheo_n]['sig12']
    
    press = data[rheo_n]['press0']
    
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
    
    
def plot_stresses(data, rheo_n, arrows = False, 
                  figdir = '/aos/home/fstdenis/Seaice_Rheologies/'):
    
    
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

    plt.scatter(sigIn.flatten(), sigIIn.flatten(), color = 'royalblue')
    
    if rheo_n == 'svp':
        elli=Ellipse(xy=((-1)/2.,0), width=(1)/ellip_ratio, height=(1), angle=-90,edgecolor='b', fc='None', lw=0.5, color = 'k')
        ax.add_patch(elli)
        # plt.ylim(-0.25,0.4)
        # plt.xlim(-1.2,0.2)
    
    plt.xlabel(r'$\sigma_\mathrm{I}/P$')
    plt.ylabel(r'$\sigma_\mathrm{II}/P$')
        
    if arrows :
        qpfac=10

        eu=np.hypot(eI[::qpfac,::qpfac],eII[::qpfac,::qpfac])
        
        ax.quiver(sigIn[::qpfac,::qpfac],sigIIn[::qpfac,::qpfac],
                  eI[::qpfac,::qpfac]/eu,eII[::qpfac,::qpfac]/eu, 
                  scale=10, color='darkgreen')
    
    
    plt.grid(alpha = 0.5)
    
    plt.savefig(figdir+'FIGURES/stress_inv_{}.png'.format(rheo_n), dpi = 500)
    
    fig, ax = plt.subplots(figsize = (8, 6))
    
    plt.scatter(sig1.flatten(), sig2.flatten(), color = 'royalblue',zorder=1)
    plt.grid(zorder = 0)
    
    plt.xlabel(r'$\sigma_1$')
    plt.ylabel(r'$\sigma_2$')
    
    plt.savefig(figdir+'FIGURES/stress_princ_{}.png'.format(rheo_n), dpi = 500)
