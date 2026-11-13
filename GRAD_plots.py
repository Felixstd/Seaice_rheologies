import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

plt.style.use('./VP_GRAM/science.mplstyle')

def friction_coefficient(I, mu0 = 0.1, muinf = 0.9, I0 = 1e-3):
    
    delta_mu = muinf-mu0
    
    return mu0+delta_mu/(I0/I+1)


def yield_curve(sigIn, kt, mu, c, hyperbolic=False):
    gamma = 2*(1+kt)*np.sqrt((1+kt)/3)
    
    if hyperbolic:

        y = np.sqrt(-c*kt + (mu/gamma**2)*(-sigIn+kt)**(2)*(sigIn+1))
        
    else:
        y = (mu/(gamma**2))*(-sigIn+kt)*(sigIn+1)**(1/2)
    
    return y 


def plastic_potential(sigIn, kt, delta):
    
    gamma = 2*(1+kt)*np.sqrt((1+kt)/3)
    
    y = delta/(gamma**2)*(-sigIn+kt)*(sigIn+1)**(1/2)
    
    return y 

#------- Constants -----

mu0   = 0.1  # static friction
muinf = 0.9  # dynamic friction
mu    = 0.7  # constant friction
delta = 0.5
kt    = 0.1 # tensile
c     = 0

# inertial number
iner_num = np.logspace(-7,0, 1000)

sigIn = np.linspace(-1, 1, 10000)

#------- Variables -------
mu_I = friction_coefficient(iner_num, 
                            mu0, muinf)

F = yield_curve(sigIn, kt, mu, 0, hyperbolic=False)
G = plastic_potential(sigIn, kt, delta)
F[F<0] = np.nan
F[sigIn > kt] = np.nan

G[G<0] = np.nan
G[sigIn > kt] = np.nan



#-----------------------------------------------------
#---------------------- Figures ---------------------- 

#---------------- Friction Coefficient ---------------
plt.figure()
ax = plt.axes()

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.plot(iner_num, mu_I, color = 'k')

plt.axvspan(1e-7, 1e-3, color='dodgerblue', alpha=0.4, lw=0)
plt.axvspan(1e-3, 1, color='goldenrod', alpha=0.2, lw=0)

plt.xscale('log')
plt.yticks([mu0, muinf], labels = [r'$\mu_0$', r'$\mu_\infty$'])

plt.xlabel(r'$I$')
plt.ylabel(r'$\mu$', rotation = 0)

plt.savefig('inertial_num.png')


#---------------- Yield Curve ---------------

plt.figure(figsize = (6,2.5))
ax = plt.axes()
ax.spines['left'].set_position(('data', 0))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
ax.set_ylabel(r'$\frac{\sigma_\mathrm{II}}{P}$', rotation=0)
ax.yaxis.set_label_coords(0.05, 0.95, transform=trans)

plt.yticks([0.10, 0.2])
plt.ylim(0, 0.2)

plt.plot(sigIn, F, color = 'royalblue')

plt.plot(sigIn, G, color = 'darkred')

plt.grid(alpha = 0.5)
plt.xlabel(r'$\frac{\sigma_\mathrm{I}}{P}$')

plt.savefig('yield_curve.png')