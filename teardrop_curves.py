"""
This is code to plot the different teardrop yield curves. 

There are a few variats possible:
- with mu
- with hyperbolicity 
- with max not changing

All of the yield curves will be plotted in normalized invariant space:
sigI/P, sigII/P
"""

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

from matplotlib.lines import Line2D

plt.style.use('./VP_GRAM/science.mplstyle')


def teardrop_yc(type, sigIn, mu, kt):
    
    if type == 'base':
        signIIn = (-mu*sigIn+kt)*(sigIn+1)**(1/2)
        
    elif type == 'kt_fixed':
        signIIn = (-mu*sigIn+mu*kt)*((sigIn+1))**(1/2)
        
    elif type == 'max_fixed':
        gamma = (2+2*kt)*((1+kt)/3)**(1/2)
        signIIn = ((-mu*sigIn+mu*kt)*((sigIn+1))**(1/2))/gamma
    
    elif type == 'hyperbolic': 
        c = 1e-3
        signIIn = (-c*kt+((-mu*sigIn+mu*kt)**2*((sigIn+1))))**(1/2)

    signIIn[signIIn<0] = np.nan
    return signIIn

#----- Constants -------#

n = 10000  #number of data points
mu = 0.7  #friction coefficient 
friction_coefficents = np.array([0.2, 0.7])[:, None]  
# friction_coefficents[0,:] *= 0.2
# friction_coefficents[1,:] *= 0.7
kts = np.array([0.05, 0.1])[:, None]   #tension coefficient 
kt = kts[0]

#setting the compressive stress
sigIn = np.linspace(-1, 1, n)

#----- Yield curves -------#

signIIn_base = teardrop_yc('base', sigIn, friction_coefficents, kt)
signIIn_kt = teardrop_yc('kt_fixed', sigIn, friction_coefficents, kt)
signIIn_max = teardrop_yc('max_fixed', sigIn, friction_coefficents, kt)
signIIn_hyper = teardrop_yc('hyperbolic', sigIn, friction_coefficents, kt)

signIIn_base_kt = teardrop_yc('base', sigIn, friction_coefficents[1], kts)
signIIn_kt_kt = teardrop_yc('kt_fixed', sigIn, friction_coefficents[1], kts)
signIIn_max_kt = teardrop_yc('max_fixed', sigIn, friction_coefficents[1], kts)

#----- Figures -------#

plt.figure(figsize = (6,2))
ax = plt.axes()

# for mu in friction_coefficents:
for row, m, linestyle in zip(signIIn_base, friction_coefficents.ravel(), ['-', '--']):
    ax.plot(sigIn, row, color='goldenrod', linestyle = linestyle)
    
# for row, m, linestyle in zip(signIIn_hyper, friction_coefficents.ravel(), ['-', '--']):
#     ax.plot(sigIn, row, color='darkgreen', linestyle = linestyle)

for row, m, linestyle in zip(signIIn_kt, friction_coefficents.ravel(), ['-', '--']):
    ax.plot(sigIn, row, color='blue', linestyle = linestyle)
    
for row, m, linestyle in zip(signIIn_max, friction_coefficents.ravel(), ['-', '--']):
    ax.plot(sigIn, row, color='crimson', linestyle = linestyle)



# put the horizontal axis line at y = 0
ax.spines['left'].set_position(('data', 0))
# hide the top and right box edges
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
ax.set_ylabel(r'$\frac{\sigma_\mathrm{II}}{P}$', rotation=0)
ax.yaxis.set_label_coords(0.05, 0.95, transform=trans)
ax.set_xticklabels([])
ax.set_yticklabels([])


handles = [
    Line2D([0], [0], color='goldenrod', lw=1),
    Line2D([0], [0], color='blue',      lw=1),
    Line2D([0], [0], color='crimson',   lw=1),
]
labels = ['Zhang', r'var-$k_t$', 'Max fixed']

ax.legend(handles, labels,bbox_to_anchor=(0.8, 0.25))

plt.grid(alpha = 0.5)
plt.xlabel(r'$\frac{\sigma_\mathrm{I}}{P}$')
# plt.title('Sensitivity to friction coefficient $\mu$')

plt.savefig('teardrop_ycs_musensitivity.png')


plt.figure(figsize = (6,2))
ax = plt.axes()

# for mu in friction_coefficents:
for row, m, linestyle in zip(signIIn_base_kt, friction_coefficents.ravel(), ['-', '--']):
    ax.plot(sigIn, row, color='goldenrod', linestyle = linestyle)

for row, m, linestyle in zip(signIIn_kt_kt, friction_coefficents.ravel(), ['-', '--']):
    ax.plot(sigIn, row, color='blue', linestyle = linestyle)
    
for row, m, linestyle in zip(signIIn_max_kt, friction_coefficents.ravel(), ['-', '--']):
    ax.plot(sigIn, row, color='crimson', linestyle = linestyle)



# put the horizontal axis line at y = 0
ax.spines['left'].set_position(('data', 0))
# hide the top and right box edges
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
ax.set_ylabel(r'$\frac{\sigma_\mathrm{II}}{P}$', rotation=0)
ax.yaxis.set_label_coords(0.05, 0.95, transform=trans)
ax.set_xticklabels([])
ax.set_yticklabels([])


handles = [
    Line2D([0], [0], color='goldenrod', lw=1),
    Line2D([0], [0], color='blue',      lw=1),
    Line2D([0], [0], color='crimson',   lw=1),
]
labels = ['Zhang', r'var-$k_t$', 'Max fixed']

ax.legend(handles, labels,bbox_to_anchor=(0.9, 0.75))

plt.grid(alpha = 0.5)
plt.xlabel(r'$\frac{\sigma_\mathrm{I}}{P}$')
# plt.title('Sensitivity to friction coefficient $\mu$')

plt.savefig('teardrop_ycs_ktsensitivity.png')