import numpy as np
import matplotlib.pyplot as plt
import xarray as xr 

datadir = '/aos/home/fstdenis/Seaice_Rheologies/Data/' 

ds = xr.open_dataset(datadir+'RGPS_composite_deformation_2008.nc')

dvdx = np.array(ds['dvdx'])*3600*24
dvdy = np.array(ds['dvdy'])*3600*24
dudx = np.array(ds['dudx'])*3600*24
dudy = np.array(ds['dudy'])*3600*24

print(np.shape(dudx))

div = dudx+dvdy
shear = np.sqrt((dudx - dvdy)**2 + (dudy +dvdx)**2)

eps_tot = np.sqrt(div**2 + shear **2)

eps = 1e-8

mask = np.isfinite(shear) & (eps_tot > eps) & (eps_tot < 4)

print(np.max(shear[mask]))

dilat = np.logspace(-9, 1, 100000)
print(dilat)

ratio = div[mask].flatten()/shear[mask].flatten()
idx = np.where(ratio < 100)[0]

plt.figure(figsize = (6,4))
cm = plt.scatter(shear[mask].flatten(), div[mask].flatten(), c=eps_tot[mask].flatten(),
    cmap="viridis")
plt.plot(dilat, dilat)
plt.plot(dilat, -dilat)
plt.xlabel(r'$\dot{\epsilon}_\mathrm{II}$ (day$^{-1}$)')
plt.ylabel(r'$\dot{\epsilon}_\mathrm{I}$ (day$^{-1}$)')
plt.grid(alpha = 0.5, zorder = 0)
plt.xscale('log')
plt.yscale('symlog', linthresh = 1e-9)
plt.colorbar(cm)
plt.savefig('dilat.png', dpi = 300)



# print(np.shape(ratio[index(x)]))
plt.figure()
plt.hist(ratio[idx], bins = 1000, density = True)
# plt.xlim(-20, 20)
# plt.scatter(shear[mask].flatten(), div[mask].flatten())
# plt.plot(dilat, dilat)
# plt.plot(dilat, -dilat)
# plt.xscale('log')
# plt.yscale('log')
plt.savefig('ratio.png', dpi = 300)
# plt.figure()
# plt.plot(ratio[mask_ratio])
# # plt.xlim(-20, 20)
# # plt.scatter(shear[mask].flatten(), div[mask].flatten())
# # plt.plot(dilat, dilat)
# # plt.plot(dilat, -dilat)
# # plt.xscale('log')
# # plt.yscale('log')
# plt.savefig('ratio.png', dpi = 300)

