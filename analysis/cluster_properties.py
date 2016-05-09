"""
Results:
Implied stellar mass of e1 cluster, assuming 7 stars are OB stars: 683.1740479924961 +/- 216.68608280894213 M_sun
Implied stellar mass of e1 cluster, assuming 7 stars are O stars: 2325.0125309614355 +/- 807.0686206062493 M_sun
Implied stellar mass of e1e2 cluster, assuming 12 stars are O stars: 4187.004606488459 +/- 1172.7498994742245 M_sun
Implied stellar mass of e1e2IRS2 cluster, assuming 18 stars are O stars: 5845.065779607007 +/- 1266.3844454611049 M_sun
emcee estimates given L=2e7+/-5e6
Cluster mass: 6544.185646762511 \pm 2370.8187171691534
N(Ostars): 18.92140625 \pm 6.2906031193991465
"""

import numpy as np
from imf import imf
import pylab as pl
from astropy.utils.console import ProgressBar
import paths
pl.matplotlib.rc_file('pubfiguresrc')

nclusters_per_bin = 100
maxmass = 3e4
mmin_bstar = 8
mmin_ostar = 20

cl_masses = np.logspace(2.5,np.log10(maxmass))
clusters = np.array([[imf.make_cluster(mass*(np.random.randn()/20.+1.),silent=True) for ii in range(nclusters_per_bin)] for mass in ProgressBar(cl_masses)])
luminosities = [imf.lum_of_cluster(cl) for section in ProgressBar(clusters) for cl in section]
#cl_masses_flat = np.array([x for mass in cl_masses for x in [mass for ii in range(nclusters_per_bin)]])
cl_masses_flat = np.array([sum(cl) for section in clusters for cl in section])
n_obstars = np.array([(cl>mmin_bstar).sum() for section in clusters for cl in section])
n_ostars = np.array([(cl>mmin_ostar).sum() for section in clusters for cl in section])

# guess, not measurement
nobstars_w51 = n_obstars[(cl_masses_flat > 4e3) & (cl_masses_flat < 9e3)]
nostars_w51 = n_ostars[(cl_masses_flat > 4e3) & (cl_masses_flat < 9e3)]


print("Implied stellar mass of e1 cluster, assuming 7 stars are OB stars: {0} +/- {1} M_sun".format(cl_masses_flat[n_obstars==7].mean(), cl_masses_flat[n_obstars==7].std()))
print("Implied stellar mass of e1 cluster, assuming 7 stars are O stars: {0} +/- {1} M_sun".format(cl_masses_flat[n_ostars==7].mean(), cl_masses_flat[n_ostars==7].std()))
print("Implied stellar mass of e1 cluster, assuming 6 stars are OB stars: {0} +/- {1} M_sun".format(cl_masses_flat[n_obstars==6].mean(), cl_masses_flat[n_obstars==6].std()))
print("Implied stellar mass of e1 cluster, assuming 6 stars are O stars: {0} +/- {1} M_sun".format(cl_masses_flat[n_ostars==6].mean(), cl_masses_flat[n_ostars==6].std()))
print("Implied stellar mass of e1e2 cluster, assuming 12 stars are O stars: {0} +/- {1} M_sun".format(cl_masses_flat[n_ostars==12].mean(), cl_masses_flat[n_ostars==12].std()))
print("Implied stellar mass of e1e2IRS2 cluster, assuming 18 stars are O stars: {0} +/- {1} M_sun".format(cl_masses_flat[n_ostars==18].mean(), cl_masses_flat[n_ostars==18].std()))

pl.figure(6).clf()
pl.semilogx(cl_masses_flat, luminosities, '.', alpha=0.5)
pl.axhline(np.log10(8.3e6), linewidth=10, alpha=0.1, zorder=-5)
pl.axhline(np.log10(2e7), linewidth=10, alpha=0.1, zorder=-5)
pl.xlim(10**2.5, (maxmass))
pl.ylim(5.5, 8.1)
pl.xlabel("Cluster mass ($M_\odot$)")
pl.ylabel("Cluster luminosity (log $L_\odot$)")
pl.savefig(paths.fpath("clusters/cluster_mass_vs_luminosity.png"))

pl.figure(9).clf()
pl.loglog(cl_masses_flat, n_ostars, '.', alpha=0.5, label='O-stars $(M>{0} M_\odot)$'.format(mmin_ostar))
pl.loglog(cl_masses_flat, n_obstars, '.', alpha=0.5, label='B-stars $(M>{0} M_\odot)$'.format(mmin_bstar))
pl.xlim(10**2.5, (maxmass))
pl.ylim(0, n_ostars.max())
pl.ylim(0, n_obstars.max())
pl.xlabel("Cluster mass ($M_\odot$)")
pl.ylabel("Number of OB-stars")
pl.legend(loc='best')
pl.savefig(paths.fpath("clusters/cluster_mass_vs_n_ostars.png"))


def lnprob(mass, luminosity=2e7, luminosity_error=5e6):
    if mass < 1:
        return -np.inf, 0
    if mass > 1e5:
        # use a best-fit relation; scatter is minimal
        lum = 10**(np.log10(mass) + 3.55)
        # n > 20 msun
        nostars = 10**(np.log10(mass) - 2.5)
    else:
        cluster = imf.make_cluster(mass, silent=True)
        lum = 10**imf.lum_of_cluster(cluster)
        nostars = (cluster > mmin_ostar).sum()
    return -(luminosity-lum)**2/(2*luminosity_error**2), nostars

import emcee
ndim = 1
nwalkers = 48
p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))*1000 + 5000
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=4)
pos,prob,state,blobs = sampler.run_mcmc(p0, 400)




pl.figure(7).clf()
pl.hist(nostars_w51, bins=30)

pl.figure(8)
pl.clf()
pl.hist(sampler.flatchain, bins=100)
pl.xlabel("Cluster Mass at $L=2 \pm 0.5\\times10^7$ L$_\odot$")

print("emcee estimates given L=2e7+/-5e6")
print("Cluster mass: {0} \pm {1}".format(sampler.flatchain.mean(),
                                         sampler.flatchain.std()))

pl.figure(9)
pl.clf()
flatblobs = np.array([x for y in sampler.blobs for x in y])
a,b,c = pl.hist(flatblobs, bins=np.arange(0,50,1))
x = np.linspace(0,50,500)
pl.plot(x, a.max()*np.exp(-(x-flatblobs.mean())**2 / (2*flatblobs.var())), '-', linewidth=2, alpha=0.5)
pl.xlabel("$N(M>{0}$M$_\odot)$ at $L=2 \pm 0.5\\times10^7$ L$_\odot$".format(mmin_ostar))

print("N(Ostars): {0} \pm {1}".format(flatblobs.mean(),
                                      flatblobs.std()))
