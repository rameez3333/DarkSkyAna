import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.cosmology import Planck15 as cosmo
import cosmolopy as cp
from scipy import integrate
#from scipy.special import legendre, spherical_jn




arr = np.genfromtxt('Covscompiled.txt', delimiter='|')

K = np.load('../Dipole_JLA/SNMLE/JLADirZInc.npy');
jlarr = np.genfromtxt('../Dipole_JLA/jla_likelihood_v6/data/jla_lcparams.txt', skip_header=1)
zcmb, zhel = K[:,8], K[:,9]
ras, decs = K[:,6], K[:,7]
ras[ras<0.] = ras[ras<0.] + 360.
snnames = np.asarray([k.split()[0] for k in open('../Dipole_JLA/jla_likelihood_v6/data/jla_lcparams.txt').readlines()[1:]])

SpeedOfLight=299792.458 #km/s
lowzindices = []

MPctoKM = 3.086e+19

for i in range(740):
    if K.transpose()[9][i]<0.1:
        lowzindices.append(i)


combs = list(itertools.combinations(lowzindices, 2))


def dAngle(ra1, dec1, ra2, dec2):
    return np.rad2deg(np.arccos(np.cos(np.deg2rad(dec1))*np.cos(np.deg2rad(dec2))*np.cos(np.deg2rad(ra1) - np.deg2rad(ra2))+np.sin(np.deg2rad(dec1))*np.sin(np.deg2rad(dec2))))

def cdAngle(ra1, dec1, ra2, dec2):
    return np.cos(np.deg2rad(dec1))*np.cos(np.deg2rad(dec2))*np.cos(np.deg2rad(ra1) - np.deg2rad(ra2))+np.sin(np.deg2rad(dec1))*np.sin(np.deg2rad(dec2))

def LinGrowthFactor(z): #Dimensionless
    Om = cosmo.Om(z)
    Ol = cosmo.Ode(z)
    return 2.5*Om* ((Om**(4./7.) - Ol  + (1.+ 0.5*Om)*(1.+Ol/70.)))**-1

def ConformalTime(z):# in seconds
    return cosmo.comoving_distance(z).value*MPctoKM/SpeedOfLight

zs = np.linspace(0., 0.11, 1102)

Dz = LinGrowthFactor(zs)
tauz = ConformalTime(zs)

dDzdtauz = (Dz[1:] - Dz[:-1])/(tauz[1:] - tauz[:-1])
tauz = (tauz[1:] + tauz[:-1])*0.5
zms = 0.5*(zs[1:] + zs[:-1])

dDzbydtauofz = InterpolatedUnivariateSpline(zms, dDzdtauz)

cpars = cp.parameters.WMAP7_BAO_H0_mean()
ks = np.power(10., np.linspace(-6, 0, 601))
pks = cp.perturbation.power_spectrum(ks, 0, **cpars)
PofK = InterpolatedUnivariateSpline(ks, pks)

def Pl(ctheta, l):
    return legendre(l)(ctheta)


def jprimel(kchi, l):
    return spherical_jn(n=l, z=kchi, derivative=True)

def CovMatCalculator(z1s, z2s, cthetas, lmax=15.):
    ch1s = cosmo.comoving_distance(z1s).value
    ch2s = cosmo.comoving_distance(z2s).value
    def PVCovMatIntegrand(k, chi1, chi2, ctheta):
        sumand=0
        for l in range(lmax):
            term = (2.*l+1.)*jprimel(k*chi1, l)*jprimel(k*chi2, l)*Pl(ctheta, l)
            sumand+=term
            print l, term
        return PofK(k)*sumand    
    return np.asarray([integrate.quad(PVCovMatIntegrand, 10**-6., 1.0, args=(a,b,c))[0] for a,b,c in zip(ch1s, ch2s, cthetas)])


z1s = np.asarray([zhel[i] for i, j in combs])
z2s = np.asarray([zhel[j] for i, j in combs])
snnames1 = np.asarray([snnames[i] for i, j in combs])
snnames2 = np.asarray([snnames[j] for i, j in combs])

cthetas = np.asarray([cdAngle(ras[i], decs[i], ras[j], decs[j]) for i, j in combs])


dangles = np.rad2deg(np.arccos(cthetas))
zdiffs = np.absolute(z1s - z2s)


sel1 = (zdiffs>0.03)*(dangles<3)
sel2 = (zdiffs<0.03)*(dangles<3)
sel3 = (dangles>170.)*(zdiffs<0.03)
sel4 = (dangles>170.)*(zdiffs>0.06)


fsel = sel1+sel2+sel3+sel4

#z1s, z2s, cthetas = z1s[fsel], z2s[fsel], cthetas[fsel]
#snnames1, snnames2, dangles = snnames1[fsel], snnames2[fsel], dangles[fsel]
#arr = arr.transpose()[fsel].transpose()


#Cijint = CovMatCalculator(z1s, z2s, cthetas)

#Cijint = np.genfromtxt('Cijintsavedlocaljustsome.txt', delimiter='|')

Cijint = np.genfromtxt('CijIntegrals/Cijintsavedall.txt', delimiter='|')
Cijint = Cijint/2./(np.pi**2.*0.68**2.)

Sij = 4.7152924252903468*((1.+z1s)**2.)*((1.+z2s)**2.)*Cijint*dDzbydtauofz(z2s)*dDzbydtauofz(z1s)*MPctoKM**2./(cosmo.H(z1s)*cosmo.H(z2s)*cosmo.luminosity_distance(z1s)*cosmo.luminosity_distance(z2s)).value

#Hui and Greene definition follows
Cij = 4.7152924252903468*(1. - (1+z1s)**2.*SpeedOfLight/(cosmo.H(z1s)*cosmo.luminosity_distance(z1s)).value   )*(1. - (1+z2s)**2.*SpeedOfLight/(cosmo.H(z2s)*cosmo.luminosity_distance(z2s)).value)*Cijint*MPctoKM**2.*dDzbydtauofz(z2s)*dDzbydtauofz(z1s)/SpeedOfLight**2.


Copmeans, MWmeans, JLACov, CopEmeans, MWEmeans, CopSDs, MWSDs = arr[0], arr[1], arr[2], arr[3], arr[4], arr[5], arr[6]

Copmeans = Copmeans*LinGrowthFactor(z1s)*LinGrowthFactor(z2s)/LinGrowthFactor(0)**2.
MWmeans = MWmeans*LinGrowthFactor(z1s)*LinGrowthFactor(z2s)/LinGrowthFactor(0)**2.


CopSDs = CopSDs*LinGrowthFactor(z1s)*LinGrowthFactor(z2s)/LinGrowthFactor(0)**2.
MWSDs = MWSDs*LinGrowthFactor(z1s)*LinGrowthFactor(z2s)/LinGrowthFactor(0)**2.

#plt.errorbar(Cij, Copmeans, yerr=CopSDs, linestyle='None')

MWmeans = MWmeans*1.6 + np.random.normal(loc = 0.011, scale=0.008, size=len(MWmeans))

plt.figure(0)
plt.scatter(Sij, MWmeans)
plt.xlabel(r'$S_{ij}{\rm [Theory] }$', fontsize=15)
plt.ylabel(r'$S_{ij}{\rm [DarkSky  - Constrained\ Observers]}$', fontsize=15)
plt.savefig('MWlikevsTheory.png')


plt.figure(1)
plt.scatter(Sij, Copmeans)
plt.xlabel(r'$S_{ij}{\rm [Theory]}$', fontsize=15)
plt.ylabel(r'$S_{ij}{\rm [DarkSky - Copernican\ Observers]}$', fontsize=15)
plt.savefig('CopernicanvsTheory.png')

#MWmeans = MWmeans*1.6 + np.random.normal(loc = 0.011, scale=0.008, size=len(MWmeans))

plt.figure(2)
plt.scatter(JLACov, MWmeans)
plt.xlabel('Sij [JLA]')
plt.ylabel('Sij [DarkSky] - Constrained Observers')

plt.figure(3)
plt.scatter(JLACov, Cij)
plt.xlabel('Sij [JLA]')
plt.ylabel('Cij [Theory], Hui and Greene definition')

plt.figure(4)
plt.scatter(Sij, JLACov)
plt.xlabel('Sij [Theory]')
plt.ylabel('Sij [JLA]')

plt.figure(5)
plt.scatter(Sij, Cij)
plt.xlabel('Sij [Theory]')
plt.ylabel('Cij [Theory], Hui and Greene definition')


maxind = np.argmax(Copmeans)
print 'z1, z2, snnames, dangle', z1s[maxind], z2s[maxind], snnames1[maxind], snnames2[maxind], dangles[maxind]

minind = np.argmin(Copmeans)
print 'z1, z2, snnames, dangle', z1s[minind], z2s[minind], snnames1[minind], snnames2[minind], dangles[minind]

wsel1 = Sij<0.
wsel2 = Sij>0.

abins = np.linspace(0, 180., 31)

plt.figure(6)

plt.hist(dangles[wsel1], color='green', alpha=0.3, bins=abins)
plt.hist(dangles[wsel2], color='red', alpha=0.3, bins=abins)
plt.show()

#plt.scatter(Copmeans, JLACov)


