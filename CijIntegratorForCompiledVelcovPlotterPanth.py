import cosmolopy as cp
import numpy as np
from scipy import integrate
from scipy.special import legendre, spherical_jn
import itertools
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.cosmology import Planck15 as cosmo
import sys


start, go = int(sys.argv[1]), int(sys.argv[2])

#K = np.load('../Dipole_JLA/SNMLE/JLADirZInc.npy');
#jlarr = np.genfromtxt('../Dipole_JLA/jla_likelihood_v6/data/jla_lcparams.txt', skip_header=1)

jlarr = np.genfromtxt('../Dipole_JLA/Pantheon/Pantheon_G10.txt')
snnames = np.asarray([k.split()[1] for k in open('../Dipole_JLA/Pantheon/Pantheon_G10.txt').readlines()])


#zcmb, zhel = K[:,8], K[:,9]
#ras, decs = K[:,6], K[:,7]
#ras[ras<0.] = ras[ras<0.] + 360.

zcmb = jlarr.transpose()[7]
ras, decs = jlarr.transpose()[34], jlarr.transpose()[35]

#snnames = snnames[zcmb<0.5]
ras, decs = ras[zcmb<0.5], decs[zcmb<0.5]


lowzindices = []

for i in range(len(zcmb[zcmb<0.5])):
        lowzindices.append(i)

zcmb = zcmb[zcmb<0.5]
zhel = zcmb

combs = np.asarray(list(itertools.combinations(lowzindices, 2)))

def cdAngle(ra1, dec1, ra2, dec2):
    return np.cos(np.deg2rad(dec1))*np.cos(np.deg2rad(dec2))*np.cos(np.deg2rad(ra1) - np.deg2rad(ra2))+np.sin(np.deg2rad(dec1))*np.sin(np.deg2rad(dec2))



cpars = cp.parameters.WMAP7_BAO_H0_mean()

#cpars['N_nu'] = cosmo.Neff
cpars['omega_M_0'] = cosmo.Om0
cpars['omega_lambda_0'] = cosmo.Ode0

ks = np.power(10., np.linspace(-6, 0, 601))
pks = cp.perturbation.power_spectrum(ks, 0, **cpars)
PofK = InterpolatedUnivariateSpline(ks, pks)

def Pl(ctheta, l):
    return legendre(l)(ctheta)


def jprimel(kchi, l):
    return spherical_jn(n=l, z=kchi, derivative=True)

def CovMatCalculator(z1s, z2s, cthetas, lmax=15):
    ch1s = cosmo.comoving_distance(z1s).value
    ch2s = cosmo.comoving_distance(z2s).value
    def PVCovMatIntegrand(k, chi1, chi2, ctheta):
        sumand=0
        for l in range(lmax):
            term = (2.*l+1.)*jprimel(k*chi1, l)*jprimel(k*chi2, l)*Pl(ctheta, l)
            sumand+=term
            #print l, term
        return PofK(k)*sumand
    retarr=[]
    for a,b,c in zip(ch1s, ch2s, cthetas):
        retarr.append(integrate.quad(PVCovMatIntegrand, 10**-6., 1.0, args=(a,b,c))[0])
        print 'Done', len(retarr)
    return np.asarray(retarr)
    #return np.asarray([integrate.quad(PVCovMatIntegrand, 10**-6., 1.0, args=(a,b,c))[0] for a,b,c in zip(ch1s, ch2s, cthetas)])

z1s = np.asarray([zhel[i] for i, j in combs])
z2s = np.asarray([zhel[j] for i, j in combs])
cthetas = np.asarray([cdAngle(ras[i], decs[i], ras[j], decs[j]) for i, j in combs])

dangles = np.rad2deg(np.arccos(cthetas))
zdiffs = np.absolute(z1s - z2s)


#sel1 = (zdiffs>0.03)*(dangles<3)
#sel2 = (zdiffs<0.03)*(dangles<3)
#sel3 = (dangles>170.)*(zdiffs<0.03)
#sel4 = (dangles>170.)*(zdiffs>0.06)


#fsel = sel1+sel2+sel3+sel4

#z1s, z2s, cthetas = z1s[fsel], z2s[fsel], cthetas[fsel]

print 'Grand Total', len(combs), 'couples'

print 'Total', len(z1s), 'couples'

Cijint = CovMatCalculator(z1s[start:start+go], z2s[start:start+go], cthetas[start:start+go])

np.savetxt('CijIntegralsPanthBF/Cijintsaved_'+str(start)+'_'+str(start+go)+'.txt', Cijint, delimiter='|')
