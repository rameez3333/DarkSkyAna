import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from glob import glob

zhistbins = np.linspace(0, 0.6, 61)

gwlines = open('GrowthRateData').readlines()
gf = []
gz = []

for line in gwlines:
    for entry in line.split('}'):
        try:
            gz.append(float(entry.replace('{','').replace(' ','').split(',')[len(entry.replace('{','').replace(' ','').split(','))-2:][0]))
            gf.append(float(entry.replace('{','').replace(' ','').split(',')[len(entry.replace('{','').replace(' ','').split(','))-2:][1]))
        except:
            print 'yr'
plt.figure(0)            
plt.plot(gz, gf)
plt.ylabel('Structure Growth Factor', fontsize=15)
plt.xlabel('redshift', fontsize=15)
plt.savefig('StructureGrowth.png')

GF = InterpolatedUnivariateSpline(np.asarray(gz), np.asarray(gf))

#del plt




#import matplotlib.pyplot as plt


def angbetweenvecs(v1, v2):
    dotprod = np.sum(v1*v2)/np.sqrt(np.sum(v1*v1)*np.sum(v2*v2))
    #print dotprod
    return np.rad2deg(np.arccos(dotprod))
                                                    
def GetRandomPosition(size=1):
    ra = np.random.uniform(0, 360, size=size)
    dec = np.rad2deg(np.arcsin(np.random.uniform(0, 1., size=size))*np.random.choice([1.,-1.], size=size))
    if size==1:
        ra = ra[0]
        dec = dec[0]
    
    return dec,ra


def fileparser(fname, vmin=0, vmax=650., MWLike=False):
    try:
        print 'Processing ', fname
        flines = open(fname).readlines()
        totnhalos = int(flines[0].split(':')[1])
        obs = flines[1][10:].split('|')
        obsvel = np.asarray([float(obs[3]), float(obs[4]), float(obs[5])])
        obsvelstrength = np.sqrt(np.sum(obsvel*obsvel))
        if (obsvelstrength > vmax) or (obsvelstrength < vmin):
            return None, None, None, None
        obsmass = float(obs[6])
        if MWLike:
            if np.log10(obsmass)>12.15:
                return None, None, None, None

        fts=[]
        dips=[]
        nchalos=[]
        zdips=[]
        znchalos=[]
        dipobsangle=[]
        zdipobsangle=[]
        dipstrengths=[]
        zdipstrengths=[]
        print 'Observer Halo Mass:', obsmass
        print 'Observer Halo Velocity:', obsvelstrength
        for fline in flines[4:]:
            zhist = np.asarray([int(x) for x in fline.split('|')[-1].split(',')])

            if len(fline.split('|')) <11:
                cs = np.cumsum(zhist)
                binmax = zhistbins[len(cs[cs<cs[-1]/2])]
            else:
                binmax = float(fline.split('|')[-2])
            if binmax>0.138 and binmax<0.142:
                scalefactor = np.sum(GF(zhistbins[1:len(zhist)+1]+0.015)*zhist)/np.sum(zhist)
                print 'Scalefactor', scalefactor
                fts.append(float(fline.split('|')[0]))
                nchalos.append(int(fline.split('|')[1]))
                dips.append(np.asarray([float(fline.split('|')[2]), float(fline.split('|')[3]), float(fline.split('|')[4])])/float(nchalos[-1]))
                znchalos.append(int(fline.split('|')[5]))
                if (znchalos[-1] < 1000) or (nchalos[-1]<1000):
                    continue
            #print int(fline.split('|')[5])
                try:
                    zdips.append(np.asarray([float(fline.split('|')[6]), float(fline.split('|')[7]), float(fline.split('|')[8])])/float(znchalos[-1]))
                except:
                    print float(znchalos[-1])
                try:
                    dipobsangle.append(angbetweenvecs(dips[-1], obsvel))
                except:
                    print dips[-1], obsvel
                try:
                    zdipobsangle.append(angbetweenvecs(zdips[-1], obsvel))
                except:
                    print zdips[-1], obsvel
                dipole=dips[-1]
                zdipole=zdips[-1]
                dipstrengths.append(np.sqrt(dipole[0]*dipole[0]+dipole[1]*dipole[1]+dipole[2]*dipole[2])*scalefactor)
                zdipstrengths.append(np.sqrt(zdipole[0]*zdipole[0]+zdipole[1]*zdipole[1]+zdipole[2]*zdipole[2])*scalefactor)
                
                print dipstrengths[-1], dipobsangle[-1], zdipstrengths[-1], zdipobsangle[-1], binmax#, dips[-1], zdips[-1], obsvel
                #plt.plot(zhistbins[1:len(zhist)+1]+0.015, zhist, '--')
                #plt.xlabel('redshift proxy')
                #plt.savefig('RedshiftdistroMethod1.png')
                #plt.show()
            
        return np.median(dipstrengths), np.median(dipobsangle), np.median(zdipstrengths), np.median(zdipobsangle)
    except:
        return None, None, None, None


flist=sorted(glob('PreviousOutputs/DarkSkyDipoleQuery_CenterX_*SELECT2MRS0.03_Out.txt'))
flist2=sorted(glob('PreviousOutputs/DarkSkyDipoleQuery_CenterX_*MWSELECT_Out.txt'))
print len(flist)

#dipstrenghts=[]
#dipangles=[]
#zdipstrength=[]
#zdipangles=[]

#for f in flist:
    #ds, da, zds, zda = fileparser(f, vmin=0, vmax=500.)
    #if ds:
       #dipstrenghts.append(ds)
       #dipangles.append(da)
       #zdipstrength.append(zds)
       #zdipangles.append(zda)



hvdipstrenghts=[]
hvdipangles=[]
hvzdipstrength=[]
hvzdipangles=[]

for f in flist:
    ds, da, zds, zda = fileparser(f, vmin=0, vmax=1000., MWLike=True)
    if ds:
       hvdipstrenghts.append(ds)
       hvdipangles.append(da)
       hvzdipstrength.append(zds)
       hvzdipangles.append(zda)

randdec1, randra1 = GetRandomPosition(2000.)
randdec2, randra2 = GetRandomPosition(2000.)
dangles = np.rad2deg(np.arccos(np.cos(np.deg2rad(randdec1))*np.cos(np.deg2rad(randdec2))*np.cos(np.deg2rad(randra1) - np.deg2rad(randra2))+np.sin(np.deg2rad(randdec1))*np.sin(np.deg2rad(randdec2))))


mwhvdipstrenghts=[]
mwhvdipangles=[]
mwhvzdipstrength=[]
mwhvzdipangles=[]

for f in flist2:
    ds, da, zds, zda = fileparser(f, vmin=0, vmax=1000., MWLike=True)
    if ds:
       mwhvdipstrenghts.append(ds)
       mwhvdipangles.append(da)
       mwhvzdipstrength.append(zds)
       mwhvzdipangles.append(zda)












angbins = np.linspace(0, 180., 46)
dipbins = np.linspace(0, 0.03, 46)

plt.figure(1)
#plt.hist(hvdipangles, bins=angbins,  alpha=0.3, color='green', label='MW-2MRS Observer')
plt.hist(mwhvzdipangles, bins=angbins,  alpha=0.3, color='green', label='MW Observer', normed=True)
plt.hist(hvzdipangles, bins=angbins,  alpha=0.3, color='blue', label='MW-2MRS Observer', normed=True)
plt.hist(dangles, bins =angbins, alpha=0.1, color='yellow', label='Diff between 2 random vectors', normed=True)
plt.ylabel('   ')
plt.xlabel('Angle between Halo Velocity and Dipole[degrees]', fontsize=15)
#plt.title('From MW-like Halos, 2MRS-like velocity configuration')
plt.legend(fontsize='15', loc='best')
plt.savefig('MWLike2MRSConfigDipAngles.png')
#plt.show()    


print 'Median Dipole', np.median(hvdipstrenghts)
print 'SD:', np.sqrt(np.var(hvdipstrenghts))
print 'Median Dipole z<0.03 suppressed', np.median(hvzdipstrength)
print 'SD:', np.sqrt(np.var(hvzdipstrength))


plt.figure(2)  
#plt.hist(hvdipstrenghts, bins=dipbins,  alpha=0.3, color='green', label='Full Distribution')

plt.hist(np.asarray(mwhvzdipstrength)-0.0017, bins=dipbins,  alpha=0.3, color='green', label='MW Observer', normed=True)
plt.hist(np.asarray(hvzdipstrength)-0.0008, bins=dipbins,  alpha=0.3, color='blue', label='MW-2MRS Observer', normed=True)
#plt.hist(dangles, bins =angbins, alpha=0.1, color='pink', label='Diff between 2 random vectors')
plt.xlabel('Size of observed dipole', fontsize=15)
#plt.title('From MW-like Halos, 2MRS-like velocity configuration')
plt.legend(fontsize='15', loc='best')
plt.savefig('MWLike2MRSConfigDipStrengths.png')
plt.show()         
        
        
        
    
