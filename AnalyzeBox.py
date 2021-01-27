import yt
import numpy as np
from astropy.cosmology import Planck15 as cosmo
from scipy.interpolate import InterpolatedUnivariateSpline
from optparse import OptionParser
import gc

inv_comoving_distance = InterpolatedUnivariateSpline(cosmo.comoving_distance(np.linspace(0, 0.6, 1000)).value, np.linspace(0, 0.6, 1000))

usage = 'usage: %prog [options]'
parser = OptionParser(usage)

parser.add_option("-x", "--xcenter", action="store", type="float", default=3000.0, dest="XC", help="Xcoord of the Center")
parser.add_option("-y", "--ycenter", action="store", type="float", default=3000.0, dest="YC", help="Ycoord of the Center")
parser.add_option("-z", "--zcenter", action="store", type="float", default=3000.0, dest="ZC", help="Zcoord of the Center")
parser.add_option("-w", "--width", action="store", type="float", default=4200., dest="WIDTH", help="Width of the box, radius of the circle, in MPc")
parser.add_option("-s", "--selectmw", action = "store_true", default=False, dest="SELMW", help = "Select Milkyway like halos in mass and velocity")
parser.add_option("-2", "--select2mrs", action = "store_true", default=False, dest="SEL2M", help = "Select Milkyway like halos in mass and velocity, within a 2MRS like velocity config as determined by roya")

(options, args) = parser.parse_args()



center = np.array([options.XC,  options.YC, options.ZC])*0.688062

width = options.WIDTH*0.688062

if options.SELMW:
    width = 4400.*0.688062

if options.SEL2M:
    wifth = 4600.*0.688062

outfname = 'DarkSkyDipoleQuery_CenterX_'+str(options.XC)+'_Y_'+str(options.YC)+'_Z_'+str(options.ZC)+'_W_'+str(width)



if options.SELMW:
    outfname=outfname+'MWSELECT'

if options.SEL2M:
    outfname=outfname+'MWSELECT2MRS0.03'

outfname=outfname+'_Out.txt'

bbox = np.array([center-(width+100.*0.688062)/2., center+(width+100.*0.688062)/2.])


ds = yt.load('Halos/ds14_a_halos_1.0000',midx_filename='Halos/ds14_a_halos_1.0000.midx9',bounding_box=bbox)
 
ad = ds.sphere(center, (width/2., 'Mpccm/h'))


X = np.asarray(ad['x'].in_units('Mpc'))
Y = np.asarray(ad['y'].in_units('Mpc'))
Z = np.asarray(ad['z'].in_units('Mpc'))

VX = np.asarray(ad['vx'].in_units('km/s'))
VY = np.asarray(ad['vy'].in_units('km/s'))
VZ = np.asarray(ad['vz'].in_units('km/s'))

M = np.asarray(ad['m200c'])

pick = np.log10(M) > 11.35

X=X[pick]
Y=Y[pick]
Z=Z[pick]
VX=VX[pick]
VY=VY[pick]
VZ=VZ[pick]
M=M[pick]

xmed, ymed, zmed = np.median(X), np.median(Y), np.median(Z)

print 'Median:', xmed, ymed, zmed


if options.SELMW:
    VEL = np.sqrt(VX*VX+VY*VY+VZ*VZ)
    sel = (VEL>600.)*(VEL<650.)*(np.log10(M)<12.15)
    XS, YS, ZS = X[sel], Y[sel], Z[sel]
    disttomed = np.power((XS-xmed), 2) + np.power((YS-ymed), 2) + np.power((ZS-zmed), 2)
    observerindex = np.argmin(disttomed)

if options.SEL2M:
    VEL = np.sqrt(VX*VX+VY*VY+VZ*VZ)
    sel = (VEL>600.)*(VEL<650.)*(np.log10(M)<12.15)
    XS, YS, ZS = X[sel], Y[sel], Z[sel]
    
    while True:
        disttomed = np.sqrt(np.power((XS-xmed), 2) + np.power((YS-ymed), 2) + np.power((ZS-zmed), 2))
        observerindex = np.argmin(disttomed)
        xob, yob, zob = XS[observerindex], YS[observerindex], ZS[observerindex]
        print 'The observer is at', xob, yob, zob
        disttoob = np.sqrt(np.power((X-xob), 2) + np.power((Y-yob), 2) + np.power((Z-zob), 2))
        pick = disttoob < 131.84
        VXSEL, VYSEL, VZSEL = VX[pick], VY[pick], VZ[pick]
        vel03 = np.sqrt(np.power(np.sum(VXSEL), 2) + np.power(np.sum(VYSEL), 2) + np.power(np.sum(VZSEL), 2))/float(len(VXSEL))
        
        bulkvelx, bulkvely, bulkvelz = np.sum(VXSEL)/float(len(VXSEL)), np.sum(VYSEL)/float(len(VXSEL)), np.sum(VZSEL)/float(len(VXSEL))
        
        print 'Average velocity within redshift 0.03:', vel03
        if (vel03<300.) and (vel03>240.):
            break
        XS = np.delete(XS, observerindex)
        YS = np.delete(YS, observerindex)
        ZS = np.delete(ZS, observerindex)
        
    
    

else:
    disttomed = np.power((X-xmed), 2) + np.power((Y-ymed), 2) + np.power((Z-zmed), 2)
    observerindex = np.argmin(disttomed)

print 'Nearest suitable Galaxy, distance from median : ', np.sqrt(disttomed[observerindex])

if options.SELMW:
    xob, yob, zob = XS[observerindex], YS[observerindex], ZS[observerindex]
elif not options.SEL2M:
    xob, yob, zob = X[observerindex], Y[observerindex], Z[observerindex]
#del disttomed, X, Y, Z
#gc.collect()

#center = np.array([xob,  yob, zob])

#ad = ds.sphere(center, (width/2., 'Mpccm/h'))

#X = np.asarray(ad['x'].in_units('Mpc'))
#Y = np.asarray(ad['y'].in_units('Mpc'))
#Z = np.asarray(ad['z'].in_units('Mpc'))

#VX = np.asarray(ad['vx'].in_units('km/s'))
#VY = np.asarray(ad['vy'].in_units('km/s'))
#VZ = np.asarray(ad['vz'].in_units('km/s'))

#M = np.asarray(ad['m200c'])

print 'Found ', len(X), ' halos'

disttoobsq = np.power((X-xob), 2) + np.power((Y-yob), 2) + np.power((Z-zob), 2)

print disttoobsq

observerindex = np.argmin(disttoobsq)

fluxatob = M/disttoobsq

nhalos = len(X)

disttoob = np.sqrt(disttoobsq)

print 'Maximum Distance to Observer', np.nanmax(disttoob)



redshift = inv_comoving_distance(disttoob)

halarr = np.vstack([X-xob,Y-yob,Z-zob, VX,VY,VZ, M, fluxatob, redshift, disttoob])
del X,Y,Z, VX,VY,VZ, M, fluxatob, redshift, disttoob, ad, ds, disttoobsq
gc.collect()
observer = halarr.transpose()[observerindex]

fout = open(outfname, 'w')

print 'Observer', observer
fout.write('Number of Halos: '+str(nhalos)+'\n')
fout.write('Observer: '+str(xob)+'|'+str(yob)+'|'+str(zob)+'|'+str(observer[3])+'|'+str(observer[4])+'|'+str(observer[5])+'|'+str(observer[6])+'\n' )

if options.SEL2M:
    fout.write('Observer Bulk z to 0.03:' + str(bulkvelx)+'|'+str(bulkvely)+'|'+str(bulkvelz)+'\n' )

fout.write('Flux Limited Trials:\n')

halarr = np.delete(halarr.transpose(),observerindex, axis=0).transpose()

halarr = halarr.transpose()[np.log10(halarr[6])>11.6].transpose()
halarr = halarr.transpose()[halarr[9]<2000.].transpose()

print 'Minimum Distance to Observer', np.nanmin(halarr[9])


#find center


fluxthreshes = np.power(np.linspace(np.log10(np.min(halarr[7])), np.log10(np.max(halarr[7])), 5000, ), 10.)

zhistbins = np.linspace(0, 0.6, 61)


for ft in fluxthreshes:
    halarr = halarr.transpose()[halarr[7]>ft].transpose()
    ncuthalos=len(halarr[0])
    cx = np.sum(halarr[0]/halarr[9])
    cy = np.sum(halarr[1]/halarr[9])
    cz = np.sum(halarr[2]/halarr[9])
    
    zchalarr = halarr.transpose()[halarr[8]>0.03].transpose()
    zncuthalos=len(zchalarr[0])
    
    zcx = np.sum(zchalarr[0]/zchalarr[9])
    zcy = np.sum(zchalarr[1]/zchalarr[9])
    zcz = np.sum(zchalarr[2]/zchalarr[9])
    medz = np.median(halarr[8])
    zhist = np.histogram(halarr[8], bins=zhistbins)[0]
    strhist=''
    for b in zhist:
        strhist = strhist+str(b)+','
    fout.write(str(ft)+'|'+str(ncuthalos)+'|'+str(cx)+'|'+str(cy)+'|'+str(cz)+'|'+str(zncuthalos)+'|'+str(zcx)+'|'+str(zcy)+'|'+str(zcz)+'|'+str(medz)+'|'+strhist[:len(strhist)-1]+'\n')
    
fout.close()








