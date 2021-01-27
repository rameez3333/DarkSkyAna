import yt
import numpy as np
from astropy.cosmology import Planck15 as cosmo
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy import units as u
from astropy.coordinates import SkyCoord
from optparse import OptionParser
import gc

inv_comoving_distance = InterpolatedUnivariateSpline(cosmo.comoving_distance(np.linspace(0, 0.6, 1000)).value, np.linspace(0, 0.6, 1000))

usage = 'usage: %prog [options]'
parser = OptionParser(usage)

parser.add_option("-x", "--xcenter", action="store", type="float", default=3000.0, dest="XC", help="Xcoord of the Center")
parser.add_option("-y", "--ycenter", action="store", type="float", default=3000.0, dest="YC", help="Ycoord of the Center")
parser.add_option("-z", "--zcenter", action="store", type="float", default=3000.0, dest="ZC", help="Zcoord of the Center")
parser.add_option("-w", "--width", action="store", type="float", default=3500., dest="WIDTH", help="Width of the box, radius of the circle, in MPc")
parser.add_option("-s", "--selecttype", action = "store", type="int", default=4, dest="SELTP", help = "Type of observer. 0: Copernican, 1: MWlike in mass and velocity, 2. MWlike in mass and velocity, Virgo Cluster nearby. 3. MWlike in mass and velocity, Virgo cluster nearby, LG Flow-Virgo angle between 40 and 50 degrees  4. MWlike in mass and velocity, Virgo cluster nearby, bulk flow out to z=0.03 5. 6dfGSv. 6. SMAC1. 7. SMAC2. 8. SMAC3. 9. 2M++")
parser.add_option("-o", "--orient", action="store", type="string", default='3Flow:Virgo', dest="ORIENT", help="What reference directions to orient by, in order. Options are Virgo, LGFlow and 3Flow.")
parser.add_option("-d", "--diag", action = "store_true", default=False, dest="DIAG", help = "diagnostics")

(options, args) = parser.parse_args()

if options.SELTP==5:
    options.ORIENT='6Flow:Virgo'
if options.SELTP==6:
    options.ORIENT='SMAC1Flow:Virgo'
if options.SELTP==7:
    options.ORIENT='SMAC2Flow:Virgo'
if options.SELTP==8:
    options.ORIENT='SMAC3Flow:Virgo'
if options.SELTP==9:
    options.ORIENT='2MppFlow:Virgo'

SELC = options.SELTP

if options.DIAG:
    import os
    import psutil
    import sys
    process = psutil.Process(os.getpid())
    del os


bfscale = {5:240.27, 6:120./0.688062, 7:200./0.688062, 8:60./0.688062, 9:292.}
bfvels = {5: [329, np.inf], 6:[687.-203, 687.+203], 7:[372.-127. , 372.+127.], 8:[175., 225.], 9:[136., 182.]}
angrange = {5: [51.4,66.4], 6:[67.9,82.9], 7:[61.0, 76.0], 8:[58.0, 73.0], 9:[62.4, 77.4]}
angrange2 = {5: [20.4,35.4], 6:[26.2,41.2], 7:[16.7, 31.7], 8:[22.5, 37.5], 9:[28.1, 43.1]}
virglgangrange=[39.5, 49.5]
velLGrange=[622.-150., 622.+150.]

#bfvels = {5: [329, np.inf], 6:[687.-203, 687.+203], 7:[372.-127. , 372.+127.], 8:[175., 225.], 9:[40., np.inf]}
#angrange = {5: [51.4,66.4], 6:[67.9,82.9], 7:[61.0, 76.0], 8:[58.0, 73.0], 9:[10, 170]}
#angrange2 = {5: [20.4,35.4], 6:[26.2,41.2], 7:[16.7, 31.7], 8:[22.5, 37.5], 9:[10, 170]}
#virglgangrange=[10, 170]
#velLGrange=[10, 1000.]


c = 299792.458 # km/s
H0 = 68.8062 #(km/s) / Mpc


CMBdipra = 168.0
CMBdipdec = -7.0
ActVirgra = 186.63
ActVirgdec = 12.72
LGMra = 166.67
LGMdec = -27.33
SixDFGSvra = 193.0
SixDFGSvdec = -45.9
TwoMppra = 194.802
TwoMppdec = -56.855

SMAC1ra = 128.97 #120h^-1 bulk flow, 687+/-203
SMAC1dec = -40.7

SMAC2ra = 147.3 # 200h^-1 bulk flow, 372+/-127
SMAC2dec = -46.0

SMAC3ra = 188.1 # 60^-1 bulk flow, 225
SMAC3dec = -52.8

Virgc = SkyCoord(ra = ActVirgra*u.degree, dec=ActVirgdec*u.degree).cartesian
LGMc = SkyCoord(ra = LGMra*u.degree, dec=LGMdec*u.degree).cartesian
CMBdipc = SkyCoord(ra = CMBdipra*u.degree, dec=CMBdipdec*u.degree).cartesian
SixDFGSvc = SkyCoord(ra = SixDFGSvra*u.degree, dec=SixDFGSvdec*u.degree).cartesian
SMAC1c = SkyCoord(ra = SMAC1ra*u.degree, dec=SMAC1dec*u.degree).cartesian
SMAC2c = SkyCoord(ra = SMAC2ra*u.degree, dec=SMAC2dec*u.degree).cartesian
SMAC3c = SkyCoord(ra = SMAC3ra*u.degree, dec=SMAC3dec*u.degree).cartesian
TwoMppc = SkyCoord(ra = TwoMppra*u.degree, dec=TwoMppdec*u.degree).cartesian

#RealVdir = np.asarray([Virgc.x.value, Virgc.y.value, Virgc.z.value])
#RealLGMdir = np.asarray([LGMc.x.value, LGMc.y.value, LGMc.z.value])
#RealCMBDipdir = np.asarray([CMBdipc.x.value, CMBdipc.y.value, CMBdipc.z.value])

RealDirs={}
RealDirs['Virgo'] = np.asarray([Virgc.x.value, Virgc.y.value, Virgc.z.value])
RealDirs['LGFlow'] = np.asarray([LGMc.x.value, LGMc.y.value, LGMc.z.value])
RealDirs['3Flow'] = np.asarray([CMBdipc.x.value, CMBdipc.y.value, CMBdipc.z.value])
RealDirs['SMAC1Flow'] = np.asarray([SMAC1c.x.value, SMAC1c.y.value, SMAC1c.z.value])
RealDirs['SMAC2Flow'] = np.asarray([SMAC2c.x.value, SMAC2c.y.value, SMAC2c.z.value])
RealDirs['SMAC3Flow'] = np.asarray([SMAC3c.x.value, SMAC3c.y.value, SMAC3c.z.value])
RealDirs['2MppFlow'] = np.asarray([TwoMppc.x.value, TwoMppc.y.value, TwoMppc.z.value])

orients = options.ORIENT.split(':')


lgbfdangle = np.rad2deg(np.arccos(np.cos(np.deg2rad(ActVirgdec))*np.cos(np.deg2rad(LGMdec))*np.cos(np.deg2rad(ActVirgra) - np.deg2rad(LGMra))+np.sin(np.deg2rad(ActVirgdec))*np.sin(np.deg2rad(LGMdec))))


center = np.array([options.XC,  options.YC, options.ZC])*0.688062


width = options.WIDTH*0.688062


if options.SELTP>1:
    width = 3500.*0.688062

#if options.SELTP>2:
    #width = 4000.*0.688062

outfname = 'VelCovResultsBF/DarkSkyVelCovPanthQuery_CenterX_'+str(options.XC)+'_Y_'+str(options.YC)+'_Z_'+str(options.ZC)+'_W_'+str(width)



outfname = outfname+'_OBS_'+str(options.SELTP)

outfname = outfname+'_'+orients[0]+'_'+orients[1]

outfname=outfname+'_Out.txt'

bbox = np.array([center-(width+100.*0.688062)/2., center+(width+100.*0.688062)/2.])


ds = yt.load('Halos/ds14_a_halos_1.0000',midx_filename='Halos/ds14_a_halos_1.0000.midx9',bounding_box=bbox)
 
ad = ds.sphere(center, (width/2., 'Mpccm/h'))

if options.DIAG:
    print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)


X = np.asarray(ad['x'].in_units('Mpc'))
Y = np.asarray(ad['y'].in_units('Mpc'))
Z = np.asarray(ad['z'].in_units('Mpc'))

VX = np.asarray(ad['vx'].in_units('km/s'))
VY = np.asarray(ad['vy'].in_units('km/s'))
VZ = np.asarray(ad['vz'].in_units('km/s'))

M = np.asarray(ad['m200c'])

totno = len(M)

pick = np.log10(M) > 11.35

X=X[pick]
Y=Y[pick]
Z=Z[pick]
VX=VX[pick]
VY=VY[pick]
VZ=VZ[pick]
M=M[pick]

if options.DIAG:
    print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)

print 'Max, Min, halo masses', np.min(M), np.max(M)

xmed, ymed, zmed = np.median(X), np.median(Y), np.median(Z)

print 'Median:', xmed, ymed, zmed



del ad
del ds
gc.collect()



nobs=[]
hcount = []
if options.SELTP==0:
    disttomed = np.sqrt(np.power((X-xmed), 2) + np.power((Y-ymed), 2) + np.power((Z-zmed), 2))
    observerindex = np.argmin(disttomed)
    xob, yob, zob = X[observerindex], Y[observerindex], Z[observerindex]
    

if options.SELTP==1:
    VEL = np.sqrt(VX*VX+VY*VY+VZ*VZ)
    sel = np.log10(M)<12.15
    XS, YS, ZS = X[sel], Y[sel], Z[sel]
    mwmassno = len(XS)
    countobs=0
    while True:
        disttomed = np.sqrt(np.power((XS-xmed), 2) + np.power((YS-ymed), 2) + np.power((ZS-zmed), 2))
        observerindex = np.argmin(disttomed)
        xob, yob, zob = XS[observerindex], YS[observerindex], ZS[observerindex]
        print 'The observer is at', xob, yob, zob
        disttoob = np.sqrt(np.power((X-xob), 2) + np.power((Y-yob), 2) + np.power((Z-zob), 2))
        pick = disttoob < 3.125/0.688062
        VXSEL, VYSEL, VZSEL = VX[pick], VY[pick], VZ[pick]
        velLG = np.sqrt(np.power(np.sum(VXSEL), 2) + np.power(np.sum(VYSEL), 2) + np.power(np.sum(VZSEL), 2))/float(len(VXSEL))
        
        bulkvelx, bulkvely, bulkvelz = np.sum(VXSEL)/float(len(VXSEL)), np.sum(VYSEL)/float(len(VXSEL)), np.sum(VZSEL)/float(len(VXSEL))
        if (velLG > (622.-150.)) and (velLG < (622.+150.)):
            print 'Found LG like observer', velLG
            break
        XS = np.delete(XS, observerindex)
        YS = np.delete(YS, observerindex)
        ZS = np.delete(ZS, observerindex)
        countobs+=1

        
elif options.SELTP==2:
    vmsel = (M > 6.e14/0.688062) * (M < 1.8e15/0.688062)
    vhno = len(M[vmsel])
    print len(M[vmsel]), 'Virgolike halos found'
    XVl, YVl, ZVl = X[vmsel], Y[vmsel], Z[vmsel]
    vdisttomed = np.sqrt(np.power((XVl-xmed), 2) + np.power((YVl-ymed), 2) + np.power((ZVl-zmed), 2))
    vindex = np.argmin(vdisttomed)
    Xv, Yv, Zv = XVl[vindex], YVl[vindex], ZVl[vindex]
    disttov = np.sqrt(np.power((X-Xv), 2) + np.power((Y-Yv), 2) + np.power((Z-Zv), 2))
    #VEL = np.sqrt(VX*VX+VY*VY+VZ*VZ)
    sel = np.log10(M)<12.15*(disttov < 16./0.688062) * (disttov > 8./0.688062)
    XS, YS, ZS = X[sel], Y[sel], Z[sel]
    mwmassno = len(XS)
    countobs=0
    while True:
        print 'Trial no', countobs
        disttomed = np.sqrt(np.power((XS-xmed), 2) + np.power((YS-ymed), 2) + np.power((ZS-zmed), 2))
        observerindex = np.argmin(disttomed)
        xob, yob, zob = XS[observerindex], YS[observerindex], ZS[observerindex]
        print 'The observer is at', xob, yob, zob
        disttoob = np.sqrt(np.power((X-xob), 2) + np.power((Y-yob), 2) + np.power((Z-zob), 2))
        pick = disttoob < 3.125/0.688062
        VXSEL, VYSEL, VZSEL = VX[pick], VY[pick], VZ[pick]
        velLG = np.sqrt(np.power(np.sum(VXSEL), 2) + np.power(np.sum(VYSEL), 2) + np.power(np.sum(VZSEL), 2))/float(len(VXSEL))
        bulkvelx, bulkvely, bulkvelz = np.sum(VXSEL)/float(len(VXSEL)), np.sum(VYSEL)/float(len(VXSEL)), np.sum(VZSEL)/float(len(VXSEL))
        if (velLG > (622.-150.)) and (velLG < (622.+150.)):
            print 'Found LG like observer', velLG
            pick2 = (disttoob < 16./0.688062) * (disttoob > 8./0.688062) * (M > 6.e14/0.688062) * (M < 1.8e15/0.688062)
            if len(M[pick2])==1:
                print 'Found MW like observer in LG like environment with Virgolike cluster nearby. Cluster is at:'
                VirgX, VirgY, VirgZ = X[pick2], Y[pick2], Z[pick2]
                print VirgX, VirgY, VirgZ
                cosangLGBV = ((VirgX - xob)*bulkvelx + (VirgY-yob)*bulkvely + (VirgZ-zob)*bulkvelz)/(np.sqrt((VirgX - xob)*(VirgX - xob) + (VirgY - yob)*(VirgY - yob) + (VirgZ - zob)*(VirgZ - zob))*np.sqrt(bulkvelx*bulkvelx + bulkvely*bulkvely + bulkvelz*bulkvelz))
                print 'Angle between LG bulk flow and Virgolike cluster direction', np.rad2deg(np.arccos(cosangLGBV))
                break
        XS = np.delete(XS, observerindex)
        YS = np.delete(YS, observerindex)
        ZS = np.delete(ZS, observerindex)
        countobs+=1
                    
elif options.SELTP==3:
    vmsel = (M > 6.e14/0.688062) * (M < 1.8e15/0.688062)
    vhno = len(M[vmsel])
    print len(M[vmsel]), 'Virgolike halos found'
    XVl, YVl, ZVl = X[vmsel], Y[vmsel], Z[vmsel]
    vdisttomed = np.sqrt(np.power((XVl-xmed), 2) + np.power((YVl-ymed), 2) + np.power((ZVl-zmed), 2))
    vindex = np.argmin(vdisttomed)
    Xv, Yv, Zv = XVl[vindex], YVl[vindex], ZVl[vindex]
    disttov = np.sqrt(np.power((X-Xv), 2) + np.power((Y-Yv), 2) + np.power((Z-Zv), 2))
    #VEL = np.sqrt(VX*VX+VY*VY+VZ*VZ)
    sel = np.log10(M)<12.15*(disttov < 16./0.688062) * (disttov > 8./0.688062)
    XS, YS, ZS = X[sel], Y[sel], Z[sel]
    mwmassno = len(XS)
    countobs=0
    while True:
        print 'Trial no', countobs
        disttomed = np.sqrt(np.power((XS-xmed), 2) + np.power((YS-ymed), 2) + np.power((ZS-zmed), 2))
        observerindex = np.argmin(disttomed)
        xob, yob, zob = XS[observerindex], YS[observerindex], ZS[observerindex]
        print 'The observer is at', xob, yob, zob
        disttoob = np.sqrt(np.power((X-xob), 2) + np.power((Y-yob), 2) + np.power((Z-zob), 2))
        pick = disttoob < 3.125/0.688062
        VXSEL, VYSEL, VZSEL = VX[pick], VY[pick], VZ[pick]
        velLG = np.sqrt(np.power(np.sum(VXSEL), 2) + np.power(np.sum(VYSEL), 2) + np.power(np.sum(VZSEL), 2))/float(len(VXSEL))
        bulkvelx, bulkvely, bulkvelz = np.sum(VXSEL)/float(len(VXSEL)), np.sum(VYSEL)/float(len(VXSEL)), np.sum(VZSEL)/float(len(VXSEL))
        if (velLG > (622.-150.)) and (velLG < (622.+150.)):
            print 'Found LG like observer', velLG
            pick2 = (disttoob < 16./0.688062) * (disttoob > 8./0.688062) * (M > 6.e14/0.688062) * (M < 1.8e15/0.688062)
            if len(M[pick2])==1:
                print 'Found MW like observer in LG like environment with Virgolike cluster nearby. Cluster is at:'
                VirgX, VirgY, VirgZ = X[pick2], Y[pick2], Z[pick2]
                print VirgX, VirgY, VirgZ
                cosangLGBV = ((VirgX - xob)*bulkvelx + (VirgY-yob)*bulkvely + (VirgZ-zob)*bulkvelz)/(np.sqrt((VirgX - xob)*(VirgX - xob) + (VirgY - yob)*(VirgY - yob) + (VirgZ - zob)*(VirgZ - zob))*np.sqrt(bulkvelx*bulkvelx + bulkvely*bulkvely + bulkvelz*bulkvelz))
                print 'Angle between LG bulk flow and Virgolike cluster direction', np.rad2deg(np.arccos(cosangLGBV))
                if (np.rad2deg(np.arccos(cosangLGBV)) > 39.5) and (np.rad2deg(np.arccos(cosangLGBV)) < 49.5):
                    break
        XS = np.delete(XS, observerindex)
        YS = np.delete(YS, observerindex)
        ZS = np.delete(ZS, observerindex)
        countobs+=1    
    
    #(VEL>(622.-150.))*(VEL<(622. + 150.))


elif options.SELTP==4:
    notfound=True
    vmsel = (M > 6.e14/0.688062) * (M < 1.8e15/0.688062)
    print len(M[vmsel]), 'Virgolike halos found'
    vhno = len(M[vmsel])
    vindcount=0
    while notfound:
        vdisttomed = np.sqrt(np.power((X[vmsel]-xmed), 2) + np.power((Y[vmsel]-ymed), 2) + np.power((Z[vmsel]-zmed), 2))
        vindex = np.argsort(vdisttomed)[vindcount]
        Xv, Yv, Zv = X[vmsel][vindex], Y[vmsel][vindex], Z[vmsel][vindex]
        disttov = np.sqrt(np.power((X-Xv), 2) + np.power((Y-Yv), 2) + np.power((Z-Zv), 2))
        #VEL = np.sqrt(VX*VX+VY*VY+VZ*VZ)
        sel = np.log10(M)<12.15*(disttov < 16./0.688062) * (disttov > 8./0.688062)
        #XS, YS, ZS = X[sel], Y[sel], Z[sel]
        mwmassno = len(X[sel])
        print mwmassno, 'MW mass halos found in shell'
        countobs=0
        while True:
            print 'Trial no', countobs
            disttomed = np.sqrt(np.power((X[sel]-xmed), 2) + np.power((Y[sel]-ymed), 2) + np.power((Z[sel]-zmed), 2))
            try:
                observerindex = np.argsort(disttomed)[countobs]
            except:
                print 'Uh oh, it apepars we have run out of MW like halos around this supercluster.'
                nobs.append(countobs)
                hcount.append(mwmassno)
                break
            xob, yob, zob = X[sel][observerindex], Y[sel][observerindex], Z[sel][observerindex]
            print 'The observer is at', xob, yob, zob
            disttoob = np.sqrt(np.power((X-xob), 2) + np.power((Y-yob), 2) + np.power((Z-zob), 2))
            pick = disttoob < 3.125/0.688062
            #VXSEL, VYSEL, VZSEL = VX[pick], VY[pick], VZ[pick]
            velLG = np.sqrt(np.power(np.sum(VX[pick]), 2) + np.power(np.sum(VY[pick]), 2) + np.power(np.sum(VZ[pick]), 2))/float(len(VX[pick]))
            bulkvelx, bulkvely, bulkvelz = np.sum(VX[pick])/float(len(VX[pick])), np.sum(VY[pick])/float(len(VX[pick])), np.sum(VZ[pick])/float(len(VX[pick]))
            if (velLG > (622.-150.)) and (velLG < (622.+150.)):
                print 'Found LG like observer', velLG
                pick2 = (disttoob < 16./0.688062) * (disttoob > 8./0.688062) * (M > 6.e14/0.688062) * (M < 1.8e15/0.688062)
                if len(M[pick2])==1:
                    print 'Found MW like observer in LG like environment with Virgolike cluster nearby. Cluster is at:'
                    VirgX, VirgY, VirgZ = X[pick2][0], Y[pick2][0], Z[pick2][0]
                    print VirgX, VirgY, VirgZ
                    cosangLGBV = ((VirgX - xob)*bulkvelx + (VirgY-yob)*bulkvely + (VirgZ-zob)*bulkvelz)/(np.sqrt((VirgX - xob)*(VirgX - xob) + (VirgY - yob)*(VirgY - yob) + (VirgZ - zob)*(VirgZ - zob))*np.sqrt(bulkvelx*bulkvelx + bulkvely*bulkvely + bulkvelz*bulkvelz))
                    print 'Angle between LG bulk flow and Virgolike cluster direction', np.rad2deg(np.arccos(cosangLGBV))
                    if (np.rad2deg(np.arccos(cosangLGBV)) > 39.5) and (np.rad2deg(np.arccos(cosangLGBV)) < 49.5):
                        pick3 = (disttoob < 131.84)
                        #VXSEL3, VYSEL3, VZSEL3 = VX[pick3], VY[pick3], VZ[pick3]
                        vel03 = np.sqrt(np.power(np.sum(VX[pick3]), 2) + np.power(np.sum(VY[pick3]), 2) + np.power(np.sum(VZ[pick3]), 2))/float(len(VX[pick3]))
                        bulkvelx3, bulkvely3, bulkvelz3 = np.sum(VX[pick3])/float(len(VX[pick3])), np.sum(VY[pick3])/float(len(VX[pick3])), np.sum(VZ[pick3])/float(len(VX[pick3]))
                        print 'Bulk velocity out to z=0.03', vel03
                        if vel03>240.:
                            notfound=False
                            break
            #XS = np.delete(XS, observerindex)
            #YS = np.delete(YS, observerindex)
            #ZS = np.delete(ZS, observerindex)
            countobs+=1
        if notfound:
            print 'Moving on to shell around the next Virgolike Halo'
            vindcount = vindcount+1



elif options.SELTP>=5:
    notfound=True
    vmsel = (M > 6.e14/0.688062) * (M < 1.8e15/0.688062)
    print len(M[vmsel]), 'Virgolike halos found'
    vhno = len(M[vmsel])
    #XVl, YVl, ZVl = X[vmsel], Y[vmsel], Z[vmsel]
    vindcount=0
    while notfound:
        vdisttomed = np.sqrt(np.power((X[vmsel]-xmed), 2) + np.power((Y[vmsel]-ymed), 2) + np.power((Z[vmsel]-zmed), 2))
        vindex = np.argsort(vdisttomed)[vindcount]
        Xv, Yv, Zv = X[vmsel][vindex], Y[vmsel][vindex], Z[vmsel][vindex]
        disttov = np.sqrt(np.power((X-Xv), 2) + np.power((Y-Yv), 2) + np.power((Z-Zv), 2))
        #VEL = np.sqrt(VX*VX+VY*VY+VZ*VZ)
        sel = np.log10(M)<12.15*(disttov < 16./0.688062) * (disttov > 8./0.688062)
        #XS, YS, ZS = X[sel], Y[sel], Z[sel]
        mwmassno = len(X[sel])
        print mwmassno, 'MW mass halos found in shell'
        countobs=0
        if options.DIAG:
            print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)
        while True:
            print 'Trial no', countobs
            disttomed = np.sqrt(np.power((X[sel]-xmed), 2) + np.power((Y[sel]-ymed), 2) + np.power((Z[sel]-zmed), 2))
            try:
                observerindex = np.argsort(disttomed)[countobs]
            except:
                print 'Uh oh, it apepars we have run out of MW like halos around this supercluster.'
                nobs.append(countobs)
                hcount.append(mwmassno)
                break
            xob, yob, zob = X[sel][observerindex], Y[sel][observerindex], Z[sel][observerindex]
            print 'The observer is at', xob, yob, zob
            disttoob = np.sqrt(np.power((X-xob), 2) + np.power((Y-yob), 2) + np.power((Z-zob), 2))
            if options.DIAG:
                print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)
            pick = disttoob < 3.125/0.688062
            #VXSEL, VYSEL, VZSEL = VX[pick], VY[pick], VZ[pick]
            #velLG = np.sqrt(np.power(np.sum(VX[pick]), 2) + np.power(np.sum(VY[pick]), 2) + np.power(np.sum(VZ[pick]), 2))/float(len(VX[pick]))
            bulkvelx, bulkvely, bulkvelz = np.sum(VX[pick])/float(len(VX[pick])), np.sum(VY[pick])/float(len(VX[pick])), np.sum(VZ[pick])/float(len(VX[pick]))
            velLG = np.sqrt((bulkvelx**2. + bulkvely**2. + bulkvelz**2.))
            if (velLG > velLGrange[0]) and (velLG < velLGrange[1]):
                print 'Found LG like observer', velLG
                pick2 = (disttoob < 16./0.688062) * (disttoob > 8./0.688062) * (M > 6.e14/0.688062) * (M < 1.8e15/0.688062)
                if len(M[pick2])==1:
                    print 'Found MW like observer in LG like environment with Virgolike cluster nearby. Cluster is at:'
                    VirgX, VirgY, VirgZ = X[pick2][0], Y[pick2][0], Z[pick2][0]
                    print VirgX, VirgY, VirgZ
                    if options.DIAG:
                        print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)
                    cosangLGBV = ((VirgX - xob)*bulkvelx + (VirgY-yob)*bulkvely + (VirgZ-zob)*bulkvelz)/(np.sqrt((VirgX - xob)*(VirgX - xob) + (VirgY - yob)*(VirgY - yob) + (VirgZ - zob)*(VirgZ - zob))*np.sqrt(bulkvelx*bulkvelx + bulkvely*bulkvely + bulkvelz*bulkvelz))
                    print 'Angle between LG bulk flow and Virgolike cluster direction', np.rad2deg(np.arccos(cosangLGBV))
                    if (np.rad2deg(np.arccos(cosangLGBV)) > virglgangrange[0]) and (np.rad2deg(np.arccos(cosangLGBV)) < virglgangrange[1]):
                        pick6 = (disttoob < bfscale[SELC])
                        #VXSEL6, VYSEL6, VZSEL6 = VX[pick6], VY[pick6], VZ[pick6]
                        #vel06 = np.sqrt(np.power(np.sum(VX[pick6]), 2) + np.power(np.sum(VY[pick6]), 2) + np.power(np.sum(VZ[pick6]), 2))/float(len(VX[pick6]))
                        bulkvelx6, bulkvely6, bulkvelz6 = np.sum(VX[pick6])/float(len(VX[pick6])), np.sum(VY[pick6])/float(len(VX[pick6])), np.sum(VZ[pick6])/float(len(VX[pick6]))
                        vel06 = np.sqrt((bulkvelx6**2. + bulkvely6**2. + bulkvelz6**2.))
                        print 'Bulk velocity out to ',bfscale[SELC], ':', vel06
                        if options.DIAG:
                            print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)
                        if ((vel06>bfvels[SELC][0]) and (vel06<bfvels[SELC][1])):
                            cosangVirgBV6 = ((VirgX - xob)*bulkvelx6 + (VirgY-yob)*bulkvely6 + (VirgZ-zob)*bulkvelz6)/(np.sqrt((VirgX - xob)*(VirgX - xob) + (VirgY - yob)*(VirgY - yob) + (VirgZ - zob)*(VirgZ - zob))*np.sqrt(bulkvelx6*bulkvelx6 + bulkvely6*bulkvely6 + bulkvelz6*bulkvelz6))
                            print 'Angle between the bulk flow and Virgolike cluster direction', np.rad2deg(np.arccos(cosangVirgBV6))
                            if options.DIAG:
                                print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)
                            if (np.rad2deg(np.arccos(cosangVirgBV6)) > angrange[SELC][0]) and (np.rad2deg(np.arccos(cosangVirgBV6)) < angrange[SELC][1]):
                                cosangLGBV6 = (bulkvelx*bulkvelx6 + bulkvely*bulkvely6 + bulkvelz*bulkvelz6)/(np.sqrt(bulkvelx**2 + bulkvely**2 + bulkvelz**2)*np.sqrt(bulkvelx6*bulkvelx6 + bulkvely6*bulkvely6 + bulkvelz6*bulkvelz6))
                                if (np.rad2deg(np.arccos(cosangLGBV6)) > angrange2[SELC][0]) and (np.rad2deg(np.arccos(cosangLGBV6)) < angrange2[SELC][1]):
                                    if options.DIAG:
                                        print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)
                                    notfound=False
                                    break
            #XS = np.delete(XS, observerindex)
            #YS = np.delete(YS, observerindex)
            #ZS = np.delete(ZS, observerindex)
            countobs+=1
        if notfound:
            print 'Moving on to shell around the next Virgolike Halo'
            if options.DIAG:
                print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)
            #XVl = np.delete(XVl, vindex)
            #YVl = np.delete(YVl, vindex)
            #ZVl = np.delete(ZVl, vindex)
            vindcount = vindcount+1











print 'Nearest suitable Galaxy, distance from median : ', disttomed[observerindex]
DMEDSUIT = disttomed[observerindex]

if options.DIAG:
    print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)


fout = open(outfname, 'w')
fout.write(str(DMEDSUIT)+'\n')
if options.SELTP==1:
    fout.write(str(totno)+','+str(mwmassno)+','+str(countobs)+'\n')
    fout.write(str(velLG)+'\n')

if options.SELTP==2:
    fout.write(str(totno)+','+str(vhno)+','+str(mwmassno)+','+str(countobs)+'\n')
    fout.write(str(velLG)+'\n')
    fout.write(str(np.rad2deg(np.arccos(cosangLGBV)))+'\n')

if options.SELTP==3:
    fout.write(str(totno)+','+str(vhno)+','+str(mwmassno)+','+str(countobs)+'\n')
    fout.write(str(velLG)+'\n')
    fout.write(str(np.rad2deg(np.arccos(cosangLGBV)))+'\n')

if options.SELTP==4:
    fout.write(str(totno)+','+str(vhno)+'\n')
    for w in nobs:
        fout.write(str(w)+',')
    fout.write('\n')
    for w in hcount:
        fout.write(str(w)+',')
    fout.write('\n')    
    fout.write(str(velLG)+','+str(vel03)+'\n')
    fout.write(str(np.rad2deg(np.arccos(cosangLGBV)))+'\n')

if options.SELTP>=5:
    fout.write(str(totno)+','+str(vhno)+'\n')
    for w in nobs:
        fout.write(str(w)+',')
    fout.write('\n')
    for w in hcount:
        fout.write(str(w)+',')
    fout.write('\n')    
    fout.write(str(velLG)+','+str(vel06)+'\n')
    fout.write(str(np.rad2deg(np.arccos(cosangVirgBV6)))+','+str(np.rad2deg(np.arccos(cosangLGBV6)))+'\n')



#Import JLA data
#cols are z,m,x,c,cluster mass, survey
#K = np.load('../Dipole_JLA/SNMLE/JLADirZInc.npy');
#Same as JLA.npy, but with additional columns [6] is RAdeg, [7] is DECdeg, [8] is Zcmb, [9] is Zhel
#Jeppe used Zcmb, rounded ([0])

#import pantheon now
jlarr = np.genfromtxt('../Dipole_JLA/Pantheon/Pantheon_G10.txt')
snnames = np.asarray([k.split()[1] for k in open('../Dipole_JLA/Pantheon/Pantheon_G10.txt').readlines()])



zcmb = jlarr.transpose()[7]
ras, decs = jlarr.transpose()[34], jlarr.transpose()[35]


ras, decs = ras[zcmb<0.5], decs[zcmb<0.5]
snnames = snnames[zcmb<0.5]

zcmb = zcmb[zcmb<0.5]


SNc = SkyCoord(ra = ras*u.degree, dec=decs*u.degree, distance=cosmo.comoving_distance(zcmb)).cartesian

del jlarr

gc.collect

SNx = SNc.x.value
SNy = SNc.y.value
SNz = SNc.z.value

if options.DIAG:
    print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)


MCDirs={}
if options.SELTP>=1:
    MCLGMdir = np.asarray([bulkvelx, bulkvely, bulkvelz])
    MCLGMdir = MCLGMdir/np.sqrt(np.sum(MCLGMdir*MCLGMdir))
    MCDirs['LGFlow'] = MCLGMdir

if options.SELTP>=2:
    MCVdir = np.asarray([VirgX-xob, VirgY-yob, VirgZ-zob])
    MCVdir = MCVdir/np.sqrt(np.sum(MCVdir*MCVdir))
    MCDirs['Virgo'] = MCVdir

if options.SELTP==4:
    MC3dir = np.asarray([bulkvelx3, bulkvely3, bulkvelz3])
    MC3dir = MC3dir/np.sqrt(np.sum(MC3dir*MC3dir))
    MCDirs['3Flow'] = MC3dir

if options.SELTP>=5:
    MC6dir = np.asarray([bulkvelx6, bulkvely6, bulkvelz6])
    MC6dir = MC6dir/np.sqrt(np.sum(MC6dir*MC6dir))
    MCDirs[orients[0]] = MC6dir


if options.DIAG:
    print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)

Zdir = np.asarray([0., 0., 1.])


def RotateAtoB(A, B):
    c = np.sum(A*B)
    s = np.sqrt(1.-c*c)
    v = np.cross(A, B)
    sqv = np.matrix([[0, -1.*v[2], v[1]], [v[2], 0, -1.*v[0]], [-1.*v[1], v[0], 0.]])
    R = np.identity(3) + sqv + sqv*sqv*(1.-c)/(s*s)
    return R


#remember that the Z axis is 0.0, 90.0 in right ascension, declination
def CartToDecRa(X):
    skc = SkyCoord(x=X[0], y=X[1], z=X[2], representation='cartesian').frame.represent_as('spherical')
    return skc.lon.degree, skc.lat.degree
    

if options.SELTP==1:
    R = RotateAtoB(RealDirs['LGFlow'], MCDirs['LGFlow'])
    Rinv = R**(-1)
    
if options.SELTP>=2:
    R = RotateAtoB(RealDirs[orients[1]], MCDirs[orients[1]])
    D = np.asarray(np.dot(R, RealDirs[orients[0]])).tolist()[0]
    R2 = RotateAtoB(D, MCDirs[orients[0]])
    R=R2*R
    Rinv = R**(-1)



for dkey in MCDirs.keys():
    print dkey
    r,d = CartToDecRa(np.asarray(np.dot(Rinv, MCDirs[dkey]).tolist()[0]))
    print r, d
    fout.write(dkey+','+str(r)+','+str(d)+'\n')


#Nlgmra, Nlgmdec = CartToDecRa(np.asarray(Nlgm)[0])

if options.SELTP>=1:
    SNRotCoords = np.asarray([np.asarray(np.dot(R, [p,q,r]).tolist()[0]) for p, q, r in zip(SNx, SNy, SNz)])
else:
    SNRotCoords = np.vstack([SNx, SNy, SNz]).transpose()

if options.DIAG:
    print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)

X, Y, Z = X-xob, Y-yob, Z-zob

if options.DIAG:
    print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)

#HX, HY, HZ, HXV, HYV, HZV, HSnD = [], [], [], [], [], [], []

for snc, jl, sname in zip(SNRotCoords, zcmb, snnames):
    disttoSN = np.power((X-snc[0]), 2) + np.power((Y-snc[1]), 2) + np.power((Z-snc[2]), 2)
    clssthindex = np.argmin(disttoSN)
    #HX.append(X[clssthindex])
    #HY.append(Y[clssthindex])
    #HZ.append(Z[clssthindex])
    #HXV.append(VX[clssthindex])
    #HYV.append(VY[clssthindex])
    #HZV.append(VZ[clssthindex])
    #HSnD.append(np.sqrt(disttoSN[clssthindex]))
    print 'SNcoords, sntype, zhel', snc[0], snc[1], snc[2], '::::', jl, sname
    print 'HaloCoords', X[clssthindex], Y[clssthindex], Z[clssthindex]
    print 'HaloVels', VX[clssthindex], VY[clssthindex], VZ[clssthindex]
    print 'HaloSep', np.sqrt(disttoSN[clssthindex])
    fout.write(str(jl)+','+str(sname)+','+str(X[clssthindex])+','+str(Y[clssthindex])+','+str(Z[clssthindex])+','+str(VX[clssthindex])+','+str(VY[clssthindex])+','+str(VZ[clssthindex])+','+str(np.sqrt(disttoSN[clssthindex]))+'\n')
    if options.DIAG:
         print 'Line: ', sys._getframe().f_lineno, 'Mem : ', process.memory_info().rss/(1024.**2.)
    #fout.write(str(jl)+','+str(sname)+','+str(HX[-1])+','+str(HY[-1])+','+str(HZ[-1])+','+str(HXV[-1])+','+str(HYV[-1])+','+str(HZV[-1])+','+str(HSnD[-1])+'\n')
    


fout.close()






















#Rm = RotateAtoB(MCVdir, Zdir)

#MClgm = np.dot(Rm, MCLGMdr)
#MClgmra, MClgmdec = CartToDecRa(np.asarray(MClgm)[0])

#MCraoffset = MClgmra - Nlgmra

#mccoords = np.vstack([X, Y, Z])
#mcvels = np.vstack([VX, VY, VZ])
#del X, Y, Z, VX, VY, VZ

#mcRotcoords = np.asarray([np.asarray(np.dot(Rm, [p,q,r]).tolist()[0]) for p, q, r in zip(mccoords[0], mccoords[1], mccoords[2])])

#mcras, mcdecs = CartToDecRa(mcRotcoords)
#mcras = mcras+MCraoffset

#mcRotvels = np.asarray([np.asarray(np.dot(Rm, [p,q,r]).tolist()[0]) for p, q, r in zip(mcvels[0], mcvels[1], mcvels[2])])






