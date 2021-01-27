import yt
import numpy as np
from astropy.cosmology import Planck15 as cosmo
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy import units as u
from astropy.coordinates import SkyCoord
from optparse import OptionParser
from glob import glob
import gc
import matplotlib.pyplot as plt

usage = 'usage: %prog [options]'
parser = OptionParser(usage)

#parser.add_option("-i", "--index", action="store", type="int", default=0, dest="IN", help="Index of the file")
parser.add_option("-f", "--fname", action="store", type="string", default='', dest="FN", help="Filename to extract central coordinates from")
(options, args) = parser.parse_args()


centdencontr=[]
rads=[]


def FetchRegion(index):
    global centdencontr
    global rads
    print 'This index', index
    fvoids = open('Overlaps/OverLapSummary'+str(index)+'.txt').readlines()
    if options.FN=='':
        print 'No catalogue filename found, inferring from index'
        flist = sorted(glob('ForZobovDISP2wtf/*_Out.txt'))
        fname = flist[index]
    else:
        fname = options.FN

    XC = float(fname[fname.find('X_')+2:fname.find('_Y')])
    YC = float(fname[fname.find('Y_')+2:fname.find('_Z')])
    ZC = float(fname[fname.find('Z_')+2:fname.find('_W')])

    LIGHT_SPEED = 299792.458
    center = np.array([XC,  YC, ZC])*0.688062
    width = 1300.*0.688062


    vlistdict={}
    xlistdict={}
    ylistdict={}
    zlistdict={}
    radlistdict={}

    for ind in range(0, len(fvoids)):
        if len(fvoids[ind].split(' ')) == 1:
            vname = int(fvoids[ind].split(' ')[0])
            vlistdict[vname]=[]
            xlistdict[vname]=[]
            ylistdict[vname]=[]
            zlistdict[vname]=[]
            radlistdict[vname]=[]
            for void in fvoids[ind+1:]:
                if len(void.split(' '))==1:
                    break
                else:
                    vlistdict[vname].append(int(void.split(' ')[0]))
                    xlistdict[vname].append(float(void.split(' ')[1]))
                    ylistdict[vname].append(float(void.split(' ')[2]))
                    zlistdict[vname].append(float(void.split(' ')[3]))
                    try:
                        radlistdict[vname].append((3.*float(void.split(' ')[10])/(4.*np.pi))**(1./3.))
                    except:
                        radlistdict[vname].append(0)
                    
    for vname in vlistdict.keys():
        vlistdict[vname] = np.asarray(vlistdict[vname])
        xlistdict[vname] = np.asarray(xlistdict[vname])
        ylistdict[vname] = np.asarray(ylistdict[vname])
        zlistdict[vname] = np.asarray(zlistdict[vname])
        radlistdict[vname] = np.asarray(radlistdict[vname])
        
    xmed, ymed, zmed={},{},{}

    for vname in vlistdict.keys():
        if len(vlistdict[vname]) > 1:
            xmed[vname] =  np.mean(xlistdict[vname])
            ymed[vname] =  np.mean(ylistdict[vname])
            zmed[vname] =  np.mean(zlistdict[vname])
    
    distdict={}
    #desired=False
    for vname in xmed.keys():
        dist = np.sqrt(np.power(xlistdict[vname] - xmed[vname], 2.) + np.power(ylistdict[vname] - ymed[vname], 2.) + np.power(zlistdict[vname] - zmed[vname], 2.))
        print dist
        distdict[vname] = dist
         
        #if ((len(distdict[vname]) > 2.) and (np.max(distdict[vname])<400.)):
            #print np.max(radlistdict[vname]), 'Maximum Radius'
            #maxrad = np.max(radlistdict[vname])
            #desired=True
    
    #if not desired:
        #return 0

    
    bbox = np.array([center-(width+100.*0.688062)/2., center+(width+100.*0.688062)/2.])
    ds = yt.load('Halos/ds14_a_halos_1.0000',midx_filename='Halos/ds14_a_halos_1.0000.midx9',bounding_box=bbox)
    ad = ds.all_data()
    M = np.asarray(ad['m200c'])
    
    
    
    X = np.asarray(ad['x'].in_units('Mpc'))
    Y = np.asarray(ad['y'].in_units('Mpc'))
    Z = np.asarray(ad['z'].in_units('Mpc'))
    xmin, ymin, zmin = np.min(X), np.min(Y), np.min(Z)
    xmax, ymax, zmax = np.max(X), np.max(Y), np.max(Z)
    X = X-xmin
    Y = Y-ymin
    Z = Z-zmin
    
    den1 = np.sum(M)/((xmax-xmin)*(ymax-ymin)*(zmax-zmin))
    print 'Extents', xmax-xmin, ymax-ymin, zmax-zmin
    print 'Total Density', den1  
      
    for vname in xmed.keys():
        try:
            print 'Now Void', vname
            #if len(distdict[vname]) < 3.:
                #continue
            #if np.max(distdict[vname])>400.:
                #continue
            xlbs, ylbs, zlbs, xubs, yubs, zubs = [],[],[],[],[],[]
            cdcs=[]
            for i in range(len(xlistdict[vname])):
                rad = radlistdict[vname][i]
                rads.append(rad)
                pick = (np.absolute(X-xlistdict[vname][i]) < rad)*(np.absolute(Y-ylistdict[vname][i]) < rad)*(np.absolute(Z-zlistdict[vname][i]) < rad)
                Msel = M[pick]
                Xsel, Ysel, Zsel = X[pick], Y[pick], Z[pick]
                
                xmin, ymin, zmin = np.min(Xsel), np.min(Ysel), np.min(Zsel)
                xmax, ymax, zmax = np.max(Xsel), np.max(Ysel), np.max(Zsel)
                xlbs.append(xmin -2.*rad), ylbs.append(ymin -2*rad), zlbs.append(zmin - 2.*rad)
                xubs.append(xmax +2.*rad), yubs.append(ymax +2*rad), zubs.append(zmax + 2.*rad)
                print 'Extents', xmax-xmin, ymax-ymin, zmax-zmin
                den2 = np.sum(Msel)/((xmax-xmin)*(ymax-ymin)*(zmax-zmin))
                print 'Density in region, ratio', den2, den2/den1
                cdcs.append(den2/den1 -1.)
                centdencontr.append(den2/den1 - 1.)
                #np.savetxt('FetchedVoids/'+str(index)+'_'+str(vname)+'.txt', np.vstack([Xsel, Ysel, Zsel, Msel]), delimiter='|')
                #print 'FetchedVoids/'+str(index)+'_'+str(vname)+'.txt'
            XLB, YLB, ZLB = np.min(xlbs), np.min(ylbs), np.min(zlbs)
            XUB, YUB, ZUB = np.max(xubs), np.max(yubs), np.max(zubs)
            ftail='.txt'
            if ((XLB < np.min(X)) or (YLB < np.min(Y)) or (ZLB < np.min(Z)) or (XUB > np.max(X)) or (YUB > np.max(Y)) or (ZUB > np.max(Z))):
                print 'Bounds breached'
                ftail = 'BB.txt'
            fpick = (X>XLB)*(Y>YLB)*(Z>ZLB)*(X<XUB)*(Y<YUB)*(Z<ZUB)
            Mf = M[fpick]
            Xf, Yf, Zf = X[fpick], Y[fpick], Z[fpick]
            d3 = np.sum(Mf)/((np.max(Xf)-np.min(Xf))*(np.max(Yf)-np.min(Yf))*(np.max(Zf)-np.min(Zf)))
            print 'Density in all subvoid region, ratio', d3, d3/den1
            np.savetxt('FetchedVoids2/'+str(index)+'_'+str(vname)+ftail, np.vstack([Xf, Yf, Zf, Mf]), delimiter='|')
            np.savetxt('FetchedVoids2/'+str(index)+'_'+str(vname)+'VoidMembership'+ftail, np.vstack([vlistdict[vname], xlistdict[vname], ylistdict[vname], zlistdict[vname], radlistdict[vname], np.asarray(cdcs)]), delimiter='|')
        except:
            print 'Fucked Void, index', vname, index
    
    
for index in range(0,100):
    FetchRegion(index)
        
print 'DC and Radius Array lengths', len(centdencontr), len(rads)


plt.figure(0)
plt.scatter(rads, centdencontr)
plt.xlabel('Radius[Mpc]')
plt.ylabel('Central Density Contrast')
plt.savefig('CentDensContrastRad.png')
plt.show()
