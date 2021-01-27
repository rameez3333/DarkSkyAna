from glob import glob
import os, sys
import numpy as np
import time
nside = str(32)


flines = open('submittemplatearbmem.sh').readlines()


#speclist = ['05', '06', '13', '25', '39', '43', '46' ,'48']
#flist = sorted(glob('WISE/wise-allsky-cat-part??'))
#flist = sorted(glob('AllWISE/wise-allwise-cat-part??'))
#flist = ['AllWISE/wise-allwise-cat-part46', 'AllWISE/wise-allwise-cat-part47', 'AllWISE/wise-allwise-cat-part48']
jname = "DarkSkyQueryVC"

tnjobs = 250

h = 0.688062

edgeclearance = 2200.

totallength = 8000./h

#memreq = '40G'

maxjobs = 50

memreqd={0:'20G', 1:'20G', 2:'20G', 3:'20G', 4:'20G', 5:'20G', 6:'20G', 7:'20G', 8:'20G', 9:'20G', 10:'40G'}
#memreqd[0] = '20G's

def countqueue():
    time.sleep(10)
    os.system('squeue -p icecube_guest > countqueue')
    njr = len(open('countqueue').readlines()) -1
    return njr -1



for s in [10, 0]:
    flist = glob('VelCovResultsBF/DarkSkyVelCovPanthQuery_CenterX_*_Y_*_Z_*_W_*_OBS_'+str(s)+'_2MppFlow*_Out.txt')
    print 'S:', s, len(flist)
    njobs = tnjobs
    #if njobs<0:
        #njobs=0
    i = 0
    while i < njobs:
        if countqueue() < maxjobs:
            CX = np.random.uniform(edgeclearance, totallength-edgeclearance)
            CY = np.random.uniform(edgeclearance, totallength-edgeclearance)
            CZ = np.random.uniform(edgeclearance, totallength-edgeclearance)
            jobname = jname+str(i)+str(CX)+str(CY)+str(CZ)+str(s)
            fout = open(jobname+'.slurm', "w")
            jobline = 'python AnalyzeBox_VelCovPanth_extVP_TO.py -d '+' -x '+str(CX)+' -y '+str(CY)+' -z '+str(CZ) +' -s '+str(s)
            #os.system(jobline)
            for line in flines:
                fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline).replace('__MEMR__', memreqd[s])  )
            fout.close()
            os.system('chmod +x ' + jobname+'.slurm')
            os.system('sbatch -p icecube_guest '+ jobname+'.slurm')
            i = i+1
        else:
            time.sleep(100)
        #time.sleep(10)

        #raw_input("Press Enter to Continue")
