from glob import glob
import os, sys
import numpy as np
import time
nside = str(32)


#flines = open('submittemplateCijInteg.sh').readlines()
flines = open('submittemplate.sh').readlines()

#speclist = ['05', '06', '13', '25', '39', '43', '46' ,'48']
#flist = sorted(glob('WISE/wise-allsky-cat-part??'))
#flist = sorted(glob('AllWISE/wise-allwise-cat-part??'))
#flist = ['AllWISE/wise-allwise-cat-part46', 'AllWISE/wise-allwise-cat-part47', 'AllWISE/wise-allwise-cat-part48']
#jname = "CijIntegrator"

jname = "CijCompiler"

njobs = int(sys.argv[1])

#tot = 11476 #JLA

tot = 345696 #Pantheon


nperjob = tot/njobs

rem = tot%njobs

for i in range(njobs):
        jobname = jname+str(i)
        fout = open(jobname+'.slurm', "w")
        #jobline = 'python CijIntegratorForCompiledVelcovPlotterPanth.py '+ str(int(nperjob*i)) + ' ' + str(int(nperjob))
        jobline = 'python VelCovCompilerPanthbatch.py '+ str(int(nperjob*i)) + ' ' + str(int(nperjob))
        #os.system(jobline)
        for line in flines:
            fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
        fout.close()
        os.system('chmod +x ' + jobname+'.slurm')
        os.system('sbatch -p icecube_guest '+ jobname+'.slurm')
        time.sleep(0.01)
        #raw_input("Press Enter to Continue")

if rem:
    i=i+1
    jobname = jname+str(i)
    fout = open(jobname+'.slurm', "w")
    #jobline = 'python CijIntegratorForCompiledVelcovPlotterPanth.py '+ str(int(nperjob*i)) + ' ' + str(int(rem))
    jobline = 'python VelCovCompilerPanthbatch.py '+ str(int(nperjob*i)) + ' ' + str(int(rem))
    #os.system(jobline)
    for line in flines:
        fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
    fout.close()
    os.system('chmod +x ' + jobname+'.slurm')
    os.system('sbatch -p icecube_guest '+ jobname+'.slurm')
    time.sleep(0.01)
