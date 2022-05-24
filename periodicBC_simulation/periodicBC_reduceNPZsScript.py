import numpy as np
import os
import sys
from zipfile import BadZipFile

eHydroDirectory = '/'.join(np.delete(sys.path[0].split('/'), -1))
sys.path.append(eHydroDirectory)
from driftDiffusionSimulatorPeriodicBC import driftDiffusionSimulatorPeriodicBC

def main(arguments):

    oldDirectory = arguments[1]+'/'
    
    # reduces all files in oldDirectory
    
    for SIMname in os.listdir(oldDirectory): # finds directories
        try:
            if os.path.isdir(oldDirectory+SIMname):
                fnameBase= SIMname # directory name in SIM_data
                print(fnameBase)

                Nsims = 0 # counts Nsims
                pathOld = oldDirectory+fnameBase+'/' # path to NPZs
                for file in os.listdir(pathOld):
                    if SIMname in file: # only looks at relevant NPZs
                        if os.path.isfile(os.path.join(pathOld, file)):
                            Nsims +=1

                pathNew = arguments[1]+'Converted/'+fnameBase+'/'
                os.makedirs(pathNew, exist_ok=True)
                for i in range(Nsims):
                    with np.load(pathOld+fnameBase+'_%03d.npz'%i,allow_pickle=True) as mat: # loads in large NPZ
                        # and generates smaller (excludes 4 variables)
                        fname = pathNew+fnameBase+'_%03d.npz'%i
                        saveFunction = "np.savez(fname, "
                        for variable in mat:
                                if variable != 'borderPath' and variable != 'i_lookup' and variable != 'j_lookup' and variable != 'overlaps':
                                    array = 'mat['+'"'+variable+'"'+']'
                                    saveFunction = saveFunction+variable+ '='+array+','
                        saveFunction = saveFunction[:-2]+'])'
                        exec(saveFunction)
                    
        except BadZipFile:
            print("Zip file error handled")
            os.system('rm '+oldDirectory+SIMname+'/*')

if __name__== "__main__":
    main(sys.argv)
