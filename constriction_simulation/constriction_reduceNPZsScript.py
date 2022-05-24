import numpy as np
import os
import sys
from zipfile import BadZipFile

eHydroDirectory = '/'.join(np.delete(sys.path[0].split('/'), -1))
sys.path.append(eHydroDirectory)
from driftDiffusionSimulatorConstriction import driftDiffusionSimulatorConstriction

def main(arguments):

    oldDirectory = arguments[1]+'/'
    
    # reduces all files in oldDirectory
    
    for fnameBase in os.listdir(oldDirectory): # finds directories
        if os.path.isdir(oldDirectory+fnameBase):
            print(fnameBase)

            Nsims = 0 # counts Nsims
            pathOld = oldDirectory+fnameBase+'/' # path to NPZs
            for file in os.listdir(pathOld):
                if fnameBase in file: # only looks at relevant NPZs
                    if os.path.isfile(os.path.join(pathOld, file)):
                        Nsims +=1

            pathNew = arguments[1]+'_converted/'+fnameBase+'/'
            os.makedirs(pathNew, exist_ok=True)
            
            flag = 1
            for i in range(Nsims):
                try:
                    old_fname = pathOld+fnameBase+'_%03d.npz'%i
                    with np.load(old_fname,allow_pickle=True) as mat: # loads in large NPZ
                        # and generates smaller (excludes 4 variables)
                        new_fname = pathNew+fnameBase+'_%03d.npz'%i
                        saveFunction = "np.savez(new_fname, "
                        for variable in mat:
                                if variable != 'borderPath' and variable != 'i_lookup' and variable != 'j_lookup' and variable != 'overlaps':
                                    array = 'mat['+'"'+variable+'"'+']'
                                    saveFunction = saveFunction+variable+ '='+array+','
                        saveFunction = saveFunction[:-2]+'])'
                        exec(saveFunction)
                        flag=0
                
                except BadZipFile:
                    print("Zip file error handled")
                    os.system('rm '+old_fname)
                    
            # if all simulations return errors, remove their old and new directories
            if flag==1:
                os.system('rm -r '+oldDirectory+fnameBase)
                os.system('rm -r '+pathNew)

if __name__== "__main__":
    main(sys.argv)

