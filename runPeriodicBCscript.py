import sys
import configparser as cp
import numpy as np
import scipy.io
import os
import time
import multiprocessing
from driftDiffusionSimulatorPeriodicBC import driftDiffusionSimulatorPeriodicBC
from myGeometryFunctions import *

'''
runPeriodicBCscript.py reads in a script text file and runs a number of simulations of
a "driftDiffusionSimulatorBase" simulation object. It enables starting from a 
fresh simulation, or taking the results of a previous simulation as initial cond-
itions. 
'''

_fieldResolution = 0.1
def main():
	
	initFile = sys.argv[1]
	
	config = cp.ConfigParser(interpolation=cp.ExtendedInterpolation())
	config.read(initFile)

	pScatter = config['DEFAULT'].getfloat('Scattering Probability', 0)
	Emin = config['DEFAULT'].getfloat('Minimum Scattering Energy', 0.9)
	Temp = config['DEFAULT'].getfloat('Temperature', 100)
	expPartDensity = config['DEFAULT'].getfloat('Experimental Particle Density', 1e12)
	simPartDensity = config['DEFAULT'].getint('Simulation Particle Density', 100)
	Nsteps = config['DEFAULT'].getint('time steps', 10000)
	dNsave = config['DEFAULT'].getint('save interval', 1000)
	collisionDist = config ['DEFAULT'].getfloat('collision overlap', 0.05)
	stepSize = config['DEFAULT'].getfloat('step size', 0.01)
	boxL = config['DEFAULT'].getfloat('Box Length', 3.0)
	
	geometryIn = config['DEFAULT'].get('geometry filename','')

	Ncpu = config['DEFAULT'].getint('Number of CPUs', 1)
	outPath = config['DEFAULT'].get('base output path', './SIM_data/')
	
	initCondFile = config['DEFAULT'].get('inital conditions','')
    
	
	splitName = initFile.split('.')
	splitName.remove('txt')
	fnameBase = '.'.join(splitName)
    
	print(outPath+fnameBase)
	if not os.path.isdir(outPath+fnameBase):
		os.mkdir(outPath+fnameBase)
            

	dSims = []
	for iterationNum in range(Ncpu):
		
		#set up the simulation from script file
		dSim = driftDiffusionSimulatorPeriodicBC()
		dSim.updateScatterProb(pScatter)
		dSim.Temp = Temp
		dSim.Emin = Emin
		dSim.expPartDensity = expPartDensity
		dSim.simPartDensity = simPartDensity
		dSim.setOverlapRadius(collisionDist)
		dSim.setStepSize(stepSize)
		dSim.setBoxSize(boxL)
        
		if geometryIn != '':
			mat = np.load(geometryIn)
			dSim.loadBody(mat['borderX'], mat['borderY'], mat['edgeStyle'])
			Lset = 2*np.max(mat['boxRange'])
			dSim.setBoundaryLength(Lset)
			dSim.boxRange = mat['boxRange']
			dSim.setFieldResolution(_fieldResolution)
			if dSim.tip:
			    dSim.addTip()
			
		else:
			dSim.updateBody()
   
		dSim.calcArea()
		dSim.setNpart(dSim.simPartDensity, dSim.Area)
		dSim.updateNparticles()
        
	
		for i in range(int(dSim.Npart/2)):
			thetas = np.random.rand()*2.*np.pi
			dSim.vX[i] = np.cos(thetas)
			dSim.vY[i] = np.sin(thetas)
			dSim.vX[i+int(dSim.Npart/2)] = -np.cos(thetas)
			dSim.vY[i+int(dSim.Npart/2)] = -np.sin(thetas)
		
		if initCondFile:
			if os.path.isfile(outPath+initCondFile+'/'+initCondFile+("_%03d"%iterationNum)+".npz"):
				mat = np.load(outPath+initCondFile+'/'+initCondFile+("_%03d"%iterationNum)+".npz")
				dSim.Nabsorbed = mat['Nabsorbed']
				dSim.Ninjected = mat['Ninjected']
				dSim.Px = mat['Px']
				dSim.Py = mat['Py']
				dSim.pR = mat['pR']
				dSim.Erho = mat['Erho']
				dSim.rho = mat['rho']
				dSim.Xpos = mat['Xpos']
				dSim.Ypos = mat['Ypos']
				dSim.vX = mat['vX']
				dSim.vX = mat['vX']
				dSim.overlaps = mat['overlaps']
				dSim.i_lookup = mat['i_lookup']
				dSim.j_lookup = mat['j_lookup']
				dSim.timeCount = mat['timeCount']
				dSim.NcornersCrossed = mat['NcornersCrossed']
				dSim.traces = mat['traces']
				dSim.L_NoScatter = mat['L_NoScatter']
                
				dSim.Nrestarts = mat['Nrestarts']
				dSim.Nrestarts += 1
				if iterationNum ==0:
					print("Number of Restarts:", dSim.Nrestarts)
					print("timesteps:", dSim.timeCount)

			
		dSims.append(dSim)        

	jobs=[]
	for iterationNum,dSim in enumerate(dSims):
		fnameOut = outPath+fnameBase+"/"+fnameBase+("_%03d"%iterationNum)+".npz"
	
		p = multiprocessing.Process(target=dSim.runAndSave, args=(Nsteps-dSim.timeCount,dNsave,fnameOut))
		jobs.append(p)
		p.start()
	p.join()


if __name__== "__main__":
  main()
