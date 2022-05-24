import sys
import configparser as cp
import numpy as np
import scipy.io
import os
import time
import multiprocessing
from myGeometryFunctions import *


def main():
	
	initFile = sys.argv[1]
	simType = initFile.split('_')[0]
    
	
	config = cp.ConfigParser(interpolation=cp.ExtendedInterpolation())
	config.read(initFile)
    
	if simType == 'base':
		from driftDiffusionSimulatorBase import driftDiffusionSimulatorBase
		geometryIn = config['DEFAULT'].get('geometry filename','')
        
	elif simType == 'constriction':
		from driftDiffusionSimulatorConstriction import driftDiffusionSimulatorConstriction
		boxL = config['DEFAULT'].getfloat('bounding box', 3)
		constrWidth,constrThick = [float(s) for s in config['DEFAULT'].get(
					'Constriction Width and Thickness','0.3,.05').split(',')]
        
	elif simType == 'fullCircle':
		from driftDiffusionSimulatorFullCircle import driftDiffusionSimulatorFullCircle
		constrictionWidth = config['DEFAULT'].getfloat('Constriction Width', 0.3)
		diameter = config['DEFAULT'].getfloat('diameter', 3)
		injectH,injectW = [float(s) for s in config['DEFAULT'].get(
					'injector height,width','0.5,1').split(',')]
        
	elif simType == 'teslaValve':
		from driftDiffusionSimulatorTeslaValve import driftDiffusionSimulatorTeslaValve
		currentDirection = config['DEFAULT'].get('current direction?', 'easy')
        
	elif simType == 'periodicBC':
		from driftDiffusionSimulatorTeslaValve import driftDiffusionSimulatorTeslaValve
        
	elif simType == 'Axolotl':
		from driftDiffusionSimulatorPeriodicBC import driftDiffusionSimulatorPeriodicBC

	Temp = config['DEFAULT'].getfloat('Temperature', 70)
	sourceDrainRatio = config['DEFAULT'].getfloat('Source Drain Ratio', 1.2)
	expPartDensity = config['DEFAULT'].getfloat('Experimental Particle Density', 2e12)
	simPartDensity = config['DEFAULT'].getfloat('Simulation Particle Density', 20)
	pScatter = config['DEFAULT'].getfloat('Scattering Probability', 0)
	diffusiveEdges = config['DEFAULT'].getboolean('diffusive edges?',False)
    
	tip = config['DEFAULT'].getboolean('add probe tip?', True)
	if tip:
		probeCenterX,probeCenterY = [float(s) for s in config['DEFAULT'].get(
					'probe tip location X,Y','0,0').split(',')]
		print("Probe tip location:", probeCenterX, probeCenterY)    
    
	reinjectWithE1 = config['DEFAULT'].getboolean('reinject partilces with unity energy?', False)  
	stepSize = config['DEFAULT'].getfloat('step size', 0.01)
	collisionDist = config ['DEFAULT'].getfloat('collision overlap', 0.05)
	fieldRes = config['DEFAULT'].getfloat('histogram resolution',0.1)
    
	outPath = config['DEFAULT'].get('base output path', './')
	Nsteps = config['DEFAULT'].getint('time steps', 10000)
	dNsave = config['DEFAULT'].getint('save interval', 1000)
	Ncpu = config['DEFAULT'].getint('Number of CPUs', 1)
	initCondFile = config['DEFAULT'].get('inital conditions','')
	makeMovie = diffusiveEdges = config['DEFAULT'].getboolean('make movie?',False)
    
	generateReinjectionProbsYN = config['DEFAULT'].getboolean('generate reinjection probabilities?', False)
	setProbsMethod = config['DEFAULT'].get('set probs method','by length')
	generatedProbs = config['DEFAULT'].get('generated probs by simulation', '[]')
	generatedProbs=np.array(eval(generatedProbs))
    
	splitName = initFile.split('.')
	splitName.remove('txt')
	fnameBase = '.'.join(splitName)
	
	print(outPath+fnameBase)
    
	if not os.path.isdir(outPath):
		os.mkdir(outPath)
    
	if not os.path.isdir(outPath+fnameBase):
		os.mkdir(outPath+fnameBase)
 
	if generateReinjectionProbsYN:
		print("This simulation will be used to calculate reinjection probabilities.")
		csvFile = outPath+fnameBase+'/injectionStats.csv'

                
	dSims = []
	for iterationNum in range(Ncpu):
        
		if simType == 'base':
			dSim = driftDiffusionSimulatorBase()
			if geometryIn != '':
				mat = np.load(geometryIn)
				dSim.loadBody(mat['borderX'], mat['borderY'], mat['edgeStyle'])
				Lset = 2*np.max(mat['boxRange'])
				dSim.setBoundaryLength(Lset)
				dSim.boxRange = mat['boxRange']

		elif simType == 'constriction':
			dSim = driftDiffusionSimulatorConstriction()
			dSim.setConstrictionWidth(constrWidth)
			dSim.setWallThickness(constrThick)
			dSim.setBoundaryLength(boxL)

		elif simType == 'fullCircle':
			dSim = driftDiffusionSimulatorFullCircle()
			dSim.setConstrictionWidth(constrictionWidth)
			dSim.setInjectorShape(injectW,injectH)
			dSim.setDiameter(diameter)

		elif simType == 'teslaValve':
			dSim = driftDiffusionSimulatorTeslaValve()
			dSim.currentDirection = currentDirection
			dSim.setCurrentDirection()
            
		elif simType == 'periodicBC':
			dSim = driftDiffusionSimulatorPeriodicBC()

		elif simType == 'Axolotl':
			dSim = driftDiffusionSimulatorAxolotl()

		dSim.tip = tip
		if dSim.tip:
			dSim.probeCenterX, dSim.probeCenterY = probeCenterX, probeCenterY
		dSim.updateBody()

		dSim.generateReinjectionProbsYN = generateReinjectionProbsYN
		dSim.setProbsMethod = setProbsMethod
		dSim.setSourceDrainRatio(sourceDrainRatio)
		dSim.generatedProbs = generatedProbs
        
		if iterationNum ==0:        
			dSim.makeMovie = makeMovie
        
		dSim.updateScatterProb(pScatter)
		dSim.Temp = Temp
		dSim.expPartDensity =expPartDensity
        
        # set Emin(T) to be width of Fermi dist; set n(T,nExp,nSim) to get right lee
		#dSim.setEmin(dSim.expPartDensity, dSim.Temp)
		#dSim.setSimPartDensity(dSim.expPartDensity, dSim.Temp)
        
        # choose density and set Emin(T,nExp,nSim) to get right lee
		dSim.simPartDensity= simPartDensity
		dSim.setEmin(dSim.expPartDensity, dSim.Temp,simPartDensity=dSim.simPartDensity)
        
		dSim.calcArea()
		dSim.setNpart(dSim.simPartDensity, dSim.Area)
		dSim.updateNparticles()
        
		dSim.calc_pScatter() # if input pScatter = 0, no bulk scattering, else is t-dependent
		if pScatter ==0: 
			dSim.updateScatterProb(pScatter)
            
		if reinjectWithE1:
			dSim.consumeAndReinject = dSim.consumeAndReinject_withE1
        
		dSim.setFieldResolution(fieldRes)
		dSim.setOverlapRadius(collisionDist)
		dSim.setStepSize(stepSize)
		dSim.updateNparticles()
        
		dSim.diffusive = diffusiveEdges
		if diffusiveEdges:
			dSim.diffusiveWalls()
			print("Diffusive Edges")
		else:
			dSim.mirrorWalls()
			print("Mirror Edges")
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
				dSim.vX = mat['vY']
				dSim.overlaps = mat['overlaps']
				dSim.i_lookup = mat['i_lookup']
				dSim.j_lookup = mat['j_lookup']
				dSim.timeCount = mat['timeCount']
				dSim.frameNum = mat['frameNum']
				dSim.NcornersCrossed = mat['NcornersCrossed']
                
				dSim.Nrestarts = mat['Nrestarts']
				dSim.Nrestarts += 1
				if iterationNum ==0:
					print("Number of Restarts:", dSim.Nrestarts)
					print("timesteps:", dSim.timeCount)
				dSim.symmetrizedNinjected_NEW = mat['symmetrizedNinjected_NEW']
				dSim.symmetrizedNinjected_OLD = mat['symmetrizedNinjected_OLD']
				dSim.symmetrizedNinjected_totalOLD = mat['symmetrizedNinjected_totalOLD']
				dSim.symmetrizedNinjected_DIF = mat['symmetrizedNinjected_DIF']
            
		dSims.append(dSim)

	print("Temp SourceDrain Emin expPartDensity simPartDensity pScatter")
	print(dSim.Temp, dSim.sourceDrainRatio, dSim.Emin, dSim.expPartDensity,dSim.simPartDensity, dSim.p_scatter)

	jobs=[]
	for iterationNum,dSim in enumerate(dSims):
        
		fnameOut = outPath+fnameBase+"/"+fnameBase+("_%03d"%iterationNum)+".npz"	
		dSim.saveLoc = outPath+fnameBase
	
		p = multiprocessing.Process(target=dSim.runAndSave, args=(Nsteps-dSim.timeCount,dNsave,fnameOut))
		jobs.append(p)
		p.start()
	p.join()

		


if __name__== "__main__":
  main()
