import numpy as np
import relativisticTransforms as rT
import matplotlib.path as mpltPath
from myGeometryFunctions import *	
from driftDiffusionSimulatorBase import driftDiffusionSimulatorBase
import os
import matplotlib.pyplot as plt

class driftDiffusionSimulatorTeslaValve(driftDiffusionSimulatorBase):
	'''
	driftDiffusionSimulatorTeslaValve simulates hydrodynamic flow of dirac electrons
    through a valve-like device where viscous flow gives rise to a path-dependent resistance
	'''
	def __init__(self):
	
		self.sourceDrainRatio = 1.
		self.currentDirection = 'easy'
        
		self.initializationCount = 2000000
	 
		super().__init__()
        
		self.setBoundaryLength(9)
		self.countsPerSnapshot = 1000
        
	def updateBody(self):
		border=np.array([[3.35300,-9.35600],
                        [-4.37600,-7.22400],
                         
                        [-7.76200,-6.25200], #bottom left voltage probe
                        [-7.30500,-3.62600],
                        [-5.66300,-5.51000],
                        [-4.61500,-6.00600],
                        [-4.04100,-5.97400],
                        [-3.62500,-5.89400],
                         
                        [-3.37300,-5.82800],
                        [-2.79300,-5.72100],
                        [-2.23100,-5.70500],
                        [-1.89600,-5.62100],
                        [-1.51900,-5.32800],
                        [-1.36900,-4.92200],
                        [-1.16000,-4.06000],
                        [-1.19600,-3.79700],
                        [-1.44700,-3.79700],
                        [-1.57900,-3.89300],
                        [-1.92600,-4.25800],
                        [-2.41600,-4.58700],
                        [-3.05000,-4.78400],
                        [-3.64200,-4.87400],
                        [-4.22800,-4.91000],
                        [-4.52100,-4.85000],
                        [-4.95700,-4.59300],
                        [-5.25600,-4.34200],
                        [-5.53200,-3.95300],
                        [-5.83100,-3.53400],
                        [-5.96200,-3.08000],
                        [-5.98600,-2.66700],
                        [-5.98000,-2.06300],
                        [-5.45600,-0.31600],
                        [-3.25200,1.19200],
                        [-2.03500,1.36700],
                        [-1.01900,1.47900],
                        [-0.12900,1.59300],
                        [0.10800,1.70500],
                        [0.25100,1.84100],
                        [0.28400,2.02400],
                        [0.29100,2.25000],
                        [0.17900,2.39700],
                        [-0.20100,2.49400],
                        [-0.81800,2.54500],
                        [-2.15000,2.64900],
                         
                        [-4.77200,2.82700], #top left voltage probe
                        [-5.37400,3.10800],
                        [-5.44900,3.67900],
                        [-7.05900,3.68800],
                        [-7.05900,4.85900],
                        [-4.98600,5.08500],
                        [-2.82800,4.71800],
                         
                        [-1.73700,4.16200],
                         
                        [-0.67400,3.92200], #upper middle contact
                        [0.20000,5.16600],
                        [0.69200,5.97100],
                        [1.18600,6.65700],
                        [2.59800,5.94400],
                        [1.89800,4.67300],
                        [1.29800,3.61800],
                        [1.47200,3.56300],

                        [4.89500,2.42500],
                        [6.02200,2.07400],
                        [6.77300,2.09500],
                        [7.09600,2.19300], #upper omhic contact
                        [7.00600,1.50200], #upper omhic contact
                        [6.22600,1.51200],
                        [5.70200,1.49500],
                        [5.37200,1.03600],
                        [5.24400,0.42000],
                        [4.57200,0.11200],
                        [4.15600,0.19200],
                         
                        [3.51800,0.31400],
                        [3.60800,0.73400],
                        [3.60800,1.28700],
                        [3.42800,1.46700],
                        [3.24900,1.54200],
                        [2.91900,1.51200],
                        [2.69500,1.42200],
                        [2.62000,1.21300],
                        [0.89800,-6.09300],
                        [0.89800,-6.40700],
                        [1.04800,-6.66200],
                        [1.18300,-6.67700],
                        [1.45200,-6.57200],
                        [1.54200,-6.33200],
                        [1.78100,-5.92800], 
                         
                        [2.23800,-4.28200], # right middle contact
                        [3.96500,-4.80500],
                        [3.20400,-6.36200],
                        [3.38300,-7.38200], 
                                    
                        [3.86200,-7.53000],
                        [4.95400,-7.90300], #lower ohmic contact
                        [4.53400,-9.51600], #lower omhic contact
                        [3.35300,-9.35600],
 
                        [-3.44000,-2.54200], #center cutout
                        [-3.24700,-2.51700],
                        [-2.42800,-2.15400],
                        [-1.47300,-1.54600],
                        [-1.09500,-1.24700],
                        [-0.70700,-0.94800],
                        [-0.56300,-0.82600],
                        [-0.51500,-0.63400],
                        [-0.51200,-0.44100],
                        [-0.58000,-0.33200],
                        [-0.72700,-0.30900],
                        [-1.31800,-0.52700],
                        [-2.17700,-0.92500],
                        [-2.91700,-1.34300],
                        [-3.38600,-1.64800],
                        [-3.61900,-1.86300],
                        [-3.65200,-2.12700], 
                        [-3.63000,-2.39500],[-3.44000,-2.54200]])
		
		self.borderX=border[:,0]#+np.ones(len(border))*.5
		self.borderY=border[:,1]+np.ones(len(border))*1.2
			
		s,d,m,r,f = 2,1,0,-1,-2 #'source','drain','mirror','rough'
		
		self.edgeStyle = [m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,d,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,s,m,
                          f,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m]
		self.nDeviceEdges=len(self.edgeStyle)
        
		self.setCurrentDirection()
        
		self.boxRange = [[-self.boxL, self.boxL], [-self.boxL, self.boxL]]
        
		self.setUpEdges()

		if self.tip:
			self.addTip()
            

	def setCurrentDirection(self):
		s,d,m,r,f = 2,1,0,-1,-2 #'source','drain','mirror','rough'
        
		if self.currentDirection == 'easy':
			self.edgeStyle = [m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,d,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,s,m,
                          f,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m]        

		elif self.currentDirection == 'hard':
			self.edgeStyle = [m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,s,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m,m,
                          m,d,m,
                          f,m,m,m,m,m,m,m,m,m,
                          m,m,m,m,m,m,m,m,m]

		self.setUpEdges()
		
	
        