import numpy as np
import relativisticTransforms as rT
import matplotlib.path as mpltPath
from myGeometryFunctions import *	
from driftDiffusionSimulatorBase import driftDiffusionSimulatorBase

class driftDiffusionSimulatorAxolotl(driftDiffusionSimulatorBase):
	'''
	driftDiffusionSimulatorTeslaValve simulates hydrodynamic flow of dirac electrons
    through a valve-like device where viscous flow gives rise to a path-dependent resistance
	'''
	def __init__(self):
	
		self.sourceDrainRatio = 1.2
	 
		super().__init__()
        
		self.countsPerSnapshot = 1000
		self.setBoundaryLength(6)
        
	def updateBody(self):

		self.borderX=np.array([-0.667, -0.91 , -0.891, -1.467, -1.458, -1.559, -1.568, -2.142,
       -2.155, -2.869, -3.025, -2.402, -2.392, -3.317, -3.29 , -2.365,
       -2.354, -2.974, -2.81 , -2.103, -2.116, -1.541, -1.549, -1.448,
       -1.443, -0.866, -0.848, -0.573, -0.333, -0.598, -0.617,  2.197,
        3.379,  3.416,  2.202,  2.202,  3.427,  3.435,  2.197,  2.195,
        3.428,  3.408,  2.19 ,  2.188,  3.396,  3.348,  2.171, -0.641,
       -0.659, -0.425, -0.667])
		self.borderY = np.array([-2.623, -1.653, -0.654, -0.642, -0.141, -0.141, -0.642, -0.631,
       -1.329, -1.912, -1.715, -1.208, -0.626, -0.61 ,  0.889,  0.872,
        1.456,  1.994,  2.183,  1.568,  0.869,  0.859,  0.359,  0.357,
        0.857,  0.846,  1.845,  2.808,  2.74 ,  1.81 ,  0.842,  0.794,
        0.928,  0.599,  0.463,  0.404,  0.434,  0.106,  0.073,  0.013,
       -0.063, -0.391, -0.316, -0.376, -0.555, -0.881, -0.707, -0.657,
       -1.623, -2.562, -2.623])
			
		s,d,m,r,f = 2,1,0,-1,-2 #'source','drain','mirror','rough'
		
		self.edgeStyle = [m,m,m,m,m,m,m,m,m,m,
                m,m,m,s,m,m,m,m,m,m,
                m,m,m,m,m,m,m,m,m,m,
                m,m,d,m,m,m,d,m,m,m,
                d,m,m,m,d,m,m,m,m,m]
		self.nDeviceEdges=len(self.edgeStyle)
        
		self.boxRange = [[-self.boxL, self.boxL], [-self.boxL, self.boxL]]
        
		self.setUpEdges()
        
		if self.tip:
			self.addTip()
            
            
		
