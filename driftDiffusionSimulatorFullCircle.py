import numpy as np
import relativisticTransforms as rT
import matplotlib.path as mpltPath
from myGeometryFunctions import *	
from driftDiffusionSimulatorBase import driftDiffusionSimulatorBase

class driftDiffusionSimulatorFullCircle(driftDiffusionSimulatorBase):
	'''
	driftDiffusionSimulatorFullCircle simulates hydrodynamic flow of dirac electrons
	through a constriction and into a circular device. The enitre arc of the circle
	is the drain and a source injects through a constriction at the arc's radial center
	'''
	def __init__(self):
	
		self.constrictionWidth = 1.
		self.sourceDrainRatio = 1.
		self.injectorWidth = .6
		self.injectorHeight = 1.5
		self.diameter = 3.
		self.diffusive = False
	 
		super().__init__()
        
		self.rollIndex = 17
        
		self.diameter=self.boxL
		self.setDiameter(self.boxL)
		self.setConstrictionWidth(self.constrictionWidth)
		self.setInjectorShape(self.injectorWidth, self.injectorHeight)
        
		self.setBoundaryLength(self.diameter*1.05)
        
		self.countsPerSnapshot = 500
        
	def updateBody(self):
		R = self.diameter/2.
		xc = self.constrictionWidth/2.
		xi = self.injectorWidth/2.
		yi = self.injectorHeight
		yw = R
		
		thetas = np.linspace(1.4*np.pi,-.4*np.pi,25)
		
		self.borderX = np.cos(thetas)*R
		self.borderY = np.sin(thetas)*R
		self.borderX = np.append(self.borderX,
									np.array([xc,xi,-xi,-xc,self.borderX[0]]))
		self.borderY = np.append(self.borderY,
									np.array([0,-yi,-yi,0,self.borderY[0]]))
			
		s,d,m,r,f = 2,1,0,-1,-2 #'source','drain','mirror','rough'
		
		self.edgeStyle = [d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,m,m,s,m,m]
		self.nDeviceEdges=len(self.edgeStyle)
		self.boxRange = [[-R, R], [-R, R]]
    
		self.setUpEdges()
        
		if self.tip:
			self.addTip()

	def setInjectorShape(self,width_in, height_in):
		self.injectorWidth = width_in
		self.injectorHeight = height_in
		self.updateBody()
	
	def setConstrictionWidth(self, width_in):
		self.constrictionWidth = width_in
		self.updateBody()
		
		
	def setDiameter(self, D):
		self.boxL = D
		self.diameter = D
		self.updateBody()

