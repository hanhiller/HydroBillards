import numpy as np
import relativisticTransforms as rT
import matplotlib.path as mpltPath
from myGeometryFunctions import *	
from driftDiffusionSimulatorBase import driftDiffusionSimulatorBase

class driftDiffusionSimulatorConstriction(driftDiffusionSimulatorBase):
	'''
	driftDiffusionSimulatorConstriction simulates hydrodynamic flow of dirac 
	electrons through a constriction
	'''
	
	def __init__(self):
		
		self.expPartDensity = 5e11
		self.constrictionWidth = 1.
		self.wallThickness = 0.05
        
		super().__init__()
        
		self.countsPerSnapshot = 1000
        

	def updateBody(self):
		#builds a constriction with source on the bottom and drain on the top
	
		xw = self.boxL/2.
		xc = self.constrictionWidth/2.
		yc = self.wallThickness/2.
		yw = self.boxL/2.
		
		self.borderX = np.array([-xw,-xw,-xc,-xc,-xw,-xw,xw,xw,xc,xc,xw,xw,-xw])
		self.borderY = np.array([-yw,-yc,-yc,yc,yc,yw,yw,yc,yc,-yc,-yc,-yw,-yw])
		
		s = 2#'source'
		d = 1#'drain'
		m = 0#'mirror'
		r =-1#'rough'
		
		self.edgeStyle = [s,m,m,m,d,d,d,m,m,m,s,s]
		self.nDeviceEdges=len(self.edgeStyle)
		
		self.setUpEdges()
        
		if self.tip:
		    self.addTip()

	def setConstrictionWidth(self, width_in):
		self.constrictionWidth = width_in
		self.updateBody()
		
		
	def setWallThickness(self, thickness_in):
		self.wallThickness = thickness_in
		self.updateBody()
        
	
