import numpy as np
import relativisticTransforms as rT
import matplotlib.path as mpltPath
from myGeometryFunctions import *	
from driftDiffusionSimulatorBase import driftDiffusionSimulatorBase

class driftDiffusionSimulatorAxolotl(driftDiffusionSimulatorBase):
	'''
	driftDiffusionSimulatorAxolotl simulates hydrodynamic flow of dirac electrons
    through a uniquely shaped device in which we have observed a number of affects
    attributable to hydrodynamics. We simulate with a movable tip to mimic scanning gate
    microscopy.
	'''
	def __init__(self):
	
		self.sourceDrainRatio = 1.2
	 
		super().__init__()
        
		self.countsPerSnapshot = 1000
		self.setBoundaryLength(7)
		self.setFieldResolution(.05)
        
	def updateBody(self):

		self.borderX=np.array([0.4525666316681676, 0.46998095899451153, 0.04999999999999982, 0.04999999999999982,
                               -0.050000000000000044, -0.050000000000000044, -0.625, -0.625, -1.1612311101832846,
                               -1.3219280126049195, -0.875, -0.875, -1.8, -1.8, -0.875, -0.875, -1.3219280126049195, 
                               -1.1612311101832846, -0.625, -0.625, -0.050000000000000044, -0.050000000000000044, 
                               0.04999999999999982, 0.04999999999999982, 0.46998095899451153, 0.4525666316681676,
                               0.4351142252308841, 0.6850761490199819, 0.7025285554572654, 0.7200190410054884, 0.78, 
                               0.78, 1.0388190451025208, 1.280300501674788, 1.03, 1.03, 3.939185328635486, 
                               5.1181970119322635, 5.161270655364881, 3.95, 3.95, 5.141660666982949,
                               5.1560550648135095, 3.95, 3.95, 5.1560550648135095, 5.141660666982949, 3.95, 3.95, 
                               5.161270655364881, 5.1181970119322635, 3.939185328635486, 1.03, 1.03, 
                               1.280300501674788, 1.0388190451025208, 0.78, 0.78, 0.7200190410054884, 
                               0.7025285554572654, 0.6850761490199819, 0.4351142252308841, 0.4525666316681676])
        
		self.borderY = np.array([-1.7476661443517307, -0.75, -0.75, -0.2, -0.2, -0.75, -0.75, -1.45, 
                                 -1.8999513267805774, -1.708440216000833, -1.3334230854612503, -0.75, -0.75, 0.75,
                                 0.75, 1.3334230854612503, 1.708440216000833, 1.8999513267805774, 1.45, 0.75, 0.75,
                                 0.2, 0.2, 0.75, 0.75, 1.7476661443517307, 2.7475138395081222, 2.751876941117443,
                                 1.7520292459610518, 0.75, 0.75, 1.75, 2.715925826289068, 2.6512210650134382, 
                                 1.717086875603151, 0.75, 0.75, 0.9052198327907406, 0.5780430285373831, 
                                 0.41857622150427476, 0.3601571930519477, 0.4121862213994097, 0.08250030827739663,
                                 0.029842806948052168, -0.029842806948052282, -0.08250030827739674, 
                                 -0.4121862213994098, -0.3601571930519478, -0.41857622150427476, 
                                 -0.5780430285373831, -0.9052198327907406, -0.75, -0.75, -1.717086875603151,
                                 -2.6512210650134382, -2.715925826289068, -1.75, -0.75, -0.75, -1.7520292459610518,
                                 -2.751876941117443, -2.7475138395081222, -1.7476661443517307])
        
		self.borderX = self.borderX-np.average(self.borderX)
		self.borderY = self.borderY-np.average(self.borderY)
			
		s,d,m,r,f = 2,1,0,-1,-2 #'source','drain','mirror','rough'
		
		self.edgeStyle = [m,m,m,m,m,m,m,m,m,m,
                m,m,s,m,m,m,m,m,m,m,
                m,m,m,m,m,m,m,m,m,m,
                m,m,m,m,m,m,m,d,m,m,
                m,d,m,m,m,d,m,m,m,d,
                m,m,m,m,m,m,m,m,m,m,
                m,m]
		self.nDeviceEdges=len(self.edgeStyle)
        
		self.setUpEdges()
        
		if self.tip:
			self.addTip()
        
            
            
		
