import numpy as np
import relativisticTransforms as rT
import matplotlib.path as mpltPath
from myGeometryFunctions import *	
from driftDiffusionSimulatorBase import driftDiffusionSimulatorBase
from scipy.spatial.distance import pdist,squareform

class driftDiffusionSimulatorPeriodicBC(driftDiffusionSimulatorBase):
	'''
	driftDiffusionSimulatorPeriodicBC simulates electron collisons in a box with periodic boundary conditions. We use the reuslts of this simulation to map our simulaiton parameters to temperature.
	'''
	def __init__(self):
	
		super().__init__()
        
        
	def updateBody(self):
		xw = self.boxL/2.
		yw = self.boxL/2.
		
        # left border
		self.borderX = np.array([-xw,-xw,xw,xw,-xw])
		self.borderY = np.array([-yw,yw,yw,-yw,-yw])
		
		self.edgeStyle = [1,1,1,1]

		self.nDeviceEdges=len(self.edgeStyle)
        
		self.boxRange = [[-self.boxL, self.boxL], [-self.boxL, self.boxL]]
        
		self.setUpEdges()
        
		self.setFieldResolution(self.fieldResolution)
        
		if self.tip:
			self.addTip()	
            
		self.setUpEdges()
    
	def setBoxSize(self, L):
		self.boxL = L
		self.updateBody()
	
		
	def setFieldResolution(self,dX):
		self.fieldResolution = dX
		Nr = int(np.round(2*self.boxL/dX))
		
		self.Nx = Nr
		self.Ny = Nr
		
		self.rho = np.zeros((self.Nx,self.Ny))
		self.Px = np.zeros((self.Nx,self.Ny))
		self.Py = np.zeros((self.Nx,self.Ny))
		self.Erho = np.zeros((self.Nx,self.Ny))
		
		_,self.histX,self.histY = np.histogram2d(np.zeros(10),np.zeros(10),
				bins = [self.Nx,self.Ny],range = self.boxRange)
        
	def timestep(self):
		# executes a single timestep of the simulation
		
		#propagate all particles in the direction of their velocity vector
		self.Xpos += self.vX*self.DeltaX
		self.Ypos += self.vY*self.DeltaX
        
		#update mean free path statistics
		self.L_NoScatter+=self.DeltaX

		injectedPtcs = []
        
		for i in range(self.Npart):
			
			if self.Xpos[i]>self.boxL/2.:
				self.Xpos[i] -= self.boxL
			if self.Xpos[i]<-self.boxL/2.:
				self.Xpos[i] += self.boxL
			if self.Ypos[i]>self.boxL/2.:
				self.Ypos[i] -= self.boxL
			if self.Ypos[i]<-self.boxL/2.:
				self.Ypos[i] += self.boxL
						
			elif self.p_scatter > np.random.rand(1):
				#randomize direction in bulk
				theta = np.random.rand(1)*2.*np.pi
                
				self.vX[i] = np.cos(theta)
				self.vY[i] = np.sin(theta)
		
		
		_overlaps = self.checkOverlaps()
		for i in injectedPtcs:
			self.overlaps[:,i] = True
			self.overlaps[i,:] = True
		self.overlaps[np.array(squareform(np.invert(_overlaps)), dtype=bool)] = False
		
		_idx_i = self.i_lookup[_overlaps]
		_idx_j = self.j_lookup[_overlaps]
		
		for i,j in zip(_idx_i,_idx_j):
			if not self.overlaps[i,j]:
				#perform relativistic scattering
				p1 = np.array([self.vX[i]*self.pR[i],
										self.vY[i]*self.pR[i]])
				p2 = np.array([self.vX[j]*self.pR[j],
										self.vY[j]*self.pR[j]])
				#calculate outgoing momenta
				p3,p4 = rT.scatterMasslessParticles(
												p1,p2,self.Emin)
				if not (np.array_equal(p1,p3) and 
										np.array_equal(p2,p4)):
                    
					self.traces[i].append([p1[0],p1[1],self.L_NoScatter[i]])
					self.traces[j].append([p2[0],p2[1],self.L_NoScatter[j]])
				
					#measure momentum amplitude
					self.pR[i] = np.sqrt(np.sum(p3**2))
					self.pR[j] = np.sqrt(np.sum(p4**2))
					#set velocities
					self.vX[i],self.vY[i] = p3/self.pR[i]
					self.vX[j],self.vY[j] = p4/self.pR[j]
					
					self.L_NoScatter[i]=0
					self.L_NoScatter[j]=0
				
			self.overlaps[i,j] = True
			
				
		#update global time index
		self.timeCount+=1
    
		
