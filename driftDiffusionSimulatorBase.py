import numpy as np
import relativisticTransforms as rT
from scipy.spatial.distance import pdist,squareform
import matplotlib.path as mpltPath
from myGeometryFunctions import *
import os
import csv
import matplotlib.pyplot as plt

import time
	
class driftDiffusionSimulatorBase:
	'''
	driftDiffusionSimulatorBase is a simulation object that serves to simulate
	electron flow in 2D. Specifically, it is designed to simulate dirac-like 
	quasiparticles (conical dispersion) and allows them to scatter off of each-
	other while conserving momentum and energy. An effective temperature is 
	included by setting a minimum energy for available scattering that is small-
	er (slightly) than the average energy of particles (roughly the Fermi ener-
	gy). Region boundaries are defined with edge types either being a "source",
	"drain","mirror", or "rough surface". Net currents are sourced by preferen-
	tially sourcing particles from drain compared with source electrode. 
	'''

	def __init__(self):
        
		self.Temp= 4
		self.expPartDensity = 1e12
		self.setEmin(self.expPartDensity, self.Temp)
		self.simPartDensity = 35
		#self.setSimPartDensity(self.expPartDensity, self.Temp)
        
		self.p_scatter = 0.01
		self.p_transmission = 1.
		self.sourceDrainRatio = 1.
		self.rOverlap = 0.05
		self.setOverlapRadius(self.rOverlap)
		self.initializationCount = 1000000
		self.setStepSize(.01)
        
		self.NumberDrainBins = 10
		self.NcornersCrossed = 0
        
		self.tip = False
		self.probeTipDiameter = .35
		self.probeCenterX = 0
		self.probeCenterY = 0

		self.rho=np.array([])
		self.Nrestarts = 0
		self.timeCount = 0
        
		self.generateReinjectionProbsYN = False
		self.setProbsMethod = 'by length'  
		self.generatedProbs = []
		self.rollIndex =0
		
		self.setBoundaryLength(3.)
		self.updateBody()
		self.calcArea()
		self.setFieldResolution(.1)
		self.setNpart(self.simPartDensity, self.Area)
		self.updateNparticles()
		self.findOverlaps()

		self.makeMovie = False
		self.countsPerSnapshot = 1000
		self.frameNum = 0 
		
	def updateBody(self):
		#builds the border and defines the edge properties
		xw = self.boxL/2.
		yw = self.boxL/2.
		
        # left border
		self.borderX = np.array([-xw,-xw])
		self.borderY = np.array([-yw,yw])
        
		# create binned drain edge (top)
		BinLength = self.boxL/self.NumberDrainBins
		for bin in np.arange(1, self.NumberDrainBins):
 		    self.borderX = np.append(self.borderX, self.borderX[bin]+BinLength )
 		    self.borderY = np.append(self.borderY, yw)
        
        # right and bottom borders
		self.borderX = np.append(self.borderX, [[xw,xw,-xw]] )
		self.borderY = np.append(self.borderY, [[yw,-yw,-yw]] )
        
		s = 2#'source'
		d = 1#'drain'
		m = 0#'mirror'
		r =-1#'rough'
		f =-2#'false'
		
		self.edgeStyle = [m]
		for bin in np.arange(self.NumberDrainBins):
 		    self.edgeStyle += [d]
                          
		self.edgeStyle += [m, s]
		self.nDeviceEdges=len(self.edgeStyle)
		
		self.setUpEdges()
        
		if self.tip:
			self.addTip()
              
		
	def addTip(self):
        
		# add circular probe tip

		s = 2#'source'
		d = 1#'drain'
		m = 0#'mirror'
		r =-1#'rough'
		f =-2#'false'
        
		thetas = np.linspace(np.pi, 3*np.pi, 50)
		self.borderX = np.append( self.borderX, self.probeCenterX + np.cos(thetas)*self.probeTipDiameter/2)
		self.borderX = np.append( self.borderX, self.borderX[0])
		self.borderY = np.append( self.borderY, self.probeCenterY + np.sin(thetas)*self.probeTipDiameter/2)
		self.borderY = np.append( self.borderY, self.borderY[0])
		
		# add edgeStyles for edges that form the probe
		self.edgeStyle = np.append(self.edgeStyle, f)
		for angle in thetas:
		    if angle == 3*np.pi:  # don't add edge for last point on circle
		        break
		    self.edgeStyle = np.append(self.edgeStyle, m)

		self.edgeStyle = np.append(self.edgeStyle, f)  # add last fictisous edge
        

	def loadBody(self,borderX_in,borderY_in,edgeStyle_in):
		#allows arbitrary boundry 
		self.borderX = borderX_in
		self.borderY = borderY_in
		self.edgeStyle = edgeStyle_in
		
		self.setUpEdges()
		
	def setUpEdges(self):
		#updates parameters based on the borders and edge properties
		self.edgeStyle = np.array(self.edgeStyle)
		Ncorners = len(self.edgeStyle)
		
		#calculate norms and edge lengths
		self.normX,self.normY = polyNorms(self.borderX,self.borderY)
		self.lengths = polyLengths(self.borderX,self.borderY)

		
		s = 2#'source'
		d = 1#'drain'
		m = 0#'mirror'
		r =-1#'rough'
        
		self.sources = [i for i,x in enumerate(self.edgeStyle) if x == s]
		self.drains = [i for i,x in enumerate(self.edgeStyle) if x == d]
		self.Probs = np.zeros(Ncorners)

		if self.setProbsMethod=='by simulation':
        # sets probabilites from a specified list generated from a simulation
			sourcesGeneratedProbSum = np.sum(self.generatedProbs[self.sources])
			drainsGeneratedProbSum = np.sum(self.generatedProbs[self.drains])
            
            #calculate probability of injecting from a source edge
			self.sourceProbs = (self.generatedProbs[self.sources]*self.sourceDrainRatio/
                        (drainsGeneratedProbSum+sourcesGeneratedProbSum*self.sourceDrainRatio))

            #calculate probability of injecting from a drain edge
			self.drainProbs = (self.generatedProbs[self.drains]/
                        (drainsGeneratedProbSum+sourcesGeneratedProbSum*self.sourceDrainRatio))
        
		else:
		#find lengths of source and drain edges to set probabilities
			sourceLength = np.sum(self.lengths[self.sources])
			drainLength = np.sum(self.lengths[self.drains])

            #calculate probability of injecting from a source edge
			self.sourceProbs = (self.lengths[self.sources]*self.sourceDrainRatio/
                        (drainLength+sourceLength*self.sourceDrainRatio))

            #calculate probability of injecting from a drain edge
			self.drainProbs = (self.lengths[self.drains]/
                        (drainLength+sourceLength*self.sourceDrainRatio))

		self.Probs[self.sources] = self.sourceProbs
		self.Probs[self.drains] = self.drainProbs

		self.ProbIdx = (np.where(self.Probs>0))[0]
		self.cumProb = np.cumsum(self.Probs[self.ProbIdx])
        
		
		#initialize Nabsorbed and Ninjected for counting statistics
		self.Nabsorbed = np.zeros(Ncorners)
		self.Ninjected = np.zeros(Ncorners)
		self.symmetrizedNinjected_totalOLD = np.zeros(self.nDeviceEdges)
		self.symmetrizedNinjected_NEW = np.zeros(self.nDeviceEdges)
		self.symmetrizedNinjected_OLD = np.zeros(self.nDeviceEdges)
		self.symmetrizedNinjected_DIF = np.zeros(self.nDeviceEdges)
        
		nContacts = len(self.sources)+len(self.drains) #number of ohmic contacts, M
		self.currentMatrix = np.zeros((nContacts, nContacts)) # MxM matrix between pairs of ohmic contacts; rows=inj; columns=abs
        
		self.contactLookUp = np.sort(np.concatenate([self.sources, self.drains])) # sorted list of ohmic contact indicies
		
		#sets border path for edge crossing functions
		self.borderPath = mpltPath.Path(
									np.vstack((self.borderX,self.borderY)).T)    
    
	def updatateCurrentMatrix(self, partIdx, absIdx, injIdx):
        # tabulates the (prevInjection, newAbsorbtion) statistics
        
		prevInjIdx = self.injectionLookup[partIdx] # lookup previous injection index
        
		if prevInjIdx==-1: # if ptc's first reinjection:
			pass
		else:
			prevInjContact = np.where(self.contactLookUp==prevInjIdx) # prev reinjection contant number
			absContact = np.where(self.contactLookUp==absIdx) # new absorbtion contact number
			self.currentMatrix[prevInjContact, absContact] += 1 # populate the (prevInjection, newAbsorbtion) element
            
		self.injectionLookup[partIdx]=injIdx # set injection lookup to new injection idex
            
	def diffusiveWalls(self):
        # changes diffusive edges into mirror edges
		s,d,m,r,f = 2,1,0,-1,-2 #'source','drain','mirror','rough'
        
		self.edgeStyle = [r if x==m else x for x in self.edgeStyle]
		self.setUpEdges()
        
	def mirrorWalls(self):
        # chanegs mirror edges into diffusive edgses
		s,d,m,r,f = 2,1,0,-1,-2 #'source','drain','mirror','rough'
		
		self.edgeStyle = [m if x==r else x for x in self.edgeStyle]
		self.setUpEdges()            

	def calcArea(self):
        # calculates area of the device
		f =-2 #'false'
		i_f=np.where(self.edgeStyle==f)[0]+1
		borderX_split =np.split(self.borderX,i_f)
		borderY_split =np.split(self.borderY,i_f)

        # trapezoid formula to find area
        #note: device region and cutouts need opposite orientation (ie device cw, cutouts ccw)
		A = 0
		for i,xComponent in enumerate(borderX_split):
		    yComponent = borderY_split[i]
		    if i ==0: # inside device
		        for j in range(len(xComponent)):
		            if j==len(xComponent)-1:
		                A += xComponent[j]*yComponent[0]-xComponent[0]*yComponent[j]
		            else:
		                A += xComponent[j]*yComponent[j+1]-xComponent[j+1]*yComponent[j]


		    else: # subtract off cutouts or probe tip
		        for j in range(len(xComponent)):
		            yPoint = borderY_split[i][j]
		            if j==len(xComponent)-1:
		                A += xComponent[j]*yComponent[0]-xComponent[0]*yComponent[j]
		            else:
		                A += xComponent[j]*yComponent[j+1]-xComponent[j+1]*yComponent[j]
		self.Area = np.abs(A/2)
        
		self.borderX_split = borderX_split
		self.borderY_split = borderY_split
       
	def setSimPartDensity_byLee(self, expPartDensity, Temp):
        # simulation density set to give correct e-e scattering length
		self.simPartDensity = 72.64* Temp**.7 / expPartDensity**.1
		
	def setNpart(self, simPartDensity, Area):
        # ptc. density and area set the number of particles
		self.Npart = int(simPartDensity*Area)
		self.injectionLookup = -1*np.ones(self.Npart) # prev reijenction indices for all ptcs; start with -1
    
	def setEmin(self, expPartDensity, Temp, simPartDensity=20):
        # this function is for a "sweep T, for fixed n" series
        # sim and exp densities and temp set scattering length via Emin
		self.Emin = 1 - 78312.7/((np.sqrt(expPartDensity)*simPartDensity)/Temp**1.5)**1.25
    
	def setEmin_byFermi(self, expPartDensity, Temp):
        # min scattering energy to get the correct width of the fermi distribution (kbT/Ef)
		self.Emin = 1 - 369.3* Temp /expPartDensity**.5
    
	def calc_pScatter(self):
		self.p_scatter = 0.000024791*self.Temp
    
	def setSourceDrainRatio(self,Ratio_in):
		#sets relative probability of injection from the source and drain edges
		self.sourceDrainRatio = Ratio_in
		
	def injectPosition(self):
		#first calculates reinjection edge based on set probabilities and then calculates the reinject position
		injectIdx = self.ProbIdx[np.where(self.cumProb>np.random.rand())[0][0]]
		injectFraction = np.random.rand()
		x1 = self.borderX[injectIdx+1]
		x0 = self.borderX[injectIdx]
		y1 = self.borderY[injectIdx+1]
		y0 = self.borderY[injectIdx]
		
		injectX = x0+injectFraction*(x1-x0)
		injectY = y0+injectFraction*(y1-y0)
		
		return injectX,injectY,injectIdx

	def injectFromContact(self,injectIdx):
		#calculates position to inject along a specific contact
		injectFraction = np.random.rand()
		x1 = self.borderX[injectIdx+1]
		x0 = self.borderX[injectIdx]
		y1 = self.borderY[injectIdx+1]
		y0 = self.borderY[injectIdx]
		
		injectX = x0+injectFraction*(x1-x0)
		injectY = y0+injectFraction*(y1-y0)
		
		return injectX,injectY

				
	def findNewInject(self):
		#choose reinjection edge and searches for a pos to inject without overlap with existing particles
		injectX,injectY,injectIdx = self.injectPosition()
		i =0
		while self.getRegionOccupied(injectX,injectY):
			injectX,injectY = self.injectFromContact(injectIdx)
			i=+1
			if i==5:
				break
		return injectX,injectY,injectIdx
		
	def setOverlapRadius(self, R_in):
		#sets radius of particle interactions
		self.rOverlap = R_in
		
	def updateNparticles(self):
		#updates the number of simultaneous particles in the simulation
		Npart_in = self.Npart
		
		# randomly assign position
		self.Xpos = self.boxL*(np.random.rand(Npart_in)-0.5)
		self.Ypos = self.boxL*(np.random.rand(Npart_in)-0.5)
		
		# find which particles are out of bounds
		outOfBounds = np.invert(self.borderPath.contains_points(
											np.vstack((self.Xpos,self.Ypos)).T))
		
		#keep randomly searching until all particles are in bounds
		while any(outOfBounds):
			for idx in np.where(outOfBounds)[0].tolist():
				self.Xpos[idx] = self.boxL*(np.random.rand()-0.5)
				self.Ypos[idx] = self.boxL*(np.random.rand()-0.5)
			outOfBounds = np.invert(self.borderPath.contains_points(
											np.vstack((self.Xpos,self.Ypos)).T))
		self.findOverlaps()

		#assign random velocity direction
		thetas = np.random.rand(Npart_in)*2.*np.pi
		self.vX = np.cos(thetas)
		self.vY = np.sin(thetas)
		#assign fixed momentum amplitude
		self.pR = np.ones(np.shape(thetas))
		
		#intializes particular statistics 
		self.L_NoScatter = np.zeros(Npart_in)
		self.traces = [[] for i in range(Npart_in)]
		
		#set up indices for quick retreival 
		self.i_lookup = np.zeros(int(self.Npart*(self.Npart-1)/2),dtype = int)
		self.j_lookup = np.zeros(int(self.Npart*(self.Npart-1)/2),dtype = int)
		_idx = 0
		for i in range(self.Npart):
			for j in range(i+1,self.Npart):
				self.i_lookup[_idx] = i
				self.j_lookup[_idx] = j
				_idx += 1

		
	def findOverlaps(self):
		#finds particle pairs that are within a the rOverlap radius. Successive-
		#ly checks the separation in X and then Y before calculating the cartes-
		#ian distance to save computational resources
            
		self.overlaps = squareform(self.checkOverlaps())
	
	def checkOverlaps(self):
		r = np.vstack((self.Xpos,self.Ypos)).T
		return pdist(r) < self.rOverlap

	def getRegionOccupied(self,X_in,Y_in):
		#checks if a given X_in and Y_in are in an already occupied region
		occupied = False
		for i in range(self.Npart):
			absDx = np.abs(self.Xpos[i]-X_in)
			if absDx < self.rOverlap:
				absDy = np.abs(self.Ypos[i]-Y_in)
				if absDy < self.rOverlap:
					if absDx**2+absDy**2 < self.rOverlap**2:
						occupied = True
		return occupied
		
	
	def randomPointInBounds(self):
		#finds a single coordinate-pair that is inbounds
		x,y = self.boxL*(np.random.rand()-0.5),self.boxL*(np.random.rand()-0.5)
		inBounds = self.borderPath.contains_points(np.vstack((x,y)).T)
		while not inBounds:
			x = self.boxL*(np.random.rand()-0.5)
			y = self.boxL*(np.random.rand()-0.5)
			inBounds = self.borderPath.contains_points(np.vstack((x,y)).T)
		return x,y
		
		
	def updateScatterProb(self,p_in):
		#sets the probability of scattering in the bulk
		self.p_scatter = p_in
		
	def setStepSize(self,dX):
		#sets the distance traversed in a single timestep
		self.DeltaX = dX
		
	def setBoundaryLength(self,L):
		#sets the overall boundary of the system
		self.boxL = L
		self.boxRange = [[-L/2.,L/2.],[-L/2.,L/2.]]
	
	def setFieldResolution(self,dX):
		self.fieldResolution = dX
		#sets several spatial histograms with box size "dX"
		Nr = int(np.round(self.boxL/dX))
		
		self.Nx = Nr
		self.Ny = Nr

		self.rho = np.zeros((self.Nx,self.Ny))  # particle density
		self.Px = np.zeros((self.Nx,self.Ny))   # X velocity density
		self.Py = np.zeros((self.Nx,self.Ny))   # Y velocity density
		self.Erho = np.zeros((self.Nx,self.Ny)) # Energy Density
        
        # density map for each source/drain
		self.transportMap = np.array([np.zeros((self.Nx,self.Ny)) for contact in self.contactLookUp])
		
		#sets the X,Y mesh
		_,self.histX,self.histY = np.histogram2d(
							np.zeros(self.Nx),np.zeros(self.Ny),bins=[self.Nx,self.Ny],
							range = self.boxRange)
		
	def randCos(self):
		#returns an angle between -pi/2 and pi/2 following a cosine distribution
		theta = np.pi*(np.random.rand()-.5)
		while np.random.rand()>np.cos(theta):
			theta = np.pi*(np.random.rand()-.5)
		return theta
		
	def randCosNorm(self,normX,normY):
		#returns a vector selected from a cosine distribution relative to an
		#input normalvector.
		theta = self.randCos()
		vXNew = -normX*np.cos(theta)-normY*np.sin(theta)
		vYNew = -normY*np.cos(theta)+normX*np.sin(theta)
		return vXNew, vYNew

	
	def timestep(self):
		# executes a single timestep of the simulation
		
		#propagate all particles in the direction of their velocity vector
		self.Xpos += self.vX*self.DeltaX
		self.Ypos += self.vY*self.DeltaX
		#update mean free path statistics
		self.L_NoScatter+=self.DeltaX
		
		#find which particles ended up out of bounds
		outOfBounds = np.invert(
			self.borderPath.contains_points(np.vstack((self.Xpos,self.Ypos)).T))

		injectedPtcs = []
        
		for i in range(self.Npart):
			
			if outOfBounds[i]:
				# backtracks if out of bounds
				X0 = self.Xpos[i] - self.vX[i]*self.DeltaX
				Y0 = self.Ypos[i] - self.vY[i]*self.DeltaX

				while not self.borderPath.contains_point((X0,Y0)):
					print("X0, Y0 out of bounds")
					theta = np.random.rand()*2*np.pi
					X0 = self.Xpos[i] - np.sin(theta)*self.DeltaX
					Y0 = self.Ypos[i] - np.cos(theta)*self.DeltaX
				
				#finds index of the line(s) crossed
				lineCrossed = segCrossPoly(self.borderX,self.borderY,
								np.array((X0,Y0)),
								np.array((self.Xpos[i],self.Ypos[i])))
                
				j = 0
				for line in lineCrossed:
					if self.edgeStyle[line] == -2:
						lineCrossed = np.delete(lineCrossed, j)
						j -= 1
					j += 1
                
				Xf, Yf = self.Xpos[i],self.Ypos[i]
				self.Xpos[i],self.Ypos[i] = X0, Y0
				
				lenCrossed = len(lineCrossed)
				if lenCrossed == 0:  # no segments crossed; exactly betweeen two
					self.vX[i] *= -1
					self.vY[i] *= -1
					
                    #print('ptc. crossed between two edges')
				
				elif lenCrossed == 1:  # 1 segment crossed; the expected case
					#select the first in the list of lines crossed (convenience)
					idxCross = lineCrossed[0]
					if self.edgeStyle[idxCross] == -1:  #diffusive edge
						self.scatterFromDiffusiveEdge(i,idxCross)
						
					if self.edgeStyle[idxCross] == 0:	#mirror edge
						self.reflectFromMirrorEdge(i,idxCross)
						
					if self.edgeStyle[idxCross] > 0:		#source or drain
						#re-inject particle and take statistics
						self.consumeAndReinject(i,idxCross)
						injectedPtcs.append(i)

				else:
                    
					self.NcornersCrossed += 1
					if self.NcornersCrossed%100 == 0: #every 100 corner crossings,
						print("Total Corners Crossed:", self.NcornersCrossed)
					crossedSourceDrains = np.where(self.edgeStyle[lineCrossed]>0)[0]
                    
					if len(crossedSourceDrains)>0: #if crossed a source/drain
						#consume by random selection among the source and drain
						#edges in the lineCrossed list
						randomCrossedSD = np.random.randint(len(crossedSourceDrains))
						self.consumeAndReinject(
											i,crossedSourceDrains[randomCrossedSD])
						injectedPtcs.append(i)
                         
					else: 
                        
                        #if crossed only other edge types
						#scatter from an edge whose norm has a component in the
						#direction of the velocity
						nonSourceDrainEdge = 0
						idxCross = lineCrossed[nonSourceDrainEdge]
						while (self.vX[i]*self.normX[idxCross] +
										self.vY[i]*self.normY[idxCross] < 0):
							nonSourceDrainEdge += 1
							idxCross = lineCrossed[nonSourceDrainEdge]
                        
						if self.edgeStyle[idxCross] == -1:  #diffusive edge
							self.scatterFromDiffusiveEdge(i,idxCross)
						
						if self.edgeStyle[idxCross] == 0:	#mirror edge
							self.reflectFromMirrorEdge(i,idxCross)
                            
						
			elif self.p_scatter > np.random.rand(1):
				#randomize direction in bulk
				theta = np.random.rand(1)*2.*np.pi
	
				self.vX[i] = np.cos(theta)
				self.vY[i] = np.sin(theta)
		
		_overlaps = self.checkOverlaps()
        
        # dissallow reinjected particles from scattering for 1 timestep
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
				
					#measure momentum amplitude
					self.pR[i] = np.sqrt(np.sum(p3**2))
					self.pR[j] = np.sqrt(np.sum(p4**2))
					#set velocities
					self.vX[i],self.vY[i] = p3/self.pR[i]
					self.vX[j],self.vY[j] = p4/self.pR[j]
					
					self.L_NoScatter[i]=0
					self.L_NoScatter[j]=0
				
			self.overlaps[i,j] = True
			
		for i in injectedPtcs:
			self.overlaps[:,i] = False
			self.overlaps[i,:] = False
            
		#update global time index
		self.timeCount+=1

		if self.timeCount == self.initializationCount:
			#reset statistics after initalizationCount number of timesteps
			self.Px*=0
			self.Py*=0
			self.rho*=0
			self.Erho*=0
			self.Ninjected*=0
			self.Nabsorbed*=0
			self.currentMatrix*=0
			self.transportMap*=0

    
		if np.mod(self.timeCount,5)==0:
			#update histograms every 5 timesteps
			h,_,_ = np.histogram2d(self.Xpos,self.Ypos,bins=[self.Nx,self.Ny],
										range = self.boxRange,weights = self.vX)
			self.Px+=h
			
			h,_,_ = np.histogram2d(self.Xpos,self.Ypos,bins=[self.Nx,self.Ny],
										range = self.boxRange,weights = self.vY)
			self.Py+=h
			
			h,_,_ = np.histogram2d(self.Xpos,self.Ypos,bins=[self.Nx,self.Ny],
										range = self.boxRange,weights = self.pR)
			self.Erho+=h
			
			h,_,_ = np.histogram2d(self.Xpos,self.Ypos,bins=[self.Nx,self.Ny],
										range = self.boxRange)
			self.rho+=h
            
            # update particle density maps for each ohmic contact; only reinjected particles from that contact contribute 
			for i,contact in enumerate(self.contactLookUp):
				h,_,_ = np.histogram2d(self.Xpos[self.injectionLookup==contact],
                                       self.Ypos[self.injectionLookup==contact],
                                       bins=[self.Nx,self.Ny], range = self.boxRange)
				self.transportMap[i]+=h
                
		
				
		
	def scatterFromDiffusiveEdge(self,partIdx,edgeIdx):
		self.vX[partIdx],self.vY[partIdx] = self.randCosNorm(
										self.normX[edgeIdx],self.normY[edgeIdx])
						
	def reflectFromMirrorEdge(self,partIdx,edgeIdx):
		self.vX[partIdx],self.vY[partIdx] = mirrorNorm(self.normX[edgeIdx],
						self.normY[edgeIdx],self.vX[partIdx],self.vY[partIdx])
        
	def findNewInjectAlongContact(self,edgeIdx):
	    injectX,injectY = self.injectFromContact(edgeIdx)
	    i =0
	    while self.getRegionOccupied(injectX,injectY):
	        injectX,injectY = self.injectFromContact(edgeIdx)
	        i=+1
	        if i==5:
	            break
	    return injectX, injectY
	
	def reinjectFromSameContact(self, edgeIdx):
        # reinject from same contact with probability (2-SD); otherwise reiject from random source contact
	    if np.random.rand() > self.sourceDrainRatio-1:
	        injectX, injectY = self.findNewInjectAlongContact(edgeIdx)
	        injectIdx = edgeIdx
	    else:
	        totalSourceProb = np.sum(self.sourceProbs)
	        sourceCumProb = np.cumsum(self.sourceProbs)/totalSourceProb
	        injectIdx = self.sources[np.where(sourceCumProb>np.random.rand())[0][0]]
	        injectX, injectY = self.findNewInjectAlongContact(injectIdx)
	    return injectX,injectY, injectIdx
    
	def consumeAndReinject(self,partIdx,edgeIdx):
        # checks if contact will transmit an incoming electron
        # if z-zero current sim, reinject from same contact
        # otherwise, reinject based on set probabilites; update stats
        
		if np.random.rand() > self.p_transmission:
			self.scatterOffContact(partIdx, edgeIdx)
            
		else:
			if self.generateReinjectionProbsYN or self.sourceDrainRatio==1.:
				xNew,yNew = self.injectFromContact(edgeIdx)
				idxNew=edgeIdx
			else:
				xNew,yNew,idxNew = self.injectPosition()
			self.Nabsorbed[edgeIdx]+=1		
			self.Ninjected[idxNew]+=1

			self.updatateCurrentMatrix(partIdx, edgeIdx, idxNew) # particle number, absorbtion index, new injection index

			vXNew,vYNew = self.randCosNorm(
                            self.normX[idxNew],self.normY[idxNew])

			self.vX[partIdx],self.vY[partIdx] = vXNew,vYNew

            #small offset to make sure particle is in bounds
			self.Xpos[partIdx] = xNew+vXNew*self.DeltaX*.0001
			self.Ypos[partIdx] = yNew+vYNew*self.DeltaX*.0001
	
    
	def scatterOffContact(self, partIdx, contactIdx):
        
        # particles not transmitted through a contact; scatters off
        # mirror/diffusive chosen from non contact edge styles
        
		nonContantEdges = self.edgeStyle[(self.edgeStyle==-1)+(self.edgeStyle==0)]
		i = np.random.randint(len(nonContantEdges))
        
		if nonContantEdges[i] == -1: # diffusive edge
			self.scatterFromDiffusiveEdge(partIdx,contactIdx)
            
		elif nonContantEdges[i] == 0: # mirror edge
			self.reflectFromMirrorEdge(partIdx,contactIdx)

	def consumeAndReinject_withE1(self,partIdx,edgeIdx):
        
        # reinjects particles with an energy of 1
        # replaces the original consume and reinejct function when called
    
		if np.random.rand() > self.p_transmission:
			self.scatterOffContact(partIdx, edgeIdx)
            
		else:
			if self.generateReinjectionProbsYN or self.sourceDrainRatio==1.:
				xNew,yNew = self.injectFromContact(edgeIdx)
				idxNew=edgeIdx
			else:
				xNew,yNew,idxNew = self.injectPosition()
			self.Nabsorbed[edgeIdx]+=1		
			self.Ninjected[idxNew]+=1

			self.updatateCurrentMatrix(partIdx, edgeIdx, idxNew) # particle number, absorbtion index, new injection index

			vXNew,vYNew = self.randCosNorm(
                            self.normX[idxNew],self.normY[idxNew])

			self.vX[partIdx],self.vY[partIdx] = vXNew,vYNew
			self.pR[partIdx] = 1

            #small offset to make sure particle is in bounds
			self.Xpos[partIdx] = xNew+vXNew*self.DeltaX*.0001
			self.Ypos[partIdx] = yNew+vYNew*self.DeltaX*.0001
        
        
	def saveFrame(self):
        
        # save a snapshot of the simulation

		if self.frameNum ==0:
			if not os.path.isdir(self.saveLoc+'/movieFrames'):
				os.mkdir(self.saveLoc+'/movieFrames')
		file = self.saveLoc+'/movieFrames/frame'+("_%03d"%self.frameNum)+'.jpg'
		totalAbs = sum(self.Nabsorbed)


		fig,ax= plt.subplots(1,1, figsize=(10,8))
		for i,xComponent in enumerate(self.borderX_split):
			yComponent= self.borderY_split[i]
			ax.plot(xComponent, yComponent,'black')
		ax.plot(self.Xpos[0], self.Ypos[0], 'b.', markersize=12)
		ax.plot(self.Xpos[1:], self.Ypos[1:], 'r.') 

		ax.text(self.boxL/5,self.boxL, 'total absorbed: '+("%03d"%totalAbs), fontsize=14)
		ax.text(self.boxL/5,self.boxL+1.25, 'timesteps: '+("%03d"%self.timeCount), fontsize=14)
		ax.text(-self.boxL/2,self.boxL, 'temperature: '+("%03d"%self.Temp), fontsize=14)
		ax.text(-self.boxL/2,self.boxL+1.25, 'sim density: '+("%03d"%self.simPartDensity), fontsize=14)
		ax.axis('off')
		fig.savefig(file)
		plt.close()
		self.frameNum+=1
    
	def saveState(self,fname):
		#saves all object data into an npz file
		varNames = list(self.__dict__.keys())
		#varNames.remove('i_lookup')
		#varNames.remove('j_lookup')
		#varNames.remove('overlaps')
		#print(varNames)

		
		saveFunction = "np.savez(fname, "
		for name in varNames:
			saveFunction = saveFunction+name+' = self.'+name+', '
			
		saveFunction = saveFunction[:-2] + ')'
		exec(saveFunction)
        
	def saveReinjectionStats(self, fname):
        # saves the sum of reinjections for a group of simulations per edge into a csv file
        
		splitName= fname.split('/')
		fnum = ((splitName[-1].split('_'))[-1]).split('.')[0]
		del splitName[-1]
		fnameDir='/'.join(splitName)
        
		Ninjected_DeviceEdgesOnly = np.roll(self.Ninjected[0:self.nDeviceEdges],self.rollIndex)
                # NOTE: this only applies to symmetric devices;
                # the roll index is unique for every geometry and can be specified in that device's init function
		#self.symmetrizedNinjected_NEW = (Ninjected_DeviceEdgesOnly[::-1]+Ninjected_DeviceEdgesOnly)/2
		self.symmetrizedNinjected_NEW = Ninjected_DeviceEdgesOnly-self.symmetrizedNinjected_totalOLD
		self.symmetrizedNinjected_DIF = self.symmetrizedNinjected_NEW-self.symmetrizedNinjected_OLD
        
		file = fnameDir+'/injectionStats'+fnum+'.csv'
		with open(file, 'a') as csvFile:
			writer=csv.writer(csvFile)
			writer.writerows([(self.symmetrizedNinjected_NEW).tolist()])
        
		if np.mean(self.symmetrizedNinjected_DIF)<1000:
			if self.timeCount > self.initializationCount:
				print(f'Equilibrium reached after {self.timeCount} time steps')
                
		self.symmetrizedNinjected_totalOLD =Ninjected_DeviceEdgesOnly
		self.symmetrizedNinjected_OLD = self.symmetrizedNinjected_NEW
    
	def runAndSave(self,Nsteps,dNsave,fname):
		
		t0 = time.time()
		self.saveState(fname)

		for i in range(Nsteps):
			
			self.timestep()
            
			if self.makeMovie:
				if self.timeCount%(dNsave/self.countsPerSnapshot) == 0:
 					self.saveFrame()
            
			if self.timeCount%dNsave == 0:
				t1 = time.time()
				print('%d:  %g'%(self.timeCount,(t1-t0)))
				t0 = t1
				self.saveState(fname)
            
                # store the reinjection statistics in a .csv file 
				if self.timeCount%(self.initializationCount/10)==0:
 					if self.generateReinjectionProbsYN:
                        # store the reinjection statistics once per simulation batch
 						self.saveReinjectionStats(fname) 
