"""
Author: Scott Jensen and Daniel Hartzler
This code is not optimized (though it is of one time use so it probably wont ever be).
This library will take data the detector image and calibrate 
"""

from __future__ import division
import numpy as np

def flip_spectral_direction(flip):
	"""
	Changes the output direction of the spectrum (based on orientation of detector)
	returns if the spectrum direction should flip with a negative value
	"""
	if flip ==True:
		return(-1)
	else:
		return(1)

def singlePointCalib (spectraLen, calibPoint, calibE=6493.3, \
	EoCrystal=6457, Radius=500, pixelSize=0.05, flip=True):

	"""
	#################################################
	#Makeing a spectrum from pixel distances
	#################################################

	Calibrates a spectra from a single point for von Hamos analyzer (namely MnCl2 as default)
	This takes the known energy to position function and fits to the a single calibrated point
	
	Pixel Size is in mm
	
	#Psuedo code
	Use calibration file to find Kb1,3
	from Kb1,3 space out energy of each point by using the function
	E=Eo/sin(theta)
	or
	E=E0/(2R/((2R)^2+x^2)^0.5)
	where x is the distance from the normal of the crystal lattice or detector face
	R is the radius of curvature of the spectrometer crystal (von Hamos curvature)
	
	So now I need to assign an energy to one point on the spectra and take 
	steps of size x along the path assigning the energy as I go. First I need
	the absolute x
	
	x=2R((E/Eo)^2-1)^0.5
		
	Define point in spectra with known energy to find an absolute x
	Use the absolute x to calibrate whole spectra
	
	x = xcalib +/- (ni -ncalib)*pixelSize
	
	where n corresponds to the integer step in the spectra
	"""
	flipC = flip_spectral_direction(flip)
	xDetCalib=2*Radius*((calibE/EoCrystal)**2-1)**0.5 #Calibration point
	
	newEList=[]
	for currentSpot in range(spectraLen): #New Scale from detector pixel
		xCur=xDetCalib-flipC*(calibPoint-currentSpot)*pixelSize #Current x location
		newE=EoCrystal*((xCur/(2*Radius))**2+1)**0.5 #New energy for each point in spectra
		newEList.append(newE)
		
	return(np.asarray(newEList))

def replacePixels (pixel_map, replaceme=10**20):
	"""
	Replace bad Pixels with only one bad pixel next to it.
	Finds all pixels with replaceme values and interpolates if there is only one bad value around it, then
	it will use the 8 pixels around to create an average for the marked pixel
	If this approach is too slow, numpy would likely reduce process time.
	pixel_map: array of values representing the detector image
	replaceme: the values listed as "bad pixels"
	"""
	
	pixMask=np.ones((3,3)) #Create a mask to unselect the center pixel for a sum
	pixMask[1,1]=0
	
	pixMapTemp=np.zeros(pixel_map.shape) #fill in new array (don't want to overwrite the original
	for i in range(len(pixel_map[:,0])):
		for j in range(len(pixel_map[0,:])):
			if is_surrounded_by_good_pixels(pixel_map, replaceme,i,j):
				pixMapTemp[i,j] = replace_bad_with_average(pixMapTemp,pixel_map, pixMask, i, j)
			elif is_next_to_more_bad_pixels(pixel_map,replaceme,i,j):
				pixMapTemp[i,j]=replaceme #Leave as bad
			else: # No bad pixels
				pixMapTemp[i,j]=pixel_map[i,j]
	return(pixMapTemp)

def is_surrounded_by_good_pixels(pixel_map, replaceme,i,j):
	"""Makes sure there are no other values around with a value of "replaceme" nearby
	pixel_map: pixelated image
	replacemen: value for bad pixels
	i, j: indecies for the pixel of interest in the map"""
	return(pixel_map[i,j]==replaceme and (pixel_map[i-1:i+2,j-1:j+2]==replaceme).sum()==1)

def replace_bad_with_average(pixel_map, pixMask, i, j):
	"""Replace the bad pixel with a sum of all points around it
	pixel_map: pixelated image
	pixMask: Mask to remove center pixel for average
	i, j: indecies for the pixel of interest in the map
	"""
	valOut=(np.sum(pixel_map[i-1:i+2,j-1:j+2]*pixMask))/8.

	#check to see if any 0's are used in the sum (indication
	#that a larger range of pixels needs to be used in the first run
	is_next_to_edge(pixel_map, i, j)
	return(valOut)

def is_next_to_edge(pixel_map, i, j):
	"""Determins if the i,j position is too close to the edge"""
	if (pixel_map[i-1:i+2,j-1:j+2]==0).sum()>=1:
		print('0 value used in sum for above location, likely too close to edge')

def is_next_to_more_bad_pixels(pixel_map,replaceme,i,j):
	"""is the value at i,j surrounded by points which are bad pixels?"""
	return(pixel_map[i,j]==replaceme and (pixel_map[i-1:i+2,j-1:j+2]==replaceme).sum()>=1)

def maskMaker (arrayshape, lineEdges):
	"""Makes a mask of 1 between lines and 0's outsize"""
	#Make a mask of 0's
	mask=np.zeros((arrayshape))

	#Fill in 1's between the range of values
	for i, edge in enumerate(lineEdges):
		if not edge[1]-edge[0]==0:
			mask[i,edge[0]:edge[1]]=np.ones(edge[1]-edge[0])
	return mask

def lineJEdges(dim, lineEnd1_xy,lineEnd2_xy,lineWidth):
	"""
	Generates a list of j-edge points for a user defined line... follows the direction of the spectra.
	This is for rotation of the spectrum which is not parallel to the edges of the detector.
	dim: dimension of array
	lineEnd1_xy, lineEnd2_xy: endpoints of the line that is parallel to 
	lineWidth: number of pixels wide

	Note: Also should be implimented in Numpy
	"""
	[j11,i11],[j22,i22]=lineEnd1_xy,lineEnd2_xy
	slope = ((j22-j11)/(i22-i11)) ## slope and intercept of line center
	intercept = (j22-slope*i22)

	j_upper = lambda x:slope*x+intercept+lineWidth/2	## upper line edge
	j_lower = lambda x:slope*x+intercept-lineWidth/2	## lower line edge

	lineEdge = []
	for i in range(dim[0]):
		j1 = j_lower(i)	 ## for a given 'i', between what two 'j' values is the "line" defined?
		j2 = j_upper(i)
		if j1<0:	## edge checking
			j1=0
		if j2<0:
			j2=0
		if j1>dim[1]-1:	 ## edge checking
			j1=dim[1]-1
		if j2>dim[1]-1:
			j2=dim[1]-1
		j11 = min([j1,j2])
		j22 = max([j1,j2])
		lineEdge.append([int(j11),int(j22)]) ## line definition array ('j' endpoints for each value of 'i') 

	return lineEdge, [slope, intercept+lineWidth/2], [slope,intercept-lineWidth/2]

def interpolate_2D(array,i,j):
	"""Interpolates to a point in the array"""
	dim = array.shape
	iDn = int(i)
	jDn = int(j)
	iUp = int(i+1)
	jUp = int(j+1)
	
	goodPoints = [0,0,0,0] #Used to mark in bounds and out of bounds points
	#if any point is out of bounds set a goodPoints flag
	if iDn<0 or iDn >= dim[0]:
		goodPoints[0] = 0
	else:
		goodPoints[0] = 1
	if iUp>=dim[0] or iUp < 0:
		goodPoints[1] = 0
	else:
		goodPoints[1] = 1
	if jDn<0 or jDn >= dim[1]:
		goodPoints[2] = 0
	else:
		goodPoints[2] = 1
	if jUp>=dim[1] or jUp<0:
		goodPoints[3] = 0
	else:
		goodPoints[3] = 1

	if np.sum(goodPoints) < 2:
		return 0

	points = [0,0,0,0]
	if goodPoints[0] and goodPoints[2]:
		points[0] = array[iDn][jDn]
	if goodPoints[2] and goodPoints[1]:
		points[1] = array[iUp][jDn]	
	if goodPoints[0] and goodPoints[3]:
		points[2] = array[iDn][jUp]
	if goodPoints[1] and goodPoints[3]:
		points[3] = array[iDn][jUp]	
	#####################################
	## values for detecting bad pixels ##
	## work in progress
	#points = [array[iDn][jDn], array[iUp][jDn], array[iDn][jUp], array[iUp][jUp]]

	####################################
	
	######################
	## 2D interpolation ##
	f_i_j1 = (iUp-i)*(points[0]) + (i-iDn)*(points[1])
	f_i_j2 = (iUp-i)*(points[2]) + (i-iDn)*(points[3])
	f_i_j = ((jUp-j)*f_i_j1 + (j-jDn)*f_i_j2)/np.sum(goodPoints)
	#f_i_j gets returned
	######################
	return f_i_j

def extractSpectrum(array, lineLower, lineWidth,extraSteps = 0):
	"""
	Extracts the spectrum from photon counts per pixel images, given the lower edge and the line width of the spectrum.
	array: detector array with number of photons in each pixel
	lineLower: lower edge of the line parallel to the spectrum recorded at some angle on the detector
	lineWidth: Width of the spectrum to be extracted 
	extraSteps: steps beyond the boundry to use
	"""
	#dim = array.shape
	dim = (array.shape[1],array.shape[0]) ## messed up and swapped up 'i' and 'j' at some point

	#NOTE:  'i' and 'j' are indicies that correspond to physical camera pixels.
	#	   'I' and 'J' are indicies that correspond to a grid (with pixel length steps) that is rotated to match the perp and parallel of the user defined line.
	theta = np.arctan(lineLower[0])
	lineWidthIJ = lineWidth*np.cos(theta) ## linewidth perpendicular to the user defined line. What would be considered the real linewidth (the variable "lineWidth" is width in the 'j' direction only)
	IHat = np.array([np.cos(theta),np.sin(theta)])	#dI = [di,dj]. Unit step in the "I" direction (i.e. along the user defined line)
	JHat = np.array([-np.sin(theta),np.cos(theta)])   #dJ = [di,dj]. Unit step in the "J" direction (i.e. perpendicular to the user definded line)
	dJ = .5	## step in half unit steps in "J" direction

	ij = np.array([-extraSteps,lineLower[1]]) ## (i,j) start point for stepping in the 'I' direction
	Jmax = lineWidthIJ # how far to step in the 'J' direction
	spect = []
	while(ij[0] < dim[0]+extraSteps): #step in 'I' direction
		ij_perp = ij # start point for stepping in the 'J' direction
		_sum = 0
		J = 0
		while(J < Jmax): #step in the 'J' direction (and average all points along this direction)
			point = interpolate_2D(array,ij_perp[0],ij_perp[1]) # get interpolated point. If you want to change the 2D interpolation to something else, change this line (i.e. point = yourFunction())
			_sum += point
			J += dJ
			ij_perp = ij_perp + dJ*JHat
		ij = ij + IHat
	
		#Taking out the division by dJ so that it is independent of line width
		spect.append(_sum*dJ) # Jmax/dJ == N points along 'J' direction

	return np.array(spect)

def specFromImage(detImage, calibPoint = 500.0, invert = True, lineEnds = [[20,526],[727,193]], width = 110, pixelLen = 0.05):
	"""
	Main funciton that will take the detector image, calibrate the pixels to x-ray emission energy.
	The spectral intensity is then extracted as intesity and x-ray energy.
	This function also replaces bad pixels with averages of surrounding pixels.
	detImage: array of values 
	calibPoint: energy of the calibration point
	invert: if the spectral extraction needs flipped based on detector orientation
	lineEnds: endpoints of the spectral line
	width: pixelWidth of spectrum
	pixelLen: dimension of a pixle in mm
	"""
	#Use the outmost range 1 pixel around desired range
	_, _, lineLower2 = lineJEdges(detImage.shape,lineEnds[0],lineEnds[1],width)
	
	#Remove the bad points
	pixOut = replacePixels(detImage)
	
	#Make the spectra 
	spec= extractSpectrum(pixOut, lineLower2, width)

	#Import Energy
	eRange = singlePointCalib(len(spec), calibPoint, flip = invert, pixelSize=pixelLen)

	#Cut out high points (bad pixels)
	spec1 = spec[spec < 10.**16]
	eRange1 = eRange[spec < 10.**16]

	return (eRange1, spec1)

