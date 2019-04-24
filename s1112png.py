#!/usr/bin/env python

import sys
import os
#import scipy.misc.imsave
#import skimage.io.imsave
#depricated #from skimage import io
import imageio as io 
import numpy as np
#import numpy.zeros as npz
#import numpy.uint8 as npui8
import math
import json
import geo
import h5py

#this will convert s111 v 1.0.0 hdf5 data format to png format that can be used to transfer the data within a web browser using a supported json file for the metadata to the PNG.

#TODO: test the size of hdf5 s111 and compare against the png/json combos

def convertS111(s111,bname,s111meta,sfill):
	createPNG(s111,bname,sfill)	
	metadata = createMetadata(s111,bname,s111meta)	
	return metadata
	
#NOTE: if this is generic for all datatypes will have to change 1st parameter name (and make sure accessing values is still possible...or better yet, set up a data structure that contains all needed values so the input variable to the shared functions is standardized.
def createPNG(s111,bname,sfill):			
	speedInd = 0
	dirInd   = 1

	#numpy.zeros(shape, dtype=float, order='C') -> Return a new array of given shape and type, filled with zeros.
	#i assume this is X,Y,4 (or Y,X,4) and dtype uint8 = Byte (-128 to 127)
	row = s111.shape[0]
	col = s111.shape[1]

	dim = 4 #s111 - 4 here for RGBA
	out = np.zeros((row,col,dim),np.uint8)
	
	#what is this doing exactly? converting speed/direction values to u and v values
	#This will loop through each cell in the dataset and take the value that is found
	for i in range(row-1,-1,-1):
		for j in range(0,col):			
		#this checks if the 'mask' bit is on...i.e. if there is a speed/direction value in this location
			if s111[i,j][speedInd] != sfill:#-9999.0 is not the s111 standard for land value; cannot hard code this, must use group_f value in file 
				speed_ms = s111[i,j][speedInd]*0.514
				dir_rad = math.radians(s111[i,j][dirInd])
				fx = speed_ms*math.sin(dir_rad)
				fy = speed_ms*math.cos(dir_rad)
				fxi = int((fx+20.48)/0.01)
				fyi = int((fy+20.48)/0.01)

				#the 0,1,2,3 the RGBA values for the PNG
				out[i,j,0] = (fxi & 0x0ff0)>>4
				out[i,j,1] = ((fxi & 0x000f)<<4)|((fyi & 0x0f00)>>8)
				out[i,j,2] = fyi & 0x00ff				
				#out[i,j,3] = mask[i,j]*255 #this is looking for just a 1 or a 0 but this will only be executed if the "mask" is basically 'false' i.e. there IS a value here so it will always be 255
				out[i,j,3] = 255
			else: 
				s111[i,j][speedInd] = -9999.0
				#print(s111[i,j][speedInd])
	outfn = getOutfn(bname)
	#scipy.misc.imsave(outfn,out)
	out = np.flipud(out) #image was upside-down on save
	io.imsave(outfn,out)
	#this was depricated need to use imageio.imwrite 
	io.imwrite(outfn, out)

#NOTE: if this is generic for all datatypes will have to change 1st parameter
def createMetadata(s111,bname,s111meta):
	#if this "label" item is just for identification then it could really be anything to identify it.
	#metadata = {'label':os.path.splitext(os.path.basename(fname))[0].split('-',1)[1],'sources':[]}
	metadata = {'label':bname,'sources':[]}

	#After PNG is created from the data values, the accompanying JSON file 
	#is created to help the javascript web interface know the metadata details. 
	minLat = np.float64(s111meta.attrs.get('southBoundLatitude'))	#s111meta.southBoundLatitude ->mask.Minimum_Latitude, s111meta.northBoundLatitude ->mask.Maximum_Latitude
	maxLat = np.float64(s111meta.attrs.get('northBoundLatitude'))
	minLon = np.float64(s111meta.attrs.get('westBoundLongitude'))
	maxLon = np.float64(s111meta.attrs.get('eastBoundLongitude'))	
	
	delta_lon = s111meta.attrs.get('gridSpacingLongitudinal') #step size
	delta_lat = s111meta.attrs.get('gridSpacingLatitudinal') #step size
	
    #maxAbsLat = max(abs(float(mask.Minimum_Latitude)),abs(float(mask.Maximum_Latitude)))
	maxAbsLat = max(abs(minLat),abs(maxLat))
	p1 = geo.Point(0,maxAbsLat)
	p2 = geo.Point(delta_lon,maxAbsLat-delta_lat)
	
	metadata['density_meters'] = geo.sphericalEarth.distanceCourseFromDeg(p1,p2)['distance']
	
	s = {'src':getOutfn(bname)}	
	s['type'] = 'grid'	
	s['georeference'] = {}		
	s['georeference']['bounds']={}	
	s['georeference']['bounds']['min'] ={'x':minLon, 'y':minLat}	
	s['georeference']['bounds']['max'] ={'x':maxLon, 'y':maxLat}
	s['georeference']['projection']='WGS84'	
	s['georeference']['density']={'x':np.float64(delta_lon),'y':np.float64(delta_lat)}	
	s['georeference']['origin'] = { "x":minLon, "y":maxLat}
	s['georeference']['stepSize'] = { "x":np.float64(delta_lon), "y":-np.float64(delta_lat)}	
	s['mask'] = 'mask'
	s['definition']= {}
	s['definition']['u'] = {"highBit":31, "lowBit":20, "base":-20.48, "scale":0.01}	
	s['definition']['v'] = {"highBit":19, "lowBit":8, "base":-20.48, "scale":0.01}
	s['definition']['mask'] = {"highBit":7, "lowBit":0, "base":0, "scale":1}	
	metadata['sources'].append(s)
		
	return metadata 
#---------------------------------
def populateDataContainer(dataIn, dataType, basename):
	dataContainer = {}
	dataContainer['baseFileName'] = basename
	dataContainer['row'] = getDimCount(dataIn, dataType,'row')
	dataContainer['col'] = getDimCount(dataIn, dataType, 'col')
	dataContainer['maskVal']  = getMaskVal(dataIn, dataType)
	dataContainer['speedVal'] = getSpeedVal(dataIn, dataType)
	dataContainer['dirVal'] = getSpeedVal(dataIn, dataType)
	dataContainer['bounds'] = getBoundingBox(dataIn, dataType)
	return dataContainer
	
def getDimCount(dataIn,dataType,dim):
	if dataType == 'slgo':				
		if dim == 'row':
			return dataIn.variables['mask'].shape[0]
		else:
			return dataIn.variables['mask'].shape[1]
		
def getOutjfn(bname):
	return bname +'.json'
	
def getOutfn(bname,id = ''):
	return bname+id+'.png' 

def createFiles(dataContainer,dataType,basename):
	#dataContainer = populateDataContainer(dataContainer,dataType,basename)
	createPNG(dataContainer)
	metadata = createJSON(dataContainer)
	return metadata
#---------------------------------
  
#what are the expected arguments for this program?
#asumming: %>python s1112png.py S111US_20180529T12Z_CBOFS_TYP2.h5
jsonfn = None

#prepare the datasets array, this will be 
datasets = []

#get the command line argument after "python" 
#is this expecting multiple files of any type at anytime?
for arg in sys.argv[1:]:
	print(arg)
	#want the extension - TODO: we only accept .hd5 so this needs to change
	basename = os.path.splitext(arg)[0]
	ext = os.path.splitext(arg)[1]
	if ext == '.h5':
		hd5fn = arg
	else:
		print('this program only accepts .hd5 file types')
		break;
  
	jsonfn = getOutjfn(basename)
	infile = sys.argv[1] #get the s111 file from the command line
	print(infile)
	datasets = [] #initialize the dataset container
	h5 = h5py.File(infile,'r') #read in the surface current dataset

	surfcurGroup = h5['SurfaceCurrent'] # the group containing our data
	groupF = h5['Group_F']
	gfsc = groupF['SurfaceCurrent']
	sfill = float(gfsc[0]['fillValue'])
	
	for key in surfcurGroup:
		if key != 'axisNames': #all the other groups will be SurfaceCurrent## for the specific location it is at
			groups = surfcurGroup[key] 
	 
			#this loops through the time series data for each location
			for groupName in groups: 
				group = groups[groupName] #one time series slice 
				s111 = group['values']			
				#append each slice to the dataset 
				#datasets.append(convertS111(s111,basename,surfcurGroup)) #is this how it was meant to work with time series?
				#createPNG(s111,bname)
				metadata = convertS111(s111,basename,groups,sfill)
				datasets.append(metadata)
				break;
			break;
	print('dataset populated')
	
	if jsonfn is not None:
		jsonFile = open(jsonfn,'w')
		json.dump(datasets,jsonFile,indent=2) #this is how to do it without storing in local memory first