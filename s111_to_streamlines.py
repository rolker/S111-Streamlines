#!/usr/bin/env python

import sys
import os
import math
import json
import datetime
import dateutil.parser
import shutil
import h5py

earth_radius = 6371000.0

class Streamline:
    '''Contains the points and metadata for one streamline.
    
    Args:
        seed (Point): Starting location of the streamline. The streamline may grow in either direction,
        so this point is not necessarily at an end of the streamline.
        level (int): Level at which this streamline begins to be shown. 
    '''
    def __init__(self,seed,level):
        self.points = [seed]
        seed.level = level
        self.level = level
        self.seedIndex = 0
        self.index = None
        self.bounds = Bounds()


    def addPoint(self,p,direction):        
        if direction > 0:
            self.points.append(p)
        else:
            self.points.insert(0,p)
            self.seedIndex += 1

        self.bounds.add(p)

    
    def asDict(self):
        ret = {'level':self.level,'index':self.index,'seedIndex':self.seedIndex,'points':[],'bounds':self.bounds.asDict()}
        for p in self.points:
            dp = p.degrees()
            ret['points'].append({'x':dp.x,'y':dp.y,'level':p.level,'magnitude':p.flow.magnitude})
        return ret

class JobardLefer:
    '''Generates streamlines using the Jobard Lefer algorithm.
    
    Attributes:
        seperationFactor (float): Used in determining the base grid size. Factor is applied to field density when
            determining dSep parameter.
        testFactor (float): Factor used to calculate dTest relative to dSep.
        iSteps (int): Number of intermediate steps to use when calculating next point along the streamline.
        dSepMaxFactor (float): Maximum dSep relative to data extent when using overview layers. 
        minMag (float): Minimum magnitude to use when extending a streamline. If magnitude fall below this
            value, it is considered null.
    '''
    def __init__(self):
        self.separationFactor = 1.5
        self.testFactor = 0.5
        self.iSteps = 5
        self.dSepMaxFactor = 3.75
        self.minMag = 0.0001
        self.levelFactor = 1

    def generate(self,field):
        '''Generates the streamlines for the given field.
            Once a seed point has been selected in a field, they make a streamline growing beyond that point backward and forward.
            The growing process is stopped when the streamline reaches an edge of the surface, a singularity in the field (source or sink)
            or becomes too close to another streamline. (for rendering: The streamline is then divided into a set of small segments of contrasting color and
            projected onto the surface)

            Computing a new streamline is achieved in the following manner. A new seed point is chosen at a minimal distance apart from all existing streamlines.
            Then a new streamline is integrated beyond the seed point backward and forward until either it goes too close from another streamlines or it leaves
            the 2D domain over which the computation takes place. The algorithm ends when no more valid seed points can be found.
        
        Args:
            field (Field): The field for which streamlines will be generated.
            
        Returns:
            Dictionary with streamlines and some parameters.
            
            - streamlines (array): The generated streamlines.
            - dSep: The dSep parameter used to generate the streamlines. It is the separating distance given by the user.
                    It represents the minimal distance between seed points and streamlines.
            - iSteps: The number of intermediate points used while generating the streamlines.

        Notes: dTest: is a percentage of dSep. It corresponds to the minimal distance under which integration of the streamlines will be stopped in the current direction.
                       typically, dTest = dSep * .05 gives good visual results by producing long streamlines.
        '''
        self.field = field
        self.bounds = field.bounds
        print('field: ',self.field)
        
        print('bounds: ',self.bounds)

        density = field.getDensity()
        print('density:',density)
        self.dSep = density * self.separationFactor 
        self.dTest = self.dSep*self.testFactor
        print('dSep:',self.dSep,'dTest:',self.dTest)

        # calculate the grid size for the points grid based on dSep where
        # lat is closest to the equator and diff in long is largest
        if(self.bounds.min.y < 0 and self.bounds.max.y > 0):
            minLat = 0.0
        else:
            minLat = min(abs(self.bounds.min.y),abs(self.bounds.max.y))

        print('minLat: ',minLat)
        p0 = Point(0,minLat)
        pdx = positionFromDistanceCourse(p0, self.dSep, 1.5708)
        pdy = positionFromDistanceCourse(p0, self.dSep, 0.0)
        dx = pdx.x - p0.x
        dy = pdy.y - p0.y
        
        self.pointsGridCellSpacing = Point(dx, dy)
        print('points grid cell spacing:',self.pointsGridCellSpacing)
        
        size = self.bounds.getSize()
        print('size:',size)
        
        self.pointsGrid = {}
        
        # As the latitude varies in the points grid, its coverage relative to dSep varies,
        # this factor is kept track of in the following dict.
        self.pointsGridCellWidthFactors = {};
            
        self.minLevel = 0;
        dSepMax = min(size.x/self.pointsGridCellSpacing.x,size.y/self.pointsGridCellSpacing.y)/self.dSepMaxFactor
        self.minLevel = int(-math.floor(math.log(dSepMax,2)))
        print('min level:',self.minLevel)
        
        # In Jobard-Lefer paper, all the seeds are generated from the initial streamline.
        #http://web.cs.ucdavis.edu/~ma/SIGGRAPH02/course23/notes/papers/Jobard.pdf
        seedCache = []
        seedSpacing = Point(max(self.pointsGridCellSpacing.x*2,size.x/250),max(self.pointsGridCellSpacing.y*2,size.y/250))
        print('seedSpacing:',seedSpacing)
        center = self.bounds.getCenter()
        x = seedSpacing.x/2.0
        #print('x:',x)
        while x < size.x/2.0:
            y = seedSpacing.y/2.0
            #print('y:',y)
            while y < size.y/2.0:
                for i in range(-1,2,2):
                    for j in range(-1,2,2):
                        seed = Point(center.x+x*i,center.y+y*j)
                        #print('seed:', seed)
                        #FIX: none of the seeds are valid - how to make sure seeds are valid? Are the seeds indexes into the data flow field? The seed values rarely change
                        #NOTE: why is this one field.pointHasValue and not self.field.pointHasValue!?
                        
                        if self.field.pointHasValue(seed):                            
                            seedCache.append(seed)
                        else:
                            pass
                y += seedSpacing.y
            x += seedSpacing.x
        print('seed cache: ',seedCache)
        self.streamlines = []
        #QUESTION: why is the range only until 1? How are these values related to the incoming data?
        for level in range(self.minLevel,1):
            self.level = level
            self.levelFactor = int(2**(-level))
            print(level,self.levelFactor)
            print('  seed cache size:',len(seedCache))
            # dSepEffective is dSep scaled for the current zoom level
            dSepEffective = self.dSep*self.levelFactor
            
            for sl in self.streamlines:
                self.extend(sl)
                # following seems to be a no-op now, maybe at some point, addStreamline did something to
                # points in a streamline that already existed, such add the points to the points grid? 
                #if self.extend(sl):
                #    self.addStreamline(sl)
                    
            slStart = 0
            keptSeeds = []
            for seed in seedCache:
                sli = slStart
                while sli < len(self.streamlines):
                    sl = self.streamlines[sli]
                    sli += 1
                    # this looks for candidate seed points as described in JL paper fig 3
                    for pn in range(0,len(sl.points),self.iSteps*self.levelFactor):                      
                      dir = sl.points[pn].flow.direction
                      # check either side (perpendicularly) of the streamline for new seed points.
                      for k in range(2):
                          newSeed = positionFromDistanceCourse(sl.points[pn],dSepEffective,dir+1.57079633+k*math.pi)
                          if self.isPointGood(newSeed, dSepEffective):
                              newSl = Streamline(newSeed,self.level)
                              self.extend(newSl)
                              if len(newSl.points) >  2:
                                  self.addStreamline(newSl)
                slStart = len(self.streamlines)
                if self.isPointGood(seed, dSepEffective):
                   # print(seed,seed.flow.direction,seed.flow.u,seed.flow.v)
                    newSl = Streamline(seed,self.level)
                    self.extend(newSl)
                    if len(newSl.points)>2:
                        self.addStreamline(newSl)
                    else:
                        keptSeeds.append(seed)
                else:
                    if self.isPointGood(seed, self.dSep):
                        keptSeeds.append(seed)
            seedCache = keptSeeds
        return {'dSep':self.dSep, 'streamlines':self.streamlines, 'iSteps':self.iSteps}

    def getPointsGridIndex(self,p):   
        '''Calculates the index into the points grid. The points grid is used to make sure that
        a new points isn't too close to an existing point.
        '''
        xi = int(math.floor((p.x-self.bounds.min.x)/self.pointsGridCellSpacing.x))
        yi = int(math.floor((p.y-self.bounds.min.y)/self.pointsGridCellSpacing.y))
        return xi,yi

        
    def addPoint(self,p,streamIndex):
        '''Adds a point to the points grid, which is used to keep track of points not available
        to be added to streamlines. The points grid is implemented as a dictionary of dictionaries
        to allow it to be sparse, in other words, cells are created on demand. Each cell contains a
        list of points.
        
        Args:
            p (Point): Point to be added.
            streamIndex (int): Index of the streamline to which point p belongs.
        '''
        
        xi,yi = self.getPointsGridIndex(p)
        
        if not yi in self.pointsGrid:
            self.pointsGrid[yi] = {}
            p0 = Point(0,yi*self.pointsGridCellSpacing.y+self.bounds.min.y)
            p1 = positionFromDistanceCourse(p0, self.dSep, 1.5708) # 90 degrees in radians
            #calculates a factor to know how much dsep spans in terms of grid spacing
            self.pointsGridCellWidthFactors[yi] = p1.x/self.pointsGridCellSpacing.x
        if not xi in self.pointsGrid[yi]:
            self.pointsGrid[yi][xi] = []
        self.pointsGrid[yi][xi].append({'point':p,'streamIndex':streamIndex})

    def isPointGood(self, p, sep, streamIndex=None):
        '''Checks if point is acceptable.
        
        Args:
            p (Point): Point to be checked.
            sep (float): Distance in meters point should be from other points.
            streamIndex (int):  
        '''
        
        if not self.field.pointHasValue(p):
            return False
        
        xi,yi = self.getPointsGridIndex(p)

        for i in range(yi-self.levelFactor, yi+self.levelFactor+1):
            if i in self.pointsGrid:
                lFactor = int(math.ceil(self.levelFactor*self.pointsGridCellWidthFactors[i]))
                for j in range(xi-lFactor,xi+lFactor+1):
                    if j in self.pointsGrid[i]:
                        for gp in self.pointsGrid[i][j]:
                            if gp['streamIndex'] != streamIndex and distanceCourseFromRad(gp['point'], p)['distance']<sep:
                                return False
                                
        return True

    def isStreamPointGood(self,sl,p,extraPoints):
        if len(extraPoints) == self.levelFactor*self.iSteps:
            for i in range(0,len(sl.points),self.iSteps*self.levelFactor):
                if distanceCourseFromRad(sl.points[i], p)['distance']<self.dTest:
                    return False
        for i in range(0,len(sl.points),self.iSteps):
            if distanceCourseFromRad(sl.points[i], p)['distance']<self.dTest:
                return False
        for i in range(0,len(extraPoints)-self.iSteps,self.iSteps):
            if distanceCourseFromRad(extraPoints[i], p)['distance']<self.dTest:
                return False
        return True

    def addStreamline(self,sl):
        '''Adds a streamline to the collection of streamlines.
        
        Args:
            sl (jl.Streamline): Streamline to be added.
        '''
        if sl.index is None:
            sl.index = len(self.streamlines)
            self.streamlines.append(sl)
            for i in range(0,len(sl.points),self.iSteps):
                self.addPoint(sl.points[i], sl.index)

    def extend(self,stream):
        '''Extends a streamline in both directions from a seed point.
        
        Args:
            stream (jl.Streamline): Streamline initialized with a seed point.
        '''
        ok = True
        extended = False
        while ok:
            ok = self.step(stream,1)
            if ok:
                extended = True
        ok = True
        while ok:
            ok = self.step(stream,-1)
            if ok:
                extended = True
        return extended

    def step(self,sl,direction):
        '''Attempts to add a point to the streamline.
        
        Args:
            sl (jl.Streamline): Streamline to be extended.
            direction (int): 1 if in the flow direction, -1 if in the direction oposite the flow.
            
        Returns:
            True if step was succefull, False otherwise.
        '''
        ok = True
        if direction > 0:
            p0 = sl.points[-1]
        else:
            p0 = sl.points[0]
        
        transportOptions = {
            'distance':self.dSep*direction,
            'steps':self.iSteps,
            'minimumMagnitude':self.minMag
        }
        
        if self.field.pointHasValue(p0):
            steps = []
            pLast = p0
            while pLast is not None and len(steps) < self.iSteps*self.levelFactor:
                iSteps = self.field.transport(pLast,transportOptions)
                if len(iSteps) == self.iSteps:
                    for ip in iSteps:
                        if self.isPointGood(ip, self.dTest*self.levelFactor, sl.index):
                            pLast = ip
                        else:
                            pLast = None
                            break
                    if pLast is not None:
                        steps += iSteps
                else:
                    pLast = None
                if pLast is not None:
                    if not self.isStreamPointGood(sl, pLast, steps):
                        pLast = None
            if pLast is not None:
                for ip in steps:
                    ip.level = self.level
                    if sl.index is not None:
                        self.addPoint(ip, sl.index)
                    sl.addPoint(ip,direction)
            else:
                ok = False
        else:
            ok = False
        return ok


def interpolate(v1,v2,p):
    if v1 is None and v2 is None:
        return None
    
    if v1 is None:
        u2 = math.sin(v2.direction)*v2.magnitude*p
        v2 = math.cos(v2.direction)*v2.magnitude*p
        mag = math.sqrt(u2*u2+v2*v2)
        dir = math.atan2(u2, v2)

        return Flow(mag,dir)
    
    ip = 1.0-p
    
    if v2 is None:
        u1 = math.sin(v1.direction)*v1.magnitude*ip
        v1 = math.cos(v1.direction)*v1.magnitude*ip
        mag = math.sqrt(u1*u1+v1*v1)
        dir = math.atan2(u1, v1)

        return Flow(mag,dir)

    u1 = math.sin(v1.direction)*v1.magnitude*ip
    v1 = math.cos(v1.direction)*v1.magnitude*ip
    u2 = math.sin(v2.direction)*v2.magnitude*p
    v2 = math.cos(v2.direction)*v2.magnitude*p
    
    u = u1+u2
    v = v1+v2

    mag = math.sqrt(u*u+v*v)
    dir = math.atan2(u, v)

    return Flow(mag,dir)

class Flow():
    def __init__(self,mag,direction):
        self.magnitude = mag
        self.direction = direction

class FlowField():
    '''A field representing currents.
    '''
    def __init__(self,data,timestep, metadata):
        self.data = data
        self.metadata = metadata #TODO: eventually, copy over just the metadata portion...this is sending the entire group of the timeseries!
        self.timestep = timestep
    
        self.dx = math.radians(self.metadata.attrs['gridSpacingLongitudinal'])
        self.dy = math.radians(self.metadata.attrs['gridSpacingLatitudinal'])
        
        self.minPoint = Point(self.metadata.attrs['westBoundLongitude'],self.metadata.attrs['southBoundLatitude']).radians()
        self.maxPoint = Point(self.metadata.attrs['eastBoundLongitude'],self.metadata.attrs['northBoundLatitude']).radians()
        print("init flowField: hdf5 bounds")
        print(self.metadata.attrs['westBoundLongitude'],self.metadata.attrs['southBoundLatitude'])
        print(self.metadata.attrs['eastBoundLongitude'],self.metadata.attrs['northBoundLatitude'])
        
        self.bounds = Bounds()
        self.bounds.add(self.minPoint)
        self.bounds.add(self.maxPoint)

    def transport(self,startPoint,options):
        '''Calculates the next point in the direction of flow.

        Args:
            startPoint (Point): Initial geographic position.
            options (dict): Options used in calculating next point.

                - distance: distance in meters away from startPoint that the resulting point should be.
                - steps: number of intermediate steps to use in calculating next point.
                - minimumMagnitude: Minimum magnitude for which flow at all intermediate steps should be
                    for flow to continue.
                
        Returns:
            If transport is successful, list of intermediate points and final point resulting from the transport.
                Empty list if unable to complete the transport.
        '''
        ret = []
        lastPoint = startPoint
        stepSize = options['distance']/float(options['steps'])
        while len(ret) < options['steps']:
            if self.pointHasValue(lastPoint) and lastPoint.flow.magnitude > options['minimumMagnitude']:
                pMid = positionFromDistanceCourse(lastPoint, stepSize/2.0, lastPoint.flow.direction)
                if self.pointHasValue(pMid) and pMid.flow.magnitude > options['minimumMagnitude']:
                    pi = positionFromDistanceCourse(lastPoint, stepSize, pMid.flow.direction)
                    ret.append(pi)
                    lastPoint = pi
                else:
                    lastPoint = None
                    break
            else:
                break
        return ret

    def pointHasValue(self,p):
        '''Checks if a flow value exists for a given points, adding it to the point if it exists.
        
        Args:
            p (Point): Geographic position where flow values is desired.
            
        Returns:
            True if flow exists at given position, False otherwise. If flow exists, it is attached
            to input parameter p.
        '''
       
        if p is None:
            return None
        if not 'flow' in p.__dict__:
            p.flow = self.getFlow(p)
        return p.flow is not None

    def getFlow(self,p):
        '''Retrieves direction and magnitude at a given point.
        
        Args:
            p (Point): Position where data is requested.
        
        Returns:
            An instance of flow.Flow.
        '''
        
        if self.bounds.contains(p):
            index = self.getIndex(p)
            if index is None:
                return None
            
            x1 = math.floor(index[0])
            x2 = math.ceil(index[0])
            y1 = math.floor(index[1])
            y2 = math.ceil(index[1])

            px = index[0]-x1
            py = index[1]-y1
            
            f11 = self.getFlowAtIndex(x1,y1)
            f12 = self.getFlowAtIndex(x1,y2)
            f21 = self.getFlowAtIndex(x2,y1)
            f22 = self.getFlowAtIndex(x2,y2)
            
            ret = interpolate(interpolate(f11,f12,py),interpolate(f21,f22,py),px)
                
            return ret
            
            
            
        
    def getFlowAtIndex(self,x,y):
        #if x >= 0 and x < self.data.attrs['numberPointsLong'] and y >= 0 and y < self.data.attrs['numberPointsLat']:
        if x >= 0 and x < self.metadata.attrs['numPointsLongitudinal'] and y >= 0 and y < self.metadata.attrs['numPointsLatitudinal']:
            speed = self.data['values'][y,x,'SurfaceCurrentSpeed']
            dir = self.data['values'][y,x,'SurfaceCurrentDirection']
            #speed = self.data[self.timestep]['Speed'][y,x]
            #dir = self.data[self.timestep]['Direction'][y,x]            
            if speed >= 0.0:
                return Flow(speed,math.radians(dir))            

    def getIndex(self,p):
        return ((p.x-self.minPoint.x)/self.dx,(p.y-self.minPoint.y)/self.dy)

    
    def getDensity(self):
        '''Calculates the grid size or data density of the field.
        
        Since the grid is in latitudes and longitudes, the grid size is not consistent in distance based
        units such as meters. This methods attempts to find the section of the grid with the smallest
        density and calculates it in meters.
        
        Returns:
            Smallest grid size in meters.
        '''
        
        maxLat = max(abs(self.bounds.min.y),abs(self.bounds.max.y))
        
        return distanceCourseFromRad(Point(0.0,maxLat-self.dy),Point(self.dx,maxLat))['distance']


def distanceCourseFromRad(p1,p2):
    '''Calculates distance and course between two points.
    
    Args:
        p1 (Point): Starting point expressed in radians.
        p2 (Point): Ending point expressed in radians.
        
    Returns:
        Dictionary with distance and course.
        
        - distance: Distance in meters.
        - course_radians: Course in radians.
    '''
    dlon = p2.x-p1.x

    clat1 = math.cos(p1.y);
    clat2 = math.cos(p2.y);
    slat1 = math.sin(p1.y);
    slat2 = math.sin(p2.y);
    tlat2 = math.tan(p2.y);

    cdlon = math.cos(dlon);
    sdlon = math.sin(dlon);

    y = math.sqrt(math.pow(clat2*sdlon,2)+math.pow(clat1*slat2-slat1*clat2*cdlon,2))
    x = slat1*slat2+clat1*clat2*cdlon
    central_angle = math.atan2(y,x)

    course = math.atan2(sdlon,clat1*tlat2-slat1*cdlon)

    return {'distance':central_angle*earth_radius, 'course_radians':course}

def positionFromDistanceCourse(p1_rad,distance,course_rad):
    '''Finds position at given distance and direction from an initial point.
    
    Args:
        p1_rad (Point): Initial position expressed in radians.
        distance (float): Distance in meters.
        course_rad (float): Direction in radians.
        
    Returns:
        Position as a Point expressed in radians.
    '''
    slat1 = math.sin(p1_rad.y)
    clat1 = math.cos(p1_rad.y)
        
    central_angle = distance/earth_radius
        
    cca = math.cos(central_angle)
    sca = math.sin(central_angle)
        
    ccourse = math.cos(course_rad)
    scourse = math.sin(course_rad)
        
    y = slat1*cca+clat1*sca*ccourse
    x = math.sqrt(math.pow(clat1*cca-slat1*sca*ccourse,2)+math.pow(sca*scourse,2))
    lat2 = math.atan2(y,x)
        
    y = sca*scourse
    x = clat1*cca-slat1*sca*ccourse
        
    dlon = math.atan2(y,x)
    lon2 = p1_rad.x+dlon
        
    return Point (lon2, lat2);


class Point:
    '''Represents a geographic point.
    
    Note:
        Point may be represented in degrees or radians, but the class doesn't
        keep track of which unit is used.
    '''
    def __init__(self,x,y):
        self.x = x
        self.y = y

    def radians(self):
        '''Returns a copy of the point expressed in radians. Point must be in degrees.
        '''
        return Point(math.radians(self.x),math.radians(self.y))

    def degrees(self):
        '''Returns a copy of the point expressed in degrees. Point must be in radians.
        '''
        return Point(math.degrees(self.x),math.degrees(self.y))

    def __str__(self):
        return '{:.12f},{:.12f}'.format(self.x,self.y)
    
    def __repr__(self):
      return self.__str__()

class Bounds:
    '''2D bounding box.
    
    Keeps track of a rectangular box containing all specified points.
    May be initialized empty.
    
    Args:
        p1 (Point): Optional first point of bounding box.
        p2 (Point): Optional second point of bounding box.
    '''
    def __init__(self,p1=None,p2=None):
        if p2 is None:
            p2 = p1
        if p1 is None:
            self.min = None
            self.max = None
        else:      
            self.min = Point(min(p1.x,p2.x),min(p1.y,p2.y))
            self.max = Point(max(p1.x,p2.x),max(p1.y,p2.y))
      
    def empty(self):
        '''Returns True if empty.'''
        return self.min is None
      
    def add(self,p):
        '''Adds a point to the bounding box, expanding or initializing it if necessary.
        
        Args:
            p (Point): Point to be added.
        '''
        if self.min is None:
            self.min = p
            self.max = p
        else:
            self.min = Point(min(self.min.x,p.x),min(self.min.y,p.y))
            self.max = Point(max(self.max.x,p.x),max(self.max.y,p.y))

    def getSize(self):
        '''Returns the size of the bounding box.
        
        Note:
            Must not be empty!
        '''
        return Point(self.max.x-self.min.x,self.max.y-self.min.y)

    def __str__(self):
        return 'min:{} max:{}'.format(self.min,self.max)

    def getCenter(self):
        '''Returns the center point of the bounding box.
        
        Note:
            Must not be empty!
        '''
        size = self.getSize()
        return Point(self.min.x+size.x/2.0,self.min.y+size.y/2.0)
      
    def contains(self,p):
      '''Returns True if p is inside bounds or False otherwise or if bounds is empty.'''
      if self.min is None:
          return False
      return p.x >= self.min.x and p.x <= self.max.x and p.y >= self.min.y and p.y <= self.max.y
  
    def asDict(self):
        dmin = self.min.degrees()
        dmax = self.max.degrees()
        return {'min':{'x':dmin.x, 'y':dmin.y}, 'max':{'x':dmax.x, 'y':dmax.y}}

infile = sys.argv[1]
jsonfn = infile+'.json'
jfile=open(jsonfn,"w")

datasets = []

jlContext = JobardLefer()
    
dataset = h5py.File(infile,'r')

surfcurGroup = dataset['SurfaceCurrent']
#print(type(surfcurGroup))
#NOTE: the group "SurfaceCurrent" contains sets of surface current data, each in a different geographical location
for key in surfcurGroup:
    if key != 'axisNames':
        groups = surfcurGroup[key]
        #print(type(groups))
        i = 0
        #meta
        #this loops through the time series data for each location
        for groupName in groups:
            #print(groupName)
            group = groups[groupName]
            #each surface current group will need to have it's own flow field, more than one is a timeseries            
            #'DateTime' no longer in s111 version as of 5/22/18
            val = dateutil.parser.parse(group.attrs['timePoint'])
            #timestamp is really the name of the group in FlowField which in current s111 version is just Group_###
            #print(val)            

            dataModel = FlowField(group,groupName,groups)
            i=i+1
            #print(i)
            
##
##    # generate returns dictionary with streamlines and some parameters (dSep and iSteps, which are related to spacing of the streamlines)
            streamlines = jlContext.generate(dataModel)
            print('streamlines created')
            print(streamlines)
            
            
##        
##    # Turn the streamlines into a python dictionary representation so the output can be tuned into a json file.
##    # The json file can be used with a browser based representation overlaid on Google Maps.
            sldicts = []
            bounds = Bounds()
            for sl in streamlines['streamlines']:
                b = sl.bounds
                bounds.add(b.min)
                bounds.add(b.max)
                sldicts.append(sl.asDict())
                streamlines['streamlines'] = sldicts
                streamlines['time'] = val.isoformat()
                streamlines['label'] = val.isoformat()
                streamlines['bounds'] = bounds.asDict()
                datasets.append(streamlines)
            break
##
jsonout = json.dumps(datasets,indent=2)
#print(jsonout)

if jfile is not None:    
  #file(jsonfn,'w').write(jsonout)
    jfile.write(jsonout)
  #shutil.copy(jsonfn, jsonfn_latest)
  
    
