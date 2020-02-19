#ifndef S111_STREAMLINES_H
#define S111_STREAMLINES_H

#include <H5Cpp.h>
#include <vector>
#include <iostream>
#include <map>
#include <set>
#include <deque>
#include <cmath>
#include <algorithm>

#define S111_EARTH_RADIUS 6371000.0

namespace s111
{

    /**
     * Latitude and Longitude position in radians
     */
    struct LatLong
    {
        LatLong()
        {
        }
        LatLong(double lat, double lon):latitude(lat),longitude(lon)
        {
        }
        
        static LatLong fromDegrees(double lat, double lon)
        {
            return LatLong(lat*M_PI/180.0,lon*M_PI/180.0);
        }
        
        double latitude_degrees() const
        {
            return 180.0*latitude/M_PI;
        }

        double longitude_degrees() const
        {
            return 180.0*longitude/M_PI;
        }
        
        double latitude;
        double longitude;
    };
    
    std::ostream& operator<<(std::ostream& os, const LatLong& ll)
    {
        os << "(" << ll.latitude_degrees() << ", " << ll.longitude_degrees() << ")";
        return os;
    }

    struct SpeedDirection
    {
        float speed; 
        float direction;
    };
    
    struct Flow
    {
        Flow():valid(false){}
        Flow(SpeedDirection const &sd):magnitude(sd.speed),direction(sd.direction*M_PI/180.0),valid(true)
        {}
        Flow(double mag, double dir):magnitude(mag),direction(dir),valid(true)
        {}
        
        float magnitude;
        
        /// direction in radians
        float direction;
            
        LatLong position;
        bool valid;
    };

    Flow interpolate(Flow const &f1, Flow const &f2, float p)
    {    
        if(!f1.valid && !f2.valid)
            return Flow();
    
        if (!f1.valid)
        {
            double u2 = sin(f2.direction)*f2.magnitude*p;
            double v2 = cos(f2.direction)*f2.magnitude*p;
            double mag = sqrt(u2*u2+v2*v2);
            double dir = atan2(u2, v2);

            return Flow(mag,dir);
        }
    
        float ip = 1.0-p;
    
        if (!f2.valid)
        {
            double u1 = sin(f1.direction)*f1.magnitude*ip;
            double v1 = cos(f1.direction)*f1.magnitude*ip;
            double mag = sqrt(u1*u1+v1*v1);
            double dir = atan2(u1, v1);

            return Flow(mag,dir);
        }

        double u1 = sin(f1.direction)*f1.magnitude*ip;
        double v1 = cos(f1.direction)*f1.magnitude*ip;
        double u2 = sin(f2.direction)*f2.magnitude*p;
        double v2 = cos(f2.direction)*f2.magnitude*p;
    
        double u = u1+u2;
        double v = v1+v2;

        double mag = sqrt(u*u+v*v);
        double dir = atan2(u, v);

        return Flow(mag,dir);
    }
    
    
    
    struct DistanceCourse
    {
        /// distance in meters
        double distance;
        
        /// course in radians
        double course;
    };
    
    /**
    Calculates distance and course between two points.
    
    Args:
        p1 (Point): Starting point expressed in radians.
        p2 (Point): Ending point expressed in radians.
        
    Returns:
        Dictionary with distance and course.
        
        - distance: Distance in meters.
        - course_radians: Course in radians.
    */
    DistanceCourse distanceCourseFromRad(LatLong const &p1, LatLong const &p2)
    {
        double dlon = p2.longitude-p1.longitude;

        double clat1 = cos(p1.latitude);
        double clat2 = cos(p2.latitude);
        double slat1 = sin(p1.latitude);
        double slat2 = sin(p2.latitude);
        double tlat2 = tan(p2.latitude);

        double cdlon = cos(dlon);
        double sdlon = sin(dlon);

        double y = sqrt(pow(clat2*sdlon,2)+pow(clat1*slat2-slat1*clat2*cdlon,2));
        double x = slat1*slat2+clat1*clat2*cdlon;
        double central_angle = atan2(y,x);

        DistanceCourse ret;
        ret.course = atan2(sdlon,clat1*tlat2-slat1*cdlon);
        ret.distance = central_angle*S111_EARTH_RADIUS;
        return ret;
    }
    

    /**
    Finds position at given distance and direction from an initial point.
    
    Args:
        p1_rad (Point): Initial position expressed in radians.
        distance (float): Distance in meters.
        course_rad (float): Direction in radians.
        
    Returns:
        Position as a Point expressed in radians.
    */
    LatLong positionFromDistanceCourse(LatLong const &p1, double distance, double course_rad)
    {
        double slat1 = sin(p1.latitude);
        double clat1 = cos(p1.latitude);
        
        double central_angle = distance/S111_EARTH_RADIUS;
        
        double cca = cos(central_angle);
        double sca = sin(central_angle);
        
        double ccourse = cos(course_rad);
        double scourse = sin(course_rad);
        
        double y = slat1*cca+clat1*sca*ccourse;
        double x = sqrt(pow(clat1*cca-slat1*sca*ccourse,2)+pow(sca*scourse,2));
        double lat2 = atan2(y,x);
        
        y = sca*scourse;
        x = clat1*cca-slat1*sca*ccourse;
        
        double dlon = atan2(y,x);
        double lon2 = p1.longitude+dlon;
        
        return LatLong(lat2,lon2);
    }

    /** 2D bounding box.
    
    Keeps track of a rectangular box containing all specified points.
    **/
    class Bounds
    {
    public:
        Bounds():m_empty(true)
        {}
        Bounds(LatLong const &p):m_empty(false),m_min(p),m_max(p)
        {}
        Bounds(LatLong const &p1, LatLong const &p2):m_empty(false)
        {
            m_min.latitude = std::min(p1.latitude,p2.latitude);
            m_min.longitude = std::min(p1.longitude,p2.longitude);
            m_max.latitude = std::max(p1.latitude,p2.latitude);
            m_max.longitude = std::max(p1.longitude,p2.longitude);
        }
        
        bool empty() const
        {
            return m_empty;
        }
        
        void add(LatLong const &p)
        {
            if(m_empty)
            {
                m_min = p;
                m_max = p;
                m_empty = false;
            }
            else
            {
                m_min.latitude = std::min(p.latitude,m_min.latitude);
                m_min.longitude = std::min(p.longitude,m_min.longitude);
                m_max.latitude = std::max(p.latitude,m_max.latitude);
                m_max.longitude = std::max(p.longitude,m_max.longitude);
            }
        }
        
        void add(Bounds const &b)
        {
            if(!b.empty())
            {
                add(b.m_min);
                add(b.m_max);
            }
        }
        
        LatLong const &min() const
        {
            return m_min;
        }

        LatLong const &max() const
        {
            return m_max;
        }
        
        LatLong getSize() const
        {
            return LatLong(m_max.latitude-m_min.latitude,m_max.longitude-m_min.longitude);
        }
        
        LatLong getCenter() const
        {
            LatLong size = getSize();
            return LatLong(m_min.latitude+size.latitude/2.0,m_min.longitude+size.longitude/2.0);
        }
        
        bool contains(LatLong const &p) const
        {
            if(m_empty)
                return false;
            return p.latitude >= m_min.latitude && p.latitude <= m_max.latitude && p.longitude >= m_min.longitude && p.longitude <= m_max.longitude;
        }
        
    private:
        bool m_empty;
        LatLong m_min;
        LatLong m_max;
    };

    std::ostream& operator<<(std::ostream& os, const Bounds& b)
    {
        os << "[" << b.min() << " - " << b.max() << "]";
        return os;
    }
    
    struct StreamlinePoint: public Flow
    {
        StreamlinePoint():level(1)
        {}
        
        StreamlinePoint(Flow const &f):Flow(f),level(1)
        {}
        
        /// Level at which point is shown
        int level;
    };
    
    ///Contains the points and metadata for one streamline.
    struct Streamline
    {
        ///seed (Point): Starting location of the streamline. The streamline may grow in either direction,
        ///so this point is not necessarily at an end of the streamline.
        ///level (int): Level at which this streamline begins to be shown. 
        Streamline(StreamlinePoint seed, int streamline_level): level(streamline_level), seedIndex(0), index(-1)
        {
            seed.level = level;
            points.push_back(seed);
            bounds.add(seed.position);
        }
        
        void add(StreamlinePoint p, bool forward)
        {
            if(forward)
                points.push_back(p);
            else
            {
                points.push_front(p);
                seedIndex += 1;
            }
            
            bounds.add(p.position);
        }
        
        std::deque<StreamlinePoint> points;
        int level;
        int seedIndex;
        int index;
        Bounds bounds;
    };
    
    struct StreamlineSet
    {
        std::vector<Streamline> streamlines;
        Bounds bounds;
        float dSep;
        int iSteps;
    };
    
    struct Index
    {
        int x;
        double xr;
        int y;
        double yr;
    };
    
    class FlowField
    {
    public:
        
        typedef std::map<std::string,H5::DataSet> DataSetMap;
        
        struct SurfaceCurrentGroup
        {
            H5::H5File file;
            H5::Group group;
            LatLong origin;
            LatLong spacing;
            Bounds bounds;
            int longitudeCount, latitudeCount;
            
            DataSetMap timePoints;
            
            std::vector<SpeedDirection> dataCache;
            std::string cachedTimePoint;
            
            Index getIndex(LatLong const &p) const
            {
                double x,y;
                x = (p.longitude - origin.longitude)/spacing.longitude;
                y = (p.latitude - origin.latitude)/spacing.latitude;
                Index ret;
                ret.x = floor(x);
                ret.xr = x-ret.x;
                ret.y = floor(y);
                ret.yr = y-ret.y;
                return ret;
            }
        };
        
        FlowField(){}
        
        void addFile(std::string const &fname)
        {
            SurfaceCurrentGroup scg;
            H5::H5File f(fname, H5F_ACC_RDONLY);
            scg.file = f;

            auto scGroup = f.openGroup("/SurfaceCurrent/SurfaceCurrent.01");
            scg.group = scGroup;
            
            double lat, lon;
            scg.group.openAttribute("gridOriginLatitude").read(H5::PredType::INTEL_F64,&lat);
            scg.group.openAttribute("gridOriginLongitude").read(H5::PredType::INTEL_F64,&lon);
            scg.origin = LatLong::fromDegrees(lat,lon);
            std::cerr << fname << "\n\torigin: " << scg.origin << std::endl;
            
            scg.group.openAttribute("gridSpacingLatitudinal").read(H5::PredType::INTEL_F64,&lat);
            scg.group.openAttribute("gridSpacingLongitudinal").read(H5::PredType::INTEL_F64,&lon);
            scg.spacing = LatLong::fromDegrees(lat,lon);
            std::cerr << "\tspacing: " << scg.spacing << std::endl;

            scg.group.openAttribute("southBoundLatitude").read(H5::PredType::INTEL_F64,&lat);
            scg.group.openAttribute("westBoundLongitude").read(H5::PredType::INTEL_F64,&lon);
            scg.bounds.add(LatLong::fromDegrees(lat,lon));
            scg.group.openAttribute("northBoundLatitude").read(H5::PredType::INTEL_F64,&lat);
            scg.group.openAttribute("eastBoundLongitude").read(H5::PredType::INTEL_F64,&lon);
            scg.bounds.add(LatLong::fromDegrees(lat,lon));

            scg.group.openAttribute("numPointsLongitudinal").read(H5::PredType::STD_I32LE, &scg.longitudeCount);
            scg.group.openAttribute("numPointsLatitudinal").read(H5::PredType::STD_I32LE, &scg.latitudeCount);
            
            std::cerr << "\tbounds: " << scg.bounds << std::endl;
            
            //std::cerr << "obj count: " << scg.group.getNumObjs() << std::endl;
            for (int i = 0; i < scg.group.getNumObjs(); i++)
            {
                //std::cerr << "\t" << i << ": " << scg.group.childObjType(i) << std::endl;
                if (scg.group.childObjType(i) == H5O_type_t::H5O_TYPE_GROUP)
                {
                    auto g = scg.group.openGroup(scg.group.getObjnameByIdx(i));
                    //std::cerr << "\t\t" << g.getObjName() << std::endl;
                    if(g.attrExists("timePoint"))
                    {
                        H5std_string timePoint;
                        H5::Attribute tpAttribute = g.openAttribute("timePoint");
                        H5::StrType stype = tpAttribute.getStrType();
                        tpAttribute.read(stype, timePoint);
                        //std::cerr << "\t\t" << timePoint << std::endl;
                        auto d = g.openDataSet("values");
                        scg.timePoints[timePoint]=d;
                        m_timePoints.insert(timePoint);
                    }
                }
            }
            
            m_surfaceCurrentGroups.push_back(scg);
            m_bounds.add(scg.bounds);
        }
        
        std::set<std::string> const &timePoints() const
        {
            return m_timePoints;
        }
        
        Bounds const &bounds() const
        {
            return m_bounds;
        }
        
        /**
        Calculates the grid size or data density of the field.
        
        Since the grid is in latitudes and longitudes, the grid size is not consistent in distance based
        units such as meters. This methods attempts to find the section of the grid with the smallest
        density and calculates it in meters.
        
        Returns:
            Smallest grid size in meters.
        */
        double getDensity() const
        {
            double maxLat = std::max(abs(m_bounds.min().latitude),abs(m_bounds.max().latitude));
            auto spacing = m_surfaceCurrentGroups.front().spacing;
            return distanceCourseFromRad(LatLong(maxLat-spacing.latitude,0.0),LatLong(maxLat,spacing.longitude)).distance;
        }
        
        Flow getFlow(LatLong const &position, std::string const &timeStep)
        {
            if(m_bounds.contains(position))
            {
                Flow f11,f12,f21,f22;
                Index i11;
                LatLong p12,p21,p22;
                for(SurfaceCurrentGroup &scg: m_surfaceCurrentGroups)
                    if(scg.bounds.contains(position))
                    {
                        i11 = scg.getIndex(position);
                        f11 = getFlowAtIndex(scg,i11,timeStep);
                        // it's possible that we'll need to interpolate values from different tiles, so search for the tiles
                        // for every index set
                        p12.longitude = position.longitude;
                        p12.latitude = position.latitude + scg.spacing.latitude;
                        p21.longitude = position.longitude+scg.spacing.longitude;
                        p21.latitude = position.latitude;
                        p22.longitude = p21.longitude;
                        p22.latitude = p12.latitude;
                        break;
                    }
                for(SurfaceCurrentGroup &scg: m_surfaceCurrentGroups)
                {
                    if(scg.bounds.contains(p12))
                    {
                        Index i12 = scg.getIndex(p12);
                        f12 = getFlowAtIndex(scg,i12,timeStep);
                    }
                    if(scg.bounds.contains(p21))
                    {
                        Index i21 = scg.getIndex(p21);
                        f21 = getFlowAtIndex(scg,i21,timeStep);
                    }
                    if(scg.bounds.contains(p22))
                    {
                        Index i22 = scg.getIndex(p22);
                        f22 = getFlowAtIndex(scg,i22,timeStep);
                    }
                }
                Flow ret = interpolate(interpolate(f11,f12,1.0-i11.yr),interpolate(f21,f22,1.0-i11.yr),1.0-i11.xr);
                ret.position = position;
                return ret;
            }
            return Flow();
        }
        
        Flow getFlowAtIndex(SurfaceCurrentGroup &scg, Index const &index, std::string const &timeStep)
        {
            if(index.x >= 0 && index.x < scg.longitudeCount && index.y >= 0 && index.y < scg.latitudeCount)
            {
                if (scg.cachedTimePoint != timeStep)
                {
                    DataSetMap::const_iterator dataSet = scg.timePoints.find(timeStep);
                    if(dataSet != scg.timePoints.end())
                    {
                        auto space = dataSet->second.getSpace();
                        scg.dataCache.resize(scg.latitudeCount*scg.longitudeCount);
                        H5::CompType ctype(sizeof(SpeedDirection));
                        ctype.insertMember("surfaceCurrentSpeed",0,H5::PredType::IEEE_F32LE);
                        ctype.insertMember("surfaceCurrentDirection",sizeof(float),H5::PredType::IEEE_F32LE);
                        dataSet->second.read(scg.dataCache.data(),ctype);
                        scg.cachedTimePoint = timeStep;
                    }
                }
                // check again, in case caching didn't work.
                if (scg.cachedTimePoint == timeStep)
                {
                    if (scg.dataCache[index.y*scg.longitudeCount+index.x].speed > 0.0)
                        return Flow(scg.dataCache[index.y*scg.longitudeCount+index.x]);
                }

            }
            return Flow();
        }
        
        /**
         Calculates the next point in the direction of flow.

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
        */
        std::vector<Flow> transport(LatLong const &startPoint, double distance, int steps, double minimumMagnitude, std::string const &timeStep)
        {
            std::vector<Flow> ret;
            Flow lastPoint = getFlow(startPoint,timeStep);
            double stepSize = distance/double(steps);
            while(ret.size() < steps)
            {
                if(lastPoint.valid && lastPoint.magnitude >= minimumMagnitude)
                {
                    auto pMid = positionFromDistanceCourse(lastPoint.position, stepSize/2.0, lastPoint.direction);
                    auto fMid = getFlow(pMid, timeStep);
                    if(fMid.valid && fMid.magnitude >= minimumMagnitude)
                    {
                        auto pi = positionFromDistanceCourse(lastPoint.position, stepSize, fMid.direction);
                        auto fi = getFlow(pi, timeStep);
                        if(fi.valid && fi.magnitude >= minimumMagnitude)
                            ret.push_back(fi);
                        lastPoint = fi;
                    }
                    else
                    {
                        lastPoint = fMid;
                        break;
                    }
                }
                else
                    break;
            }
            return ret;
        }
        
    private:
        std::vector<SurfaceCurrentGroup> m_surfaceCurrentGroups;
        std::set<std::string> m_timePoints;
        Bounds m_bounds;
    };
    
    class FlowFieldTimeStep
    {
    public:
        FlowFieldTimeStep(FlowField &ff, std::string const &timeStep):m_flowField(ff),m_timeStep(timeStep) 
        {}
        
        double getDensity() const
        {
            return m_flowField.getDensity();
        }
        
        Bounds const &bounds() const
        {
            return m_flowField.bounds();
        }
        
        bool pointHasValue();

        Flow getFlow(LatLong const &position)
        {
            return m_flowField.getFlow(position, m_timeStep);
        }
        
        std::string const &timeStep() const
        {
            return m_timeStep;
        }
        
        std::vector<Flow> transport(LatLong const &position, double distance, int steps, double minimumMagnitude)
        {
            return m_flowField.transport(position, distance, steps, minimumMagnitude, m_timeStep);
        }
        
        
    private:
        FlowField &m_flowField;
        std::string m_timeStep;
    };
    
    /// Generates streamlines using the Jobard Lefer algorithm.
    class JobardLefer
    {
    public:
        JobardLefer(): m_seperationFactor(1.5), m_testFactor(0.5), m_iSteps(5), m_dSepMaxFactor(3.75), m_minMag(0.01)
        {
        }
        
        /**
         Generates the streamlines for the given field.
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
        **/
        StreamlineSet generate(FlowFieldTimeStep &field)
        {
            std::cerr << "time step: " << field.timeStep() << std::endl;
            StreamlineSet ret;
            GenerateContext context(field,ret);
            context.density = field.getDensity();
            std::cerr << "density: " << context.density << " meters" << std::endl;
            context.dSep = context.density * m_seperationFactor;
            context.dTest = context.dSep * m_testFactor;
            std::cerr << "dSep: " << context.dSep << ", dTest: " << context.dTest << std::endl;
            
            // calculate the grid size for the points grid based on dSep where
            // lat is closest to the equator and diff in long is largest
            double minLat;
            if(field.bounds().min().latitude < 0 && field.bounds().max().latitude > 0)
                minLat = 0.0;
            else
                minLat = std::min(fabs(field.bounds().min().latitude),fabs(field.bounds().max().latitude));

            std::cerr << "minLat: " << minLat*180/M_PI << std::endl;
            
            LatLong p0(minLat,0.0);
            LatLong pdx = positionFromDistanceCourse(p0, context.dSep, 1.5708);
            LatLong pdy = positionFromDistanceCourse(p0, context.dSep, 0.0);
            double dx = pdx.longitude - p0.longitude;
            double dy = pdy.latitude - p0.latitude;
        
            context.pointsGridCellSpacing = LatLong(dy,dx);
            std::cerr << "points grid cell spacing: " << context.pointsGridCellSpacing << std::endl;
            
            LatLong size = field.bounds().getSize();
            std::cerr << "size: " << size << std::endl;
            
            float dSepMax = std::min(size.longitude/context.pointsGridCellSpacing.longitude,size.latitude/context.pointsGridCellSpacing.latitude)/m_dSepMaxFactor;
            context.minLevel = int(-floor(log2(dSepMax)));
            std::cerr << "minLevel: " << context.minLevel << std::endl;
            
            std::vector<Flow> seedCache;
            LatLong seedSpacing(std::max(context.pointsGridCellSpacing.latitude*2.0,size.latitude/250.0),std::max(context.pointsGridCellSpacing.longitude*2.0,size.longitude/250.0));
            std::cerr << "seedSpacing: " << seedSpacing << std::endl;
            LatLong center = field.bounds().getCenter();
            double lon = seedSpacing.longitude/2.0;
            while(lon<size.longitude/2.0)
            {
                double lat = seedSpacing.latitude/2.0;
                while (lat < size.latitude/2.0)
                {
                    for(int i = -1; i < 2; i+=2)
                        for(int j = -1; j < 2; j+=2)
                        {
                            Flow f = field.getFlow(LatLong(center.latitude+lat*j,center.longitude+lon*i));
                            if(f.valid)
                                seedCache.push_back(f);
                        }
                    lat += seedSpacing.latitude;
                }
                lon += seedSpacing.longitude;
            }
            std::cerr << "seedCache size: " << seedCache.size() << std::endl;
            for(int level = context.minLevel; level <= 0; level++)
            {
                context.level = level;
                context.levelFactor = pow(2,-level);
                std::cerr << "  level: " << level << ", factor: " << context.levelFactor << std::endl;
                std::cerr << "  seedCache size: " << seedCache.size() << std::endl;
                double dSepEffective = context.dSep * context.levelFactor;
                for(auto sl: ret.streamlines)
                    extend(context,sl);
                int slStart = 0;
                std::vector<Flow> keptSeeds;
                for(auto seed: seedCache)
                {
                    int sli = slStart;
                    while(sli < ret.streamlines.size())
                    {
                        auto sl = ret.streamlines[sli];
                        sli++;
                        //this looks for candidate seed points as described in JL paper fig 3
                        for(int pn = 0; pn < sl.points.size(); pn += m_iSteps*context.levelFactor)
                        {
                            float dir = sl.points[pn].direction;
                            // check either side (perpendicularly) of the streamline for new seed points.
                            for(int k = 0; k < 2; k++)
                            {
                                auto newSeed = positionFromDistanceCourse(sl.points[pn].position,dSepEffective,dir+1.57079633+k*M_PI);
                                Flow newSeedFlow = field.getFlow(newSeed);
                                if(isPointGood(context,newSeedFlow,dSepEffective))
                                {
                                    Streamline newSl(newSeedFlow,context.level);
                                    extend(context,newSl);
                                    if(newSl.points.size() > 2)
                                        addStreamline(context,newSl);
                                }
                            }
                        }
                    }
                    slStart = ret.streamlines.size();
                    if(isPointGood(context,seed,dSepEffective))
                    {
                        Streamline newSl(seed,context.level);
                        extend(context,newSl);
                        if(newSl.points.size() > 2)
                            addStreamline(context, newSl);
                        else
                            keptSeeds.push_back(seed);
                    }
                    else
                        if(isPointGood(context,seed,context.dSep))
                            keptSeeds.push_back(seed);
                }
                seedCache = keptSeeds;
            }
            ret.dSep = context.dSep;
            ret.iSteps = m_iSteps;
            return ret;
        }
        
    private:
        /// Used in determining the base grid size. Factor is applied to field density when determining dSep parameter.
        float m_seperationFactor;
        
        /// Factor used to calculate dTest relative to dSep.
        float m_testFactor;
        
        /// Number of intermediate steps to use when calculating next point along the streamline.
        int m_iSteps;
        
        /// Maximum dSep relative to data extent when using overview layers.
        float m_dSepMaxFactor;
        
        /// Minimum magnitude to use when extending a streamline. If magnitude fall below this value, it is considered null.
        float m_minMag;
        
        struct PointsGridItem
        {
            LatLong point;
            int streamIndex;
        };
        
        struct GenerateContext
        {
            GenerateContext(FlowFieldTimeStep &f, StreamlineSet &streamlineset):field(f),streamlineset(streamlineset){}

            FlowFieldTimeStep &field;
            double density;
            double dSep;
            double dTest;
            
            LatLong pointsGridCellSpacing;
            
            std::map<int,std::map<int,std::vector<PointsGridItem> > > pointsGrid;
        
            /// As the latitude varies in the points grid, its coverage relative to dSep varies,
            /// this factor is kept track of in the following dict.
            std::map<int,double> pointsGridCellWidthFactors;
            
            int minLevel;
            int level;
            int levelFactor;
            
            StreamlineSet &streamlineset;
        };

        /**
         Extends a streamline in both directions from a seed point.
        
        Args:
            stream: Streamline initialized with a seed point.
        */
        bool extend(GenerateContext &context, Streamline &stream)
        {
            bool ok = true;
            bool extended = false;
            while(ok)
            {
                ok = step(context, stream,1);
                if(ok)
                    extended = true;
            }
            ok = true;
            while(ok)
            {
                ok = step(context,stream,-1);
                if(ok)
                    extended = true;
            }
            return extended;
        }

        /**
         Attempts to add a point to the streamline.
        
        Args:
            sl: Streamline to be extended.
            direction (int): 1 if in the flow direction, -1 if in the direction oposite the flow.
            
        Returns:
            True if step was succefull, False otherwise.
        */
        bool step(GenerateContext &context, Streamline &sl, int direction)
        {
            bool ok = true;
            StreamlinePoint p0;
            if(direction > 0)
                p0 = sl.points.back();
            else
                p0 = sl.points.front();
            Flow pLast = context.field.getFlow(p0.position);
            if(pLast.valid)
            {
                std::vector<Flow> steps;
                while(pLast.valid && steps.size() < m_iSteps*context.levelFactor)
                {
                    std::vector<Flow> iSteps = context.field.transport(pLast.position,context.dSep*direction,m_iSteps,m_minMag);
                    if(iSteps.size() == m_iSteps)
                    {
                        for(auto ip: iSteps)
                            if(isPointGood(context, ip, context.dTest*context.levelFactor, sl.index))
                                pLast = ip;
                            else
                            {
                                pLast = Flow();
                                break;
                            }
                        if(pLast.valid)
                            steps.insert(steps.end(),iSteps.begin(),iSteps.end());
                    }
                    else
                        pLast = Flow();
                    if(pLast.valid)
                        if(!isStreamPointGood(context,sl,pLast.position, steps))
                            pLast = Flow();
                }
                if(pLast.valid)
                    for(auto ip: steps)
                    {
                        StreamlinePoint sp(ip);
                        sp.level = context.level;
                        if(sl.index != -1)
                            addPoint(context, sp, sl.index);
                        sl.add(sp,direction==1);
                    }
                else
                    ok = false;
            }
            else
                ok = false;
            return ok;
        }
        
        
        struct PointsGridIndex
        {
            int x,y;
        };
        
        /**
        Checks if point is acceptable.
        
        Args:
            p (Point): Point to be checked.
            sep (float): Distance in meters point should be from other points.
            streamIndex (int):  
        */
        bool isPointGood(GenerateContext &context, Flow const &p, float sep, int streamIndex = -1)
        {
            if(!p.valid)
                return false;
            PointsGridIndex index = getPointsGridIndex(context, p.position);
            for(int i = index.y - context.levelFactor; i <= index.y+context.levelFactor; i++)
            {
                auto pgIterator = context.pointsGrid.find(i);
                if(pgIterator != context.pointsGrid.end())
                {
                    int lFactor = ceil(context.levelFactor*context.pointsGridCellWidthFactors[i]);
                    for(int j = index.x-lFactor; j <= index.x+lFactor; j++)
                    {
                        auto jIterator = pgIterator->second.find(j);
                        if (jIterator != pgIterator->second.end())
                        {
                            for(auto gp: jIterator->second)
                                if(gp.streamIndex != streamIndex && distanceCourseFromRad(gp.point, p.position).distance < sep)
                                    return false;
                        }
                    }
                }
            }
            return true;
        }
        
        bool isStreamPointGood(GenerateContext &context, Streamline &sl, LatLong const &p, std::vector<Flow> const &extraPoints)
        {
            if(extraPoints.size() == context.levelFactor*m_iSteps)
                for(int i = 0; i < sl.points.size(); i+=m_iSteps*context.levelFactor)
                    if(distanceCourseFromRad(sl.points[i].position,p).distance < context.dTest)
                        return false;
            for(int i = 0; i < sl.points.size(); i+=m_iSteps)
                if(distanceCourseFromRad(sl.points[i].position,p).distance < context.dTest)
                    return false;
            for(int i = 0; i<extraPoints.size()-m_iSteps; i+=m_iSteps)
                if(distanceCourseFromRad(extraPoints[i].position,p).distance < context.dTest)
                    return false;
            return true;
        }
        
        /**
        Adds a point to the points grid, which is used to keep track of points not available
        to be added to streamlines. The points grid is implemented as a dictionary of dictionaries
        to allow it to be sparse, in other words, cells are created on demand. Each cell contains a
        list of points.
        
        Args:
            p (Point): Point to be added.
            streamIndex (int): Index of the streamline to which point p belongs.
        */
        void addPoint(GenerateContext &context, StreamlinePoint const &p, int streamIndex)
        {
            PointsGridIndex index = getPointsGridIndex(context,p.position);
            if(context.pointsGrid.find(index.y) == context.pointsGrid.end())
            {
                LatLong p0(index.y*context.pointsGridCellSpacing.latitude+context.field.bounds().min().latitude,0.0);
                LatLong p1 = positionFromDistanceCourse(p0, context.dSep, 1.5708); // 90 degrees in radians
                //calculates a factor to know how much dsep spans in terms of grid spacing
                context.pointsGridCellWidthFactors[index.y] = p1.longitude/context.pointsGridCellSpacing.longitude;
            }
            PointsGridItem pgi;
            pgi.point = p.position;
            pgi.streamIndex = streamIndex;
            context.pointsGrid[index.y][index.x].push_back(pgi);
        }

        /**
        Calculates the index into the points grid. The points grid is used to make sure that
        a new points isn't too close to an existing point.
        */
        PointsGridIndex getPointsGridIndex(GenerateContext const &context, LatLong const &p)
        {
            PointsGridIndex ret;
            ret.x = floor((p.longitude-context.field.bounds().min().longitude)/context.pointsGridCellSpacing.longitude);
            ret.y = floor((p.latitude-context.field.bounds().min().latitude)/context.pointsGridCellSpacing.latitude);
            return ret;
        }
        
        ///Adds a streamline to the collection of streamlines.
        void addStreamline(GenerateContext &context, Streamline &sl)
        {
            if(sl.index == -1)
            {
                sl.index = context.streamlineset.streamlines.size();
                context.streamlineset.streamlines.push_back(sl);
                context.streamlineset.bounds.add(sl.bounds);
                for(int i = 0; i < sl.points.size(); i+=m_iSteps)
                    addPoint(context, sl.points[i], sl.index);
            }
        }
    };
}

#endif
