#ifndef S111_STREAMLINES_H
#define S111_STREAMLINES_H

#include "H5Cpp.h"
#include <vector>
#include <iostream>

namespace s111
{
    struct Flow
    {
        Flow(float mag, float dir):magnitude(mag),direction(dir)
        {}
            
        float magnitude; 
        float direction;
    };
    
    class FlowField
    {
    public:
        FlowField(){}
        
        void addFile(std::string const &fname)
        {
            H5::H5File f(fname, H5F_ACC_RDONLY);
            m_files.push_back(f);

            auto scGroup = f.openGroup("/SurfaceCurrent/SurfaceCurrent.01");
            scGroup.
            //m_datasets.push_back(f.openDataSet("SurfaceCurrent"));
        }
    private:
        std::vector<H5::H5File> m_files;
        std::vector<H5::DataSet> m_datasets;
    };
}

#endif
