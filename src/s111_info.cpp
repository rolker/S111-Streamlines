#include "s111_streamlines.h"

int main(int argc, char *argv[])
{
    s111::FlowField flowField;
    for(int i = 1; i < argc; i++)
        flowField.addFile(argv[i]);
    std::cerr << "bounds: " << flowField.bounds() << std::endl;
    
    s111::JobardLefer jl;
    s111::FlowFieldTimeStep ffts(flowField,*(flowField.timePoints().begin()));
    auto streamlines = jl.generate(ffts);
    
    std::cerr << "generated " << streamlines.streamlines.size() << " streamlines." << std::endl;
    std::cerr << "bounds: " << streamlines.bounds << std::endl;
    
    std::cout << "{\n";
    std::cout << "  \"type\": \"FeatureCollection\",\n";
    std::cout << "  \"bbox\": [" << streamlines.bounds.min().longitude_degrees() << ", " << streamlines.bounds.min().latitude_degrees() << ", " << streamlines.bounds.max().longitude_degrees() << ", " << streamlines.bounds.max().latitude_degrees() << "],\n";
    
    std::cout << "  \"features\": [\n";
    
    for(int i = 0; i < streamlines.streamlines.size(); i++)
    {
        std::cout << "    {\n";
        std::cout << "      \"type\": \"Feature\",\n";
        std::cout << "      \"geometry\": {\n";
        std::cout << "        \"type\": \"LineString\",\n";
        std::cout << "        \"coordinates\": [\n";
        for(int j = 0; j < streamlines.streamlines[i].points.size(); j++)
        {
            std::cout << "          [" << streamlines.streamlines[i].points[j].position.longitude_degrees() << ", " << streamlines.streamlines[i].points[j].position.latitude_degrees() << "]";
            if(j < streamlines.streamlines[i].points.size()-1)
                std::cout << ",";
            std::cout << "\n";
        }
        std::cout << "        ]\n";
        std::cout << "      },\n";
        std::cout << "      \"properties\": {\n";
        
        std::cout << "        \"index\": " << streamlines.streamlines[i].index << ",\n";
        std::cout << "        \"streamline_level\": " << streamlines.streamlines[i].level << ",\n";
        std::cout << "        \"seed_index\": " << streamlines.streamlines[i].seedIndex << ",\n";

        std::cout << "        \"point_levels\": [";
        for(int j = 0; j < streamlines.streamlines[i].points.size(); j++)
        {
            std::cout << streamlines.streamlines[i].points[j].level;
            if(j < streamlines.streamlines[i].points.size()-1)
                std::cout << ",";
        }
        std::cout << "],\n";

        std::cout << "        \"magnitudes\": [";
        for(int j = 0; j < streamlines.streamlines[i].points.size(); j++)
        {
            std::cout << streamlines.streamlines[i].points[j].magnitude;
            if(j < streamlines.streamlines[i].points.size()-1)
                std::cout << ",";
        }
        std::cout << "],\n";

        std::cout << "        \"directions\": [";
        for(int j = 0; j < streamlines.streamlines[i].points.size(); j++)
        {
            std::cout << streamlines.streamlines[i].points[j].direction*180.0/M_PI;
            if(j < streamlines.streamlines[i].points.size()-1)
                std::cout << ",";
        }
        std::cout << "],\n";
        
        
        std::cout << "        \"dSep\": " << streamlines.dSep << ",\n";
        std::cout << "        \"iSteps\": " << streamlines.iSteps << "\n";
        std::cout << "        }\n";
        
        std::cout << "    }";
        if(i < streamlines.streamlines.size()-1)
            std::cout << ",";
        std::cout << "\n";
    }
    
    std::cout << "  ]\n";
    
    std::cout << "}\n";
    
    
    
    return 0;
}
