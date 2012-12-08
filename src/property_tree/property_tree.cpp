#include <iostream>
#include <fstream>

#include "property_tree.hpp"

PropertyTree_jsoncpp::PropertyTree_jsoncpp(const char* conffilename)
{
    std::ifstream input(conffilename);
    Json::Reader reader;

    bool parsing_successful = reader.parse(input, jsonvalue);
    if (!parsing_successful)
    {
        std::cout << "Failed to parse configuration file "
                  << reader.getFormatedErrorMessages() << std::endl;
        // TODO
    } 
}

