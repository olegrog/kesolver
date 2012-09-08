#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <map>
#include <vector>
#include <utility>
#include <exception>
#include <stdexcept>

#include <libxml/xmlreader.h>

#include "Loader.hpp"
#include "auxiliary.hpp"

Loader::Loader(const std::string& kinfile_name) {
	std::cout << kinfile_name << std::endl;

	xmlTextReaderPtr reader;
	reader = xmlReaderForFile(kinfile_name.c_str(), NULL, 0);

	if (reader != NULL) {
        while (xmlTextReaderRead(reader) == 1) {
			std::string name = toStr(xmlTextReaderConstName(reader));
			std::cout << name << std::endl;
			if ( name == "nodes" ) 
				readNodes(reader);
			else if ( name == "cells" )
				readElems(reader, cells);
			else if ( name == "facets" )
				readElems(reader, facets);
			else if ( name == "physical" )
				readPhysicalData(reader);
			else if ( xmlTextReaderHasAttributes(reader) ) {
				readData(reader, name);
                if ( name == "gas" )
                    readGas(reader);
            }
		}
		xmlFreeTextReader(reader);
	} else {
		std::cout << "Unable to open " << kinfile_name << std::endl;
	}

	xmlCleanupParser();
	xmlMemoryDump();

	try {
		std::string functionfilename = getData("velfunction", "filename");
        readVelFunction(functionfilename);
	}
	catch (std::invalid_argument) {}

}

const std::string Loader::toStr(const xmlChar* xml_str) {
	if (xml_str != NULL)
		return std::string( (const char*)xml_str );
	else
		return std::string("");
}

void Loader::readData(xmlTextReaderPtr& reader, const std::string& name) {
	attributes_type attributes;
	xmlTextReaderMoveToFirstAttribute(reader);
	do {
		std::string key = toStr(xmlTextReaderConstName(reader));
		std::string value =  toStr(xmlTextReaderConstValue(reader));
		attributes[key] = value;
//		cout << key << " = " << value << endl;
	} while ( xmlTextReaderMoveToNextAttribute(reader) );
	xmlTextReaderMoveToElement(reader);
	data[name] = attributes;
}


void Loader::readNodes(xmlTextReaderPtr& reader) {
	int old_depth = xmlTextReaderDepth(reader);

	int len = 0;
	xmlTextReaderMoveToFirstAttribute(reader);
	do {
		if ( toStr(xmlTextReaderConstName(reader)) == "len" )
			len = strTo<int>( toStr(xmlTextReaderConstValue(reader)) );
	} while ( xmlTextReaderMoveToNextAttribute(reader) );
	xmlTextReaderMoveToElement(reader);

	xmlTextReaderRead(reader);
	if ( ( toStr(xmlTextReaderConstName(reader)) == "#text" ) &&
			( xmlTextReaderDepth(reader) == old_depth+1 ) ) {
		std::string value = toStr(xmlTextReaderConstValue(reader));
		std::istringstream ss(value);
		nodes.resize(len);
		for (int i = 0; i < len; ++i) {
			double x, y, z;
			ss >> x >> y >> z;
			nodes[i] = V3d(x, y, z);
//			cout << V3d(x, y, z) << endl;
		}
	}
	else {	// TODO
	}
	std::cout << "number of nodes is " << nodes.size() << std::endl;
}

void Loader::readElems(xmlTextReaderPtr& reader, std::vector< std::vector<int> >& elems) {
	int old_depth = xmlTextReaderDepth(reader);
	xmlTextReaderRead(reader);
	if ( ( toStr(xmlTextReaderConstName(reader)) == "#text" ) &&
			( xmlTextReaderDepth(reader) == old_depth+1 ) ) {
		std::string value = toStr(xmlTextReaderConstValue(reader));
		std::istringstream ss(value);
		while (ss) {
			ss >> std::ws;
			char line[128];
			ss.getline(line, 128);
			std::string str(line);
			std::istringstream iss(str);
			std::vector<int> data;
			std::copy( std::istream_iterator<int>(iss), std::istream_iterator<int>(),
				std::back_inserter< std::vector<int> >(data) );
//			cout << "elems[].size = " << data.size() << endl;
			if (data.size() > 0)
				elems.push_back(data);
		}
	}
	else { // TODO
	}
	std::cout << "number of elems is " << elems.size() << std::endl;
}

void Loader::readPhysicalData(xmlTextReaderPtr& reader) {
	int old_depth = xmlTextReaderDepth(reader);
	xmlTextReaderRead(reader);
	if ( ( toStr(xmlTextReaderConstName(reader)) == "#text" ) &&
			( xmlTextReaderDepth(reader) == old_depth+1 ) ) {
		std::string value = toStr(xmlTextReaderConstValue(reader));
		std::istringstream ss(value);
		while (ss) {
			ss >> std::ws;
			char line[255];
			ss.getline(line, 255);
			std::string str(line);
			std::istringstream iss(str);
			int key;
			iss >> key;
			std::vector<std::string> data;
			std::copy( std::istream_iterator<std::string>(iss),
                       std::istream_iterator<std::string>(),
				       std::back_inserter< std::vector<std::string> >(data) );
//			cout << "key = " << key << " physical_data[].size = " << data.size() << endl;
			if (data.size() > 0)
				physical_data[key] = data;
		}
	}
	else { // TODO
	}
	std::cout << "number of physical_data is " << physical_data.size() << std::endl;
}

void Loader::readVelFunction(const std::string& functionfilename) {

	std::ifstream function( functionfilename.c_str(), std::ifstream::binary );
	std::string name("Loader::loadFunction");

	size_t size;
	function.read( reinterpret_cast<char*>(&size), sizeof(size_t) );
	velfunction.resize(size);

	for (size_t i = 0; i < size; ++i) {
		int j, k;
		function.read( reinterpret_cast<char*>(&j), sizeof(int) );
		function.read( reinterpret_cast<char*>(&k), sizeof(int) );
		velfunction[j].resize(k);
		function.read( reinterpret_cast<char*>(&velfunction[j][0]), sizeof(double)*k );
	}

}

void Loader::readGas(xmlTextReaderPtr& reader) {
	int old_depth = xmlTextReaderDepth(reader);
	xmlTextReaderRead(reader);
	if ( ( toStr(xmlTextReaderConstName(reader)) == "#text" ) &&
			( xmlTextReaderDepth(reader) == old_depth+1 ) ) {
		std::string value = toStr(xmlTextReaderConstValue(reader));
		std::istringstream ss(value);
		while (ss) {
			ss >> std::ws;
			char line[64];
			ss.getline(line, 64);
			std::string str(line);
			std::istringstream iss(str);
			std::vector<std::string> data;
			std::copy( std::istream_iterator<std::string>(iss), 
                       std::istream_iterator<std::string>(),
				       std::back_inserter< std::vector<std::string> >(data) );
//			cout << "key = " << key << " physical_data[].size = " << data.size() << endl;
			if (data.size() > 0)
				gas_data.push_back(data);
		}
	}
	else { // TODO
	}
	std::cout << "gas_data.size() = " << gas_data.size() << std::endl;
}

const std::string Loader::getData(const std::string& item, const std::string& name) const {
    xmlmap::const_iterator item_p = data.find(item);
	if (item_p == data.end())
		throw std::invalid_argument("wrong item: " + item);
	else  {
        const attributes_type& values = (*item_p).second;
        attributes_type::const_iterator value_p = values.find(name);
		if ( value_p == values.end() ) 
			throw std::invalid_argument("wrong name: " + name);
		else {
			std::string s = (*value_p).second;
			return s;
        }
    }
}


