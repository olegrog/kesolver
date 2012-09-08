#ifndef LOADER_HPP
#define LOADER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include <libxml/xmlreader.h>

#include "v.hpp"
#include "auxiliary.hpp"

class Loader {
	public:
		Loader(const std::string& kshfile_name);

		const std::vector<V3d>& getNodes() const { return nodes; }
		const std::vector< std::vector<int> >& getCells() const { return cells;	}
		const std::vector< std::vector<int> >& getFacets() const { return facets; }
		const std::vector<std::string>& getPhysicalData(int i) const {
			return (*(physical_data.find(i))).second;
		}
		const std::vector< std::vector<std::string> >& getGasData() const { return gas_data; }

		template <typename T>
		T getData(const std::string& item, const std::string& name) const {
			return strTo<T>( getData(item, name) );
		}
		const std::string getData(const std::string& item, const std::string& name) const;

        const std::vector<double>& getVelFunction(int i) const {
            return velfunction[i];
        }

	private:
		std::vector<V3d> nodes;
		std::vector< std::vector<int> > cells, facets;
		std::map<int, std::vector<std::string> > physical_data;
		std::vector< std::vector<double> > velfunction;
		std::vector< std::vector<std::string> > gas_data;

		typedef std::map< std::string, std::string > attributes_type;
		typedef std::map< std::string, attributes_type > xmlmap;
		xmlmap data;

		void readData(xmlTextReaderPtr& reader, const std::string& name);
		void readNodes(xmlTextReaderPtr& reader);
		void readElems(xmlTextReaderPtr& reader, std::vector< std::vector<int> >& elems);
		void readPhysicalData(xmlTextReaderPtr& reader);
		void readGas(xmlTextReaderPtr& reader);

		void readVelFunction(const std::string& functionfilename);

		const std::string toStr(const xmlChar* xml_str);

};

#endif /*LOADER_HPP*/
