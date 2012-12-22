#ifndef PRINTER_H
#define PRINTER_H

#include <fstream>

#include "Gas.hpp"
#include "MeshMpi.hpp"
#include "mesh/unstruct/Polygon.hpp"

#include "property_tree/property_tree.hpp"

class Printer {
	public:
		void print(int i,
                   MeshMpi& mesh, 
				   const Gas& gas);

		Printer(const PropertyTree& tree);

	private:
		std::string dirname, filename, functionfilename;
        bool save_func;
        int save_macro_point, save_func_freq;

		void saveMacroParams(int i,
                             MeshMpi& mesh, 
				             const Gas& gas);

		void saveSpeedFunction(int i,
                               MeshMpi& mesh, 
				               const Gas& gas);

};

#endif /*PRINTER_H*/
