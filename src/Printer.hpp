#ifndef PRINTER_H
#define PRINTER_H

#include <fstream>

#include "Gas.hpp"
#include "Mesh.hpp"
#include "mesh/unstruct/Polygon.hpp"

#include "property_tree/property_tree.hpp"

class Printer {
	public:
		void print(int i,
                   Mesh& mesh, 
				   const Gas& gas);

		Printer(const PropertyTree& tree);

	private:
		std::string dirname, filename, functionfilename;
        bool save_func;
        int save_macro_point, save_func_freq;

		void saveMacroParams(int i,
                             Mesh& mesh, 
				             const Gas& gas);

		void saveSpeedFunction(int i,
                               Mesh& mesh, 
				               const Gas& gas);

};

#endif /*PRINTER_H*/
