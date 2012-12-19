#ifndef PRINTER_H
#define PRINTER_H

#include <fstream>

#include "Gas.hpp"
#include "mesh/unstruct/Polygon.hpp"

#include "property_tree/property_tree.hpp"

class Printer {
	public:
		void print(int i, const std::vector<Polygon*>& spacemesh,
				const std::vector<int>& mypolys,
				const Gas& gas,
				int size, int rank);

		Printer(const PropertyTree& tree);

	private:
		std::string dirname, filename, functionfilename;
        bool save_func;
        int save_macro_point, save_func_freq;

		void saveMacroParams(int i, const std::vector<Polygon*>& spacemesh,
				const std::vector<int>& mypolys,
				const Gas& gas,
				int size, int rank);

		void saveSpeedFunction(int i, const std::vector<Polygon*>& spacemesh,
				const std::vector<int>& mypolys,
				int size, int rank);


};

#endif /*PRINTER_H*/
