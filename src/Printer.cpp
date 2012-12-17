#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <stdexcept>
#include <limits>

#include "mpi.h"

#include "Printer.hpp"

Printer::Printer(const PropertyTree& tree)
{
    save_func = tree.isMember("savefunc") ? tree["savefunc"].asBool() : false;

    save_func_freq = std::numeric_limits<int>::max();
    if (save_func) {
        functionfilename = tree["savefuncfilename"].asString();
        save_func_freq   = tree["savefuncfreq"].asInt();
    }   

    dirname =  tree["dir"].asString();
    filename = tree["file"].asString();

    save_macro_point = tree["savemacro"].asInt();

    int status;
    status = mkdir( dirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

    std::cout << "dirname = " << dirname << ' ' << status << std::endl;
}

void Printer::saveMacroParams(int i, const std::vector<Polygon*>& spacemesh,
                              const std::vector<int>& mypolys, const Gas& gas, 
                              int size, int rank)
{
    std::ostringstream ss;
    for (unsigned int j = 0; j < mypolys.size(); j++)
        ss << mypolys[j] << ' ' << gas.macro(spacemesh[mypolys[j]]->f().f()) << std::endl;

    char macroFileName[256]; 
    sprintf(macroFileName, (dirname + filename).c_str(), i);
    std::ofstream macroParamsFile;

    for(int k = 0; k < size; k++) {
        if (rank == k) {
            if (rank == 0)
                macroParamsFile.open(macroFileName, std::ios::out);
            else
                macroParamsFile.open(macroFileName, std::ios::out | std::ios::app);

            macroParamsFile << ss.str();

            macroParamsFile.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void Printer::saveSpeedFunction(int i, const std::vector<Polygon*>& spacemesh,
                                const std::vector<int>& mypolys, int size, int rank)
{
    char functionFileName[256];
    sprintf(functionFileName, (dirname + functionfilename).c_str(), i);
    std::ofstream function;
    std::string name("Printer::saveFunction");

    for (int k = 0; k < size; k++) {
        if (rank == k) {
            if (rank == 0)
                function.open(functionFileName, std::ios::out);
            else
                function.open(functionFileName, std::ios::out | std::ios::app);
            for (size_t i = 0; i < mypolys.size(); i++) {
                function.write( reinterpret_cast<const char*>(&mypolys[i]), sizeof(int) );
                SpeedFunction& f = spacemesh[mypolys[i]]->f();
                size_t fsize = f.size();
                function.write( reinterpret_cast<const char*>(&fsize), sizeof(int) );
                double* f_p = &(f.f()[0]);
                function.write( reinterpret_cast<char*>(f_p),
                        sizeof(double)*fsize);
                if (function.bad()) {
                    // TODO
                }
            }
            function.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void Printer::print(int i, const std::vector<Polygon*>& spacemesh,
        const std::vector<int>& mypolys, const Gas& gas, int size, int rank)
{
    if ( i % save_macro_point == 0 ) { 
        std::cout << "save_macro" << std::endl;
        saveMacroParams(i, spacemesh, mypolys, gas, size, rank);
    }
    if ( ( save_func ) && ( i % save_func_freq == 0 ) ) {
        std::cout << "save_func" << std::endl;
        saveSpeedFunction(i, spacemesh, mypolys, size, rank);
    }
}

