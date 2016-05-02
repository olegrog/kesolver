#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <stdexcept>
#include <limits>

#include "Printer.hpp"

#include "logger/logger.hpp"

Printer::Printer(const PropertyTree& tree)
{
    save_func = tree.isMember("savefunc") ? tree["savefunc"].asBool() : false;

    LOG(INFO) << "save_func = " << save_func;

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

void Printer::saveMacroParams(int i,
                              Mesh& mesh,
                              const Gas& gas)
{
    std::ostringstream ss;
    for (unsigned int j = 0; j < mesh.getMyCells().size(); j++)
        ss << mesh.getMyCellIndexes()[j] << ' ' 
           << gas.macro(mesh.getMyCells()[j]->f().f()) << ' '
           << mesh.getMyCells()[j]->getExcludedPoints()/gas.density(mesh.getMyCells()[j]->f().f())
           << std::endl;

    char macroFileName[256]; 
    sprintf(macroFileName, (dirname + filename).c_str(), i);
    std::ofstream macroParamsFile;

    int size = mesh.getSize();
    int rank = mesh.getRank();

    for(int k = 0; k < size; k++) {
        if (rank == k) {
            if (rank == 0)
                macroParamsFile.open(macroFileName, std::ios::out);
            else
                macroParamsFile.open(macroFileName, std::ios::out | std::ios::app);

            macroParamsFile << ss.str();

            macroParamsFile.close();
        }
        mesh.barrier();
    }
}

void Printer::saveSpeedFunction(int i,
                                Mesh& mesh,
                                const Gas& gas)
{
    char functionFileName[256];
    sprintf(functionFileName, (dirname + functionfilename).c_str(), i);
    std::ofstream function;
    std::string name("Printer::saveFunction");

    int size = mesh.getSize();
    int rank = mesh.getRank();

    for (int k = 0; k < size; k++) {
        if (rank == k) {
            if (rank == 0)
                function.open(functionFileName, std::ios::out);
            else
                function.open(functionFileName, std::ios::out | std::ios::app);
            for (size_t i = 0; i < mesh.getMyCells().size(); i++) {
                int index = mesh.getMyCellIndexes()[i];
                function.write(reinterpret_cast<const char*>(&index), sizeof(int));
                SpeedFunction& f = mesh.getMyCells()[i]->f();
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
        mesh.barrier();
    }
}

void Printer::print(int i,
                    Mesh& mesh,
                    const Gas& gas)
{
    if ( i % save_macro_point == 0 ) { 
        LOG(INFO) << "save_macro";
        saveMacroParams(i, mesh, gas);
    }
    if ( ( save_func ) && ( i % save_func_freq == 0 ) ) {
        LOG(INFO) << "save_func";
        saveSpeedFunction(i, mesh, gas);
    }
}

