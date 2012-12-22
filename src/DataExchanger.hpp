#ifndef DATAEXCHANGER_H
#define DATAEXCHANGER_H

#include <vector>

#include "mpi.h"

#include "mesh/unstruct/Polygon.hpp"

class DataExchanger {
    private:
        std::vector<int> recv_tetr;
        std::vector<int> send_tetr;
        std::vector<int> recv_process;
        std::vector<int> send_process;
        std::vector<int> _recv_tetr;
        std::vector<int> _send_tetr;
        std::vector<int> _recv_process;
        std::vector<int> _send_process;

        size_t all_size;
        MPI_Request *request;
        MPI_Status  *statuses;

        void prepareSendRecv(std::vector<Polygon*>& cells,
                             int rank);
        void prepareSendRecv2(std::vector<Polygon*>& cells, 
                              int rank);

        void cleanSendRecv();

    public:
        void swap();

        void init(std::vector<Polygon*>& cells,
                  int rank);
        void init2(std::vector<Polygon*>& cells,
                  int rank);
        void mpiInit(std::vector<Polygon*>& cells);

        const std::vector<int>& getExchangeCells() const { return recv_tetr; }
};

#endif /*DATAEXCHANGER_H*/
