#include <iostream>
#include "mpi.h"
#include "DataExchanger.hpp"
#include "Constructors.hpp"

void DataExchanger::swap()
{
   	MPI_Startall(all_size, request);
    MPI_Waitall(all_size, request, statuses);
}


void DataExchanger::prepareSendRecv(std::vector<Polygon*>& cells, int rank)
{
	for(size_t i = 0; i < cells.size(); i++)
        if (cells[i]->getRank() == rank)
            for(size_t j = 0; j < cells[i]->getNeigbors().size(); j++) {
                int m = cells[cells[i]->getNeigbors()[j]]->rank;
                if (m != rank) {
                    _recv_process.push_back(m);
                    _send_process.push_back(m);
                    _recv_tetr.push_back(cells[i]->getNeigbors()[j]);
                    _send_tetr.push_back(i);
                }
	        }
}

void DataExchanger::prepareSendRecv2(std::vector<Polygon*>& cells, int rank)
{
	for(size_t i = 0; i < cells.size(); i++)
        if (cells[i]->getRank() == rank)
            for(size_t j = 0; j < cells[i]->getNeigbors().size(); j++) {
                int nr = cells[i]->getNeigbors()[j];
                int m = cells[nr]->rank;
                if (m != rank) {
                    _recv_process.push_back(m);
                    _send_process.push_back(m);
                    _recv_tetr.push_back(nr);
                    _send_tetr.push_back(i);
                }
                for(size_t k = 0; k < cells[nr]->getNeigbors().size(); k++) {
                    int nr2 = cells[nr]->getNeigbors()[k];
                    if (cells[nr2]->rank != rank) {
                        _recv_process.push_back(cells[nr2]->rank);
                        _send_process.push_back(cells[nr2]->rank);
                        _recv_tetr.push_back(nr2);
                        _send_tetr.push_back(i);
                    }
                }
            }
}

void DataExchanger::cleanSendRecv()
{
    for(size_t i = 0; i < _send_tetr.size(); i++) {
        int k = 0;
        for(size_t j = 0; j < i; j++)
			if(_send_tetr[i] == _send_tetr[j] && _send_process[i] == _send_process[j])
				k++;
        if(k == 0) {
			send_tetr.push_back(_send_tetr[i]);
			send_process.push_back(_send_process[i]);
		}
    }
    _send_tetr.clear();
    _send_process.clear();
    for(size_t i = 0; i < _recv_tetr.size(); i++) {
        int k = 0;
        for(size_t j = 0; j < i; j++)
			if(_recv_tetr[i] == _recv_tetr[j] && _recv_process[i] == _recv_process[j])
				k++;
        if(k == 0) {
			recv_tetr.push_back(_recv_tetr[i]);
			recv_process.push_back(_recv_process[i]);
		}
    }
    _recv_tetr.clear();
    _recv_process.clear();

	std::cout << "send = " << send_tetr.size() << ' ' << 
                 "recv = " << recv_tetr.size() << std::endl;
}

void DataExchanger::mpiInit(std::vector<Polygon*>& cells)
{
    size_t rsize = recv_process.size();
    size_t ssize = send_process.size();

    all_size = rsize + ssize;
	request  = new MPI_Request[all_size];
    statuses = new MPI_Status[all_size];

    int j = 0;

    for(size_t i = 0; i < recv_tetr.size(); i++) {
        int index   = recv_tetr[i];
        double* f   = &(cells[recv_tetr[i]]->f().f()[0]);
        size_t size = cells[recv_tetr[i]]->f().size();
        MPI_Recv_init(f, size, MPI_DOUBLE,
                      recv_process[i], index,
                      MPI_COMM_WORLD, &request[j++]);
    }
    for(size_t i = 0; i < send_tetr.size(); i++) {
        int index   = send_tetr[i];
        double* f   = &(cells[send_tetr[i]]->f().f()[0]);
        size_t size = cells[send_tetr[i]]->f().size();
        MPI_Send_init(f, size, MPI_DOUBLE,
                      send_process[i], index,
                      MPI_COMM_WORLD, &request[j++]);
    }
}

void DataExchanger::init(std::vector<Polygon*>& cells, int rank)
{
	prepareSendRecv(cells, rank);
	cleanSendRecv();
}

void DataExchanger::init2(std::vector<Polygon*>& cells, int rank)
{
	prepareSendRecv2(cells, rank);
	cleanSendRecv();
}

