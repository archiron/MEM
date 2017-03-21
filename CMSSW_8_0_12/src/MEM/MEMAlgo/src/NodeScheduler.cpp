# include "MEM/MEMAlgo/interface/NodeScheduler.h"

void NodeScheduler::initNodeScheduler( RunConfig cfg, int mpi_rank ) {
  
  config = cfg;
  //cfg->print();
}

/*
 * @brief Configure Vegas state for integration with 
 *        basic data:
 *        (dim, nbrOfCalls, limits, ...)
 * @param vegasState
 */

