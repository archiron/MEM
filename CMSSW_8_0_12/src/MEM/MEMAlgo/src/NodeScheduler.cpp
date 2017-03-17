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
/*void NodeScheduler::configOCLVegasState( ocl_vegas_state_t *vegasState, size_t  nbrOfRequestedPoints ) {

  if (FCTDims > MaxNumberOfDIMs) {
     fprintf(stderr, "ABORT: too large number of dimensions [dim=%d]\n", FCTDims);
     return;
  }
  vegasState->verbose = 1;
  setVegasDefault( vegasState, FCTDims );
  vegasState->nbrOfWarmUpFctCalls    = VEGASWarmUpNbrOfPoints; 
  vegasState->nbrOfRequestedFctCalls = nbrOfRequestedPoints;  
  vegasState->maxNbrOfIIterations    = NbrMaxVegasChi2Iterations; 
  int d;
  for ( d=0; d<FCTDims; d++) {
    vegasState->xl[d] = XLowerDomain;
    vegasState->xu[d] = XUpperDomain;
  }
  vegasState->currentIt = 0;
  
  //
  vegasState->workGroupSize = -1;
  
  //
  vegasState->integralInDomainPerIt=0;
  vegasState->varianceInDomainPerIt=0;
  vegasState->CumulatedIntegral = 0;
  vegasState->CumulatedSigma = 0;
  
  // Debug
  vegasState->liForTest[0]=0;
  vegasState->liForTest[1]=0;
  vegasState->liForTest[2]=0;
  vegasState->liForTest[3]=0;
}
*/

/*
 * @brief Test if the configuration is well suited for Configure Vegas state for integration with 
 *        basic data:
 *        (dim, nbrOfCalls, limits, ...)
 * @param vegasState
 */
/*void NodeScheduler::validateOCLVegasState( ocl_vegas_state_t *vegasState, const RunConfig *cfg ) 
   {
  

  if ( vegasState->nbrBoxesInDomain > MAXNbrOfVegasBoxes) { 
    throw new domain_error("Number of boxes too large");
  }
   
}
*/