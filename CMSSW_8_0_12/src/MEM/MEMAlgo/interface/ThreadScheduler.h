/* 
 * File:   deviceScheduler.h
 * Author: grasseau
 *
 * Created on 27 août 2014, 23:25
 */

#ifndef THREADSCHEDULER_H
#define	THREADSCHEDULER_H

# include "MEM/MEMAlgo/interface/NodeScheduler.h"
# include "MEM/MEMAlgo/interface/MGIntegration.h"

# if (USE_GSL_LIB == 1)
# include "gsl/gsl_math.h"
# include "gsl/gsl_monte_vegas.h"
# include "gsl/gsl_pow_int.h"
# endif

class ThreadScheduler: public NodeScheduler {

 // GG XXX
 public:
 
  int totalNbrOfQueueAndProcess;

  long int remainingTasks;
  
  MGIntegration *integration;
  
  // GSL
//  gsl_monte_function gslFctDescr;
//  gsl_monte_vegas_state *state;
  gsl_monte_function fctDescr_ttH;
  gsl_monte_vegas_state *state_ttH;
  gsl_monte_function fctDescr_ttH_miss;
  gsl_monte_vegas_state *state_ttH_miss;
  gsl_rng *rng;
  gsl_rng *saveRNGSeed;
  gsl_monte_function fctDescr_ttZ;
  gsl_monte_vegas_state *state_ttZ;
  gsl_monte_function fctDescr_ttZ_miss;
  gsl_monte_vegas_state *state_ttZ_miss;
 
  ThreadScheduler();
  ~ThreadScheduler();
 
  virtual void initNodeScheduler( RunConfig config, int mpi_rank );
  

//  void initAllDevices ( ) {};

//  void cleanNodeScheduler( );

  void runNodeScheduler ( eventList_t evMsgList, int nbrOfEvents);

  void completeNodeScheduler( eventList_t evList );

  void terminateNodeScheduler( );

private: 
  /* GG XXX Not used up to now.
   *   But can be useful to develop p-thread implementation  
  // Integration set-up of function parameters and vega state
  // with event-independaent data
  void setupFctParameters( const RunConfig *config,  ocl_fct_parameters_t *fctParameters );
  void setupVegasStates( ocl_vegas_state_t *vegasStates);

  // Complete the set-up of function parameters and vegas states
  // with event-dependaent data
  void updateFctParametersWithIMsg( const IntegrationMsg_t *evMsg,  ocl_fct_parameters_t *fctParameters, bool fctType );
  void updateVegasStateWithIMsg( const IntegrationMsg_t *evMsg, ocl_vegas_state_t *vegasState, bool fctType);

  // Save integral computation in IntegrationMsg_t
  void saveIntegralInIMsg( const ocl_vegas_state_t *vegasState, const QueueStatus_t *qStatus, IntegrationMsg_t *iMsg, bool iType);          
*/

};

#ifdef	__cplusplus
extern "C" {
#endif
  

#ifdef	__cplusplus
}
#endif

#endif	/* THREADSCHEDULER_H */

