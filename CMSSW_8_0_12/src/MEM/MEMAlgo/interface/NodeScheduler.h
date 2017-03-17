/* 
 * File:   deviceScheduler.h
 * Author: grasseau
 *
 * Created on 27 août 2014, 23:25
 */

#ifndef NODESCHEDULER_H
#define	NODESCHEDULER_H

# include <fstream>
# include <stdexcept>

# include "MEM/MEMAlgo/interface/config.h"
# include "MEMDataFormats/MEMDataFormats/interface/RunConfig.h"
# include "MEM/MEMAlgo/interface/MGIntegration.h"

//# include "Utilities/timing.h"

const int nbrOfEventsPerBlock = NbrEventPerBlock;
typedef IntegrationMsg_t eventList_t[nbrOfEventsPerBlock];

// Queue Status
typedef struct {
  
  const int ComputationCompleted = 1;
  const int ComputationStarted   = 0;
  const int ComputationNotActive = -1;
/*  static const int ComputationCompleted = 1;
  static const int ComputationStarted   = 0;
  static const int ComputationNotActive = -1;*/

  //ticks_t tickStart;
  //ticks_t tickEnd;
  int     platform;
  int     device; 
  // Index of the processed events (startEventIndex, nbrOfEvents) in 
  // the event list vegasState
  int startEventIndex;
  int nbrOfEvents;
  int     status;
} QueueStatus_t;

class NodeScheduler {


 // GG XXX
 public:
  typedef enum {Vegas, EvalOnGrid} ComputationType_t;
  struct timespec delay;
  
  RunConfig config;

  int nbrOfQueues;
  // 
  //ocl_vegas_state_t    hostState;
  //ocl_vegas_state_t    *vegasStates;
  //ocl_fct_parameters_t *fctParameters;

  // Devices Status
  //
  QueueStatus_t *queueStatus;
  
  // For Grid evaluation
  unsigned long int nbrPointsPerAxe;
  ComputationType_t computationType;
  
  // One log file per process
  // GG XXX should be set in MPI Scheduler
  ofstream  logStream_; 
  
  NodeScheduler() { nbrPointsPerAxe = 2; computationType = Vegas; 
                    delay.tv_sec = 0L;  //.tv_sec
                    delay.tv_sec = 50L; // tv_nsec 
                  };
                  
  ~NodeScheduler() {};

  virtual void initNodeScheduler( RunConfig config, int mpi_rank );
  /*
  void initOneDevice( cl_command_queue dQueue, cl_kernel k, 
                    ocl_vegas_state_t *hostState, 
                    ocl_fct_parameters_t *hostFctParams,
                    int queueID );
  */
//  virtual void initAllDevices ( ) = 0;

//  virtual void cleanNodeScheduler( )= 0;

  virtual void runNodeScheduler ( eventList_t evMsgList, int nbrOfEvents)= 0;

  virtual void completeNodeScheduler( eventList_t evList ) = 0;

  virtual void terminateNodeScheduler( ) = 0;

  // GG XX should be in OCL functions
  //void configOCLVegasState( ocl_vegas_state_t *vegasState, size_t  nbrOfRequestedPoints );
  //void validateOCLVegasState( ocl_vegas_state_t *vegasState, const RunConfig *cfg );
 
  // Integration selection
  // EvalOnGrid
  ComputationType_t getComputationType() { return computationType;}
  void setComputationType(ComputationType_t ct) { computationType=ct;}

};

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef	__cplusplus
}
#endif

#endif	/* NODESCHEDULER_H */

