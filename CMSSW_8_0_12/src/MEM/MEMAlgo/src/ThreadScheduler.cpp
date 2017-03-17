# include <fstream>
# include <iostream>
# include <iomanip>
# include <string>
# include <sstream>

# include "MEM/MEMAlgo/interface/ThreadScheduler.h"

# include "MEM/MEMAlgo/interface/MGIntegration.h"
# include "MEM/MEMAlgo/interface/LHAPDF.h"
# include "MEM/MEMAlgo/interface/Processes.h"



/* GG XXX EvalOnGrid, future use
void
dumpDBufferAsEvalOnGrid ( ocl_vegas_state_t *vegasState) {
  double X[MaxNumberOfDIMs], dX[MaxNumberOfDIMs];
  double index[MaxNumberOfDIMs];
  int d, k;
  
  int nbrOfDims = vegasState->dim;
  
  long int nbrPoints = pow( (double) nbrPointsPerAxe, (double) nbrOfDims ) + 0.5;
  for (d=0; d < nbrOfDims; d++) {
    dX[d] =  ( vegasState->xu[d] - vegasState->xl[d]  ) 
           / nbrPointsPerAxe;
  } 
  
  printf("  Eval On Grid (ev = %6ld) \n", vegasState->eventID);
  printf("     Idx  mTau2    PTauLep   cosTLTH  PQuark1  PQuark2    wME \n" ); 
  for (k=0; k < nbrPoints; k++) {
    int remain = k;
    for( d=0; d<nbrOfDims; d++ ) {
      index[d] = remain % nbrPointsPerAxe;
      remain   = remain / nbrPointsPerAxe;
      X[d] =   vegasState->xl[d] + index[d]*dX[d];
    }   
    printf("    %4d %8e %8e %8e %8e %8e %8e \n", k, 
            X[0], X[1], X[2], X[3], X[4],
            vegasState->dBuffer[k] );
  }
}
*/
/* GG XXX Fonction mère : à changer avec les IntegrationMsg 
void
displayEvent ( QueueStatus_t *qs, ocl_vegas_state_t *vegasState, int nbrOfCPU, const ocl_vegas_state_t *vegasCpuState) {

  if ( qs->startEventIndex < 0) {
     printf("  Event  Integral  Error    ChiSq   \n");
     printf("====================================\n");
  } else {
    printf(" res %6ld %3d %3d %8e %8e %8f %10ld %2d %2d\n", vegasState->eventID, qs->startEventIndex, qs->nbrOfEvents, 
           vegasState->CumulatedIntegral, vegasState->CumulatedSigma,
           vegasState->ChiSquare, qs->tickEnd - qs->tickStart,
           qs->platform, qs->device);
    // Not used
    //
    printf(" liForTest %ld %ld %ld %ld\n", 
	   vegasState->liForTest[0],
	   vegasState->liForTest[1],
	   vegasState->liForTest[2],
	   vegasState->liForTest[3]
           );  
    printf(" dForTest %e %e %e %e \n", 
	   vegasState->dForTest[0],
	   vegasState->dForTest[1],
	   vegasState->dForTest[2],
	   vegasState->dForTest[3]
           );
    */
    
    // Debug    
    // printState(vegasState);
    // printDistSampling(vegasState);
/*
    int cpu;
    for (cpu = 0; cpu  < nbrOfCPU; cpu++) {
      printf("     cpu = %3d %.6f %.6f \n", cpu, 
         vegasCpuState[cpu].CumulatedIntegral, vegasCpuState[cpu].CumulatedSigma);
    }
  }
}
*/


ThreadScheduler::ThreadScheduler():NodeScheduler() {
  
  totalNbrOfQueueAndProcess=0;
  remainingTasks=0;
 
  integration = 0;
}

ThreadScheduler::~ThreadScheduler() {
  // GG XXX GSL deallocation ?
   // Free space
   // Free space
  // Must be done in MPIScheduler logStream_.close();  
  delete integration;
}

/*
 * @brief Init the Device Scheduler
 *        Done once, fused qcpu queus are created 
 *        (because one device by virtual core on CPU) 
 * @param cfg      Configuration
 * @param mpi_rank MPI rank to init OCL
 *  
 */
void ThreadScheduler::initNodeScheduler( RunConfig cfg, int mpi_rank ) {
 
  NodeScheduler::initNodeScheduler( cfg, mpi_rank );

  // Configure integration
//  integration = new MGIntegration( cfg->configName_.c_str(), cfg, mpi_rank ); 
//  integration = new MGIntegration( cfg->configFileName.c_str(), mpi_rank ); 
  integration = new MGIntegration( cfg, mpi_rank ); 
  
  // Log file
  int nb_digits = log10( (double)(mpi_rank) +1 ) + 1; 
  nb_digits = max( nb_digits, 1);
  int str_len = nb_digits + 8 + cfg.configName_.size();
  //
  char *logNameFile = new char[str_len];
  sprintf( logNameFile, "%s.%d.log", cfg.configName_.c_str(), mpi_rank);
  // cout << "Debug: size = " << str_len << "file_name=[" << logNameFile << "]"<< endl;
  logStream_.open( logNameFile );
  delete logNameFile;
  //
  integration->setLogStream( logStream_ );  

  // Bind function to integrate to ROOT(GSL)VEGAS MC algorithm
//  int nbrOfDim = integration->getNbrOfDim();
  
  //
  // Init Vegas with the greatest number of evaluation points
  //
  //printf("initNodeScheduler: number of boxes %ld\n", hostState.nbrBoxesInDomain);

  // GSL init
  rng = gslRNGInit( );
  saveRNGSeed = gsl_rng_clone( rng );
/*  gslFctDescr = { 
    &wrapper_evalHTauTauLepHadCollinearCosTheta, integration->nbrOfDim_, integration };
  state = gsl_monte_vegas_alloc ( integration->nbrOfDim_ );  */
  fctDescr_ttH = { 
    &wrapper_evalttH, (size_t) integration->nbrOfDim_ttH_, integration };
  state_ttH = gsl_monte_vegas_alloc ( integration->nbrOfDim_ttH_ );  
    
  // LHAPDF
  // Fortran init 
//   const int SUBSET = 1;
   const string NAME = cfg.LHAPDFFileName_.c_str();
//  const string NAME = "cteq66";
   LHAPDF::initPDFSet( PDFSet, cfg.LHAPDFFileName_.c_str());
//  LHAPDF::initPDFSet(NAME, LHAPDF::LHGRID, SUBSET);
  // pdf_t *lhapdf_ = PDF_new();
  // PDF_fillData( lhapdf_);
//    const double Q = 10.0, mz = 91.2;
  /*cout << "alphas(mz) = " << LHAPDF::alphasPDF(mz) << endl;
  cout << "qcdlam4    = " << LHAPDF::getLam4(SUBSET) << endl;
  cout << "qcdlam5    = " << LHAPDF::getLam5(SUBSET) << endl;
  cout << "orderPDF   = " << LHAPDF::getOrderPDF() << endl;
  cout << "xmin       = " << LHAPDF::getXmin(SUBSET) << endl;
  cout << "xmax       = " << LHAPDF::getXmax(SUBSET) << endl;
  cout << "q2min      = " << LHAPDF::getQ2min(SUBSET) << endl;
  cout << "q2max      = " << LHAPDF::getQ2max(SUBSET) << endl;
  cout << "orderalfas = " << LHAPDF::getOrderAlphaS() << endl;
  cout << "num flav   = " << LHAPDF::getNf() << endl;
  cout << "name       = " << NAME << endl;
  cout << "number     = " << LHAPDF::numberPDF() << endl;
  cout << endl;*/

//  const int NUMBER = LHAPDF::numberPDF();

//  const double x = (1. - 0.5) / 10. ;
//  std::cout << "fg0 : " << LHAPDF::xfx( PDFSet, x, Q,  0 ) << std::endl;

/*  for (int n = 0; n < NUMBER + 1; ++n) {
    cout << "Set number: " << n << endl;
    LHAPDF::initPDF(n);
    for (int ix = 1; ix < 11; ++ix) {
      const double x = (ix - 0.5) / 10.0;
      cout << "x=" << x << ", Q=" << Q << ", f=0: " << LHAPDF::xfx(x, Q, 1) << endl;
    }
    cout << endl;
  }*/
  LHAPDF::getDescription();
  // MadGraph
//  initHTauTauMEProcesses( cfg ); 
  //init_ttH_HTauTauMEProcesses( cfg->configName_.c_str() ); 
  init_ttH_HTauTauMEProcesses( cfg ); 

}

/*
 * @brief Clean  queueStatus to restart the ThreadScheduler
 *        (dim, nbrOfCalls, limits, ...)
 * @param ...
 */
/*
void ThreadScheduler::cleanNodeScheduler( ) {
  int qcpu = 0; 
  int qcpu_max = totalNbrOfQueueAndProcess;
  for ( qcpu = 0; qcpu < qcpu_max; qcpu++) {
    queueStatus[qcpu].status  = QueueStatus_t::ComputationNotActive;
    queueStatus[qcpu].startEventIndex = -1;
  }
}*/

void ThreadScheduler::runNodeScheduler (
                eventList_t evList,
                int nbrOfEvents) {
    
//  long int startedTasks = 0;
  long int  endedTasks=0;
//  int nbrSubmittedEvents = 0;
  
  int k; 
  cerr << "  compute "  << nbrOfEvents << endl;  
  std::cout << "runNodeScheduler :" << integration->force_missing_jet_integration_ << std::endl;

  for ( k=0; k < nbrOfEvents; k++) {
    std::cout << k << "/" << nbrOfEvents << std::endl;
    integration->setEventParameters( evList[k], integration->force_missing_jet_integration_ );
    integration->copyBoundaries( &evList[k]);
    std::cout<<"nbrOfPermut_ runNodeScheduler "<<evList[k].nbrOfPermut_<<std::endl;

/*  
# if 0  
     // Integrate DY
    integration->signalME_ = false;

    // Set boundaries
    integration->lowerValues_[mTauTau2_id] = integration->mTauTau2DY_[0];
    integration->upperValues_[mTauTau2_id] = integration->mTauTau2DY_[1];    
    integration->lowerValues_[PTauLep_id] = integration->PTauLepBoundariesDY_[0];
    integration->upperValues_[PTauLep_id] = integration->PTauLepBoundariesDY_[1];
    integration->lowerValues_[cosThetaTauLepTauHad_id] = integration->cosTauLepTauHadBoundariesDY_[0];
    integration->upperValues_[cosThetaTauLepTauHad_id] = integration->cosTauLepTauHadBoundariesDY_[1];    
   
    // Copy the initial state to the current RNG State
    if ( integration->flagSameRNG_ )
        gsl_rng_memcpy ( rng, saveRNGSeed);

    gslIntegrate( &fctDescr, integration, rng, state, true,
            &integralParameters.integralDY_, 
            &integralParameters.stderrDY_,
            &integralParameters.chiSquareDY_);
# endif    
*/    
    // Integrate VBF
/*
    integration->signalME_ = true;
    integration->lowerValues_[mTauTau2_id] = integration->mTauTau2BoundariesVBF_[0];
    integration->upperValues_[mTauTau2_id] = integration->mTauTau2BoundariesVBF_[1];
    integration->lowerValues_[PTauLep_id]  = integration->PTauLepBoundariesVBF_[0];
    integration->upperValues_[PTauLep_id]  = integration->PTauLepBoundariesVBF_[1];
    integration->lowerValues_[cosThetaTauLepTauHad_id] = integration->cosTauLepTauHadBoundariesVBF_[0];
    integration->upperValues_[cosThetaTauLepTauHad_id] = integration->cosTauLepTauHadBoundariesVBF_[1];
*/
    std::cout <<"nbrOfPermut_ 2 : " << integration->nbrOfPermut_<<std::endl;
	for(int perm=0; perm<integration->nbrOfPermut_; perm++){
	      
        integration->initVersors(perm);          
 	integration->signalME_ = true;
	integration->m_TauTau_2_ = pow(integration->mTauTau_ttH_[perm],2);
	integration->lowerValues_[PTauLep_id] = integration->PTauLep_ttH_Lower_[perm];
	integration->upperValues_[PTauLep_id] = integration->PTauLep_ttH_Upper_[perm];
	integration->lowerValues_[cosThetaTauLepTauHad_id] = integration->cosTheta_diTau_ttH_Lower_[perm];
	integration->upperValues_[cosThetaTauLepTauHad_id] = integration->cosTheta_diTau_ttH_Upper_[perm];
	integration->lowerValues_[EQuark1_id] = integration->EQuark1_Lower_[perm];
	integration->upperValues_[EQuark1_id] = integration->EQuark1_Upper_[perm];
	integration->lowerValues_[cosThetaNu_tlep_id] = integration->cosThetaNu_tlep_Boundaries_[0];
	integration->upperValues_[cosThetaNu_tlep_id] = integration->cosThetaNu_tlep_Boundaries_[1];
	integration->lowerValues_[phiNu_tlep_id] = integration->phiNu_tlep_Boundaries_[0];
	integration->upperValues_[phiNu_tlep_id] = integration->phiNu_tlep_Boundaries_[1];
		
	integration->integr_EfficiencyttH_ = 0;
	integration->tot_DrawsttH_ = 0;
		
   // Compute Integral          
   // Copy the initial state to the current RNG State
    if ( integration->flagSameRNG_ )
      gsl_rng_memcpy ( rng, saveRNGSeed);

	gslIntegrate( &fctDescr_ttH, integration,  integration->nbrOfDim_ttH_, integration->nbrOfPoints_ttH_, rng, state_ttH, true,
			      &evList[k].integralttH_[perm], 
			      &evList[k].stderrttH_[perm],
			      &evList[k].chiSquarettH_[perm]);     
/*    gslIntegrate( &gslFctDescr, integration, rng, state, true,
                &evList[k].integralVBF_,
                &evList[k].stderrVBF_,   
                &evList[k].chiSquareVBF_ ); */   

                /// GG XXX queueStatus[qcpu].tickEnd = getticks();
/*    cerr << "   VBF Integral = "   << evList[k].integralVBF_ 
         << ", stderr = "     << evList[k].stderrVBF_
         << ", Chi-square = " << evList[k].chiSquareVBF_
         << endl;
    if (integration->verbose_ >= ResultLevel) {
      *integration->stream_ << "  VBF Integral = "   << evList[k].integralVBF_ 
         << ", stderr = "     << evList[k].stderrVBF_
         << ", Chi-square = " << evList[k].chiSquareVBF_
         << endl;
      *integration->stream_ << "  DY Integral = "   << evList[k].integralDY_ 
         << ", stderr = "     << evList[k].stderrDY_
         << ", Chi-square = " << evList[k].chiSquareDY_
         << endl;
    }*/
    
    } // end perm loop 
    
    endedTasks++;
  }

}


/*
 * @brief Wait all device taks finished
 * @param ...
 */  
void ThreadScheduler::completeNodeScheduler( eventList_t evList ) {

}


void ThreadScheduler::terminateNodeScheduler( ) {
  
}
