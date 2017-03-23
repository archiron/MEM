# include <fstream>
# include <iostream>
# include <iomanip>
# include <string>
# include <sstream>

# include "MEM/MEMAlgo/interface/ThreadScheduler.h"

# include "MEM/MEMAlgo/interface/MGIntegration.h"
# include "MEM/MEMAlgo/interface/LHAPDF.h"
# include "MEM/MEMAlgo/interface/Processes.h"

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
  fctDescr_ttH_miss = { 
    &wrapper_evalttH, (size_t) integration->nbrOfDim_ttH_miss_, integration };
  state_ttH_miss = gsl_monte_vegas_alloc ( integration->nbrOfDim_ttH_miss_ ); 
//  std::cout << "retour state_ttH : " << state_ttH << std::endl; // TEMPORAIRE
    
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
  
  //std::cout << "ThreadScheduler::initNodeScheduler : appel init_ttH_HTauTauMEProcesses" << std::endl; // TEMPORAIRE
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
  cerr << "  compute " << nbrOfEvents << " events" << endl;  
  //std::cout << "runNodeScheduler force_missing_jet_integration_ : " << integration->force_missing_jet_integration_ << std::endl;

  for ( k=0; k < nbrOfEvents; k++) {
    //std::cout << k+1 << "/" << nbrOfEvents << std::endl;
    
    integration->setEventParameters( evList[0], integration->force_missing_jet_integration_ ); // seems to be OK (get all datas)
    
    std::cout << "nbrOfPermut_ runNodeScheduler1 " << evList[0].nbrOfPermut_ << std::endl; // TEMPORAIRE
    for( int perm = 0; perm < evList[0].nbrOfPermut_; perm++ ){  // TEMPORAIRE
        std::cout << "integralttH_[" << perm << "] : " << evList[0].integralttH_[perm] ; // TEMPORAIRE
        std::cout << " - include_perm_ttH_[" << perm << "] : " << evList[0].include_perm_ttH_[perm] << std::endl; // TEMPORAIRE
    }/**/ // TEMPORAIRE
    
    integration->copyBoundaries( &evList[0]);                                                  // seems to be OK (get all datas)
    
    std::cout << "nbrOfPermut_ runNodeScheduler2 " << evList[k].nbrOfPermut_ << std::endl; // TEMPORAIRE
    for( int perm = 0; perm < evList[0].nbrOfPermut_; perm++ ){  // TEMPORAIRE
        std::cout << "integralttH_[" << perm << "] : " << evList[0].integralttH_[perm] ; // TEMPORAIRE
        std::cout << " - include_perm_ttH_[" << perm << "] : " << evList[0].include_perm_ttH_[perm] << std::endl; // TEMPORAIRE
    }/**/ // TEMPORAIRE

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
    std::cout << "nbrOfPermut_ ThreadScheduler : " << integration->nbrOfPermut_ << std::endl;

    if(integration->integration_type_ == 0 && !integration->force_missing_jet_integration_){
        
        for(int perm=0; perm<integration->nbrOfPermut_; perm++){
            std::cout << "permut : " << perm + 1 << "/" << integration->nbrOfPermut_ << std::endl;
            integration->initVersors(perm);          

            if( integration->include_perm_ttH_[perm] ){
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

                std::cout << "AVANT GSL" << std::endl;
                //std::cout << "integration->nbrOfDim_ttH_ : "    << integration->nbrOfDim_ttH_    << std::endl;
                //std::cout << "integration->nbrOfPoints_ttH_ : " << integration->nbrOfPoints_ttH_ << std::endl;
                //std::cout << "state_ttH : " << state_ttH << std::endl; // TEMPORAIRE
                std::cout << "&evList["<< k << "].integralttH_["  << perm <<"] : " << evList[k].integralttH_[perm]  << std::endl;
                //std::cout << "&evList["<< k << "].stderrttH_["    << perm <<"] : " << evList[k].stderrttH_[perm]    << std::endl;
                //std::cout << "&evList["<< k << "].chiSquarettH_[" << perm <<"] : " << evList[k].chiSquarettH_[perm] << std::endl;
        
                gslIntegrate( &fctDescr_ttH, integration,  integration->nbrOfDim_ttH_, integration->nbrOfPoints_ttH_, rng, state_ttH, true,
                    &evList[k].integralttH_[perm], 
                    &evList[k].stderrttH_[perm],
                    &evList[k].chiSquarettH_[perm]);     

                std::cout << "APRES GSL" << std::endl;
                //std::cout << "integration->nbrOfDim_ttH_ : "    << integration->nbrOfDim_ttH_    << std::endl;
                //std::cout << "integration->nbrOfPoints_ttH_ : " << integration->nbrOfPoints_ttH_ << std::endl;
                //std::cout << "state_ttH : " << state_ttH << std::endl;
                std::cout << "&evList["<< k << "].integralttH_["  << perm <<"] : " << evList[k].integralttH_[perm]  << std::endl;
                //std::cout << "&evList["<< k << "].stderrttH_["    << perm <<"] : " << evList[k].stderrttH_[perm]    << std::endl;
                //std::cout << "&evList["<< k << "].chiSquarettH_[" << perm <<"] : " << evList[k].chiSquarettH_[perm] << std::endl;

                /// GG XXX queueStatus[qcpu].tickEnd = getticks();
                //evList[k].compTimettH_[perm] = ((float)t)/CLOCKS_PER_SEC;
                evList[k].totalDrawsttH_[perm] = integration->tot_DrawsttH_;
                evList[k].integrationEfficiencyttH_[perm] = integration->integr_EfficiencyttH_;
                //evList[k].weight_ttH_ += evList[k].integralttH_[perm];
                //var_ttH += pow(evList[k].stderrttH_[perm],2);
            }
        } // end perm loop 
    } // end if integration_type_ ..
    else if((integration->integration_type_ == 1 && integration->run_missing_jet_integration_) || (integration->integration_type_ == 0 && integration->force_missing_jet_integration_)){
      
        //Loop over permutations
        for(int perm=0; perm<integration->nbrOfPermut_; perm++){
            std::cout << "permut miss : " << perm + 1 << "/" << integration->nbrOfPermut_ << std::endl;
            integration->initVersors_miss(perm);          
	
        //
        // Integrate ttH
        //

            if( integration->include_perm_ttH_[perm] ){
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
                integration->lowerValues_[cosTheta_missing_jet_id] = integration->cosTheta_missing_jet_Boundaries_[0];
                integration->upperValues_[cosTheta_missing_jet_id] = integration->cosTheta_missing_jet_Boundaries_[1];
                integration->lowerValues_[phi_missing_jet_id] = integration->phi_missing_jet_Boundaries_[0];
                integration->upperValues_[phi_missing_jet_id] = integration->phi_missing_jet_Boundaries_[1];	 	 
	  
                integration->tot_DrawsttH_ = 0;
                integration->integr_EfficiencyttH_ = 0;
	  
                // Compute Integral          
                if ( integration->flagSameRNG_ )
                    // Copy the initial state to the current RNG State
                    gsl_rng_memcpy ( rng, saveRNGSeed );
	  
                gslIntegrate( &fctDescr_ttH_miss, integration, integration->nbrOfDim_ttH_miss_, integration->nbrOfPoints_ttH_miss_, rng, state_ttH_miss, true,
                    &evList[k].integralttH_[perm], 
                    &evList[k].stderrttH_[perm],
                    &evList[k].chiSquarettH_[perm]);     
	  
                //evList[k].compTimettH_[perm] =  ((float)t)/CLOCKS_PER_SEC;
                evList[k].totalDrawsttH_[perm] = integration->tot_DrawsttH_;
                evList[k].integrationEfficiencyttH_[perm] = integration->integr_EfficiencyttH_;
                //evList[k].weight_ttH_ += evList[k].integralttH_[perm]; 
                //var_ttH += pow(evList[k].stderrttH_[perm],2);
            }
        } // end perm loop 
    }
    endedTasks++;
  } // end nbrOfEvents loop 

}


/*
 * @brief Wait all device taks finished
 * @param ...
 */  
void ThreadScheduler::completeNodeScheduler( eventList_t evList ) {

}


void ThreadScheduler::terminateNodeScheduler( ) {
  
}
