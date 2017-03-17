/* 
 * File:   RunConfig.h
 * Author: grasseau
 *
 * Created on 19 mai 2015, 11:46
 */

#ifndef MEMDataFormats_MEMDataFormats_RUNCONFIG_H
#define	MEMDataFormats_MEMDataFormats_RUNCONFIG_H

# include <stdint.h>
# include <string>
# include <vector>
# include <climits>
# include <iostream> 

using namespace std;

typedef enum {
  BadEventType  = -1,
  MCEventType   =  0,
  Run1EventType =  1,
  Run1_etau_EventType = 2,
  PyRun1EventType =  3,
  PyRun2EventType = 4,
} EventType_t;

class RunConfig {
    public :
        ~RunConfig();
//        explicit RunConfig(const char *configname);
        explicit RunConfig(); // const std::string configname
        void print() const;
  
      // File name of the Python configuration
      std::string configName_;
      EventType_t eventType_;
        
      int64_t maxNbrOfEventsToRead_; // 
      // Inputs
      vector<string> inputFileNames_;
  
      // Vector 
      std::string integralsOutFileName_;
//      std::string fctValuesOutFileName_;
      string treeName_;
      string MGParamCardFileName_;
      string LHAPDFFileName_;  
//      string LUTOneProngPi0FileName_;  
//      string LUTThreeProngsFileName_;  
  
      // New class order
  
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      // Input ROOT Files
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
  
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      // Physics 
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
      double sqrtS_;
      double Q_;/**/
      double mTau_;
      double mHiggs_;
      double mZ_ ;
        
        // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // Vegas
        // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        int nbrOfDimttH_;
        int nbrOfDimttZ_;
        int nbrOfDimttW_;
        int nbrOfDimttjets_;
        int nbrOfDimttbar_SL_;
        int nbrOfDimttbar_DL_;
        int nbrOfDimttZ_Zll_;
        int nbrOfDimttH_miss_;
        int nbrOfDimttZ_miss_;
        int nbrOfDimttjets_miss_;
        int nbrOfDimttbar_SL_miss_;
        int nbrOfDimttZ_Zll_miss_;
        //
        // Integration Boundaries
        vector<double> phiNu_tlep_Boundaries_;
        vector<double> cosThetaNu_tlep_Boundaries_;
        vector<double> phiNu_ttau_ttbar_DL_ttW_Boundaries_;
        vector<double> cosThetaNu_ttau_ttbar_DL_ttW_Boundaries_;
        vector<double> phiNu_W_ttW_Boundaries_;
        vector<double> cosThetaNu_W_ttW_Boundaries_;

        vector<double> phi_missing_jet_Boundaries_ ;
        vector<double> cosTheta_missing_jet_Boundaries_ ;
 
        //
        // Random generator
        bool flagSameRNG_;
  
        // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        // Function to integrate (integrand))
        // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        //
        // Fct evaluation
        bool flagTFLepTau_; 
        bool flagTFHadTau_;
        bool flagTFFake_;
        bool flagTFMET_   ;
        bool flagTFJet1_  ;
        bool flagTFJet2_  ;
        bool flagTFBJet_leptop_  ;
        bool flagTFBJet_hadtop_  ;
        bool flagTF_fakelep_  ;
        bool flagTF_fakeleptau_  ;
        bool flagTFTop_  ;
        bool flagJac_    ;
        bool flagWME_     ;
        
        bool runttZ_integration_  ;
        bool runttW_integration_  ;
        bool runttjets_integration_  ;
        bool runttbar_SL_integration_  ;
        bool runttbar_DL_integration_  ;
        bool runttbar_DL_fakelep_integration_  ;
        bool runttZ_Zll_integration_  ;

        bool run_missing_jet_integration_  ;
        bool force_missing_jet_integration_  ;
        bool force_missing_jet_integration_ifnoperm_  ;

        int nbrOfPoints_ttZ_  ;
        int nbrOfPoints_ttH_  ;
        int nbrOfPoints_ttW_  ;
        int nbrOfPoints_ttjets_  ;
        int nbrOfPoints_ttbar_SL_  ;
        int nbrOfPoints_ttbar_DL_  ;
        int nbrOfPoints_ttZ_Zll_  ;

        int nbrOfPoints_ttZ_miss_  ;
        int nbrOfPoints_ttH_miss_  ;
        int nbrOfPoints_ttjets_miss_  ;
        int nbrOfPoints_ttbar_SL_miss_  ;
        int nbrOfPoints_ttZ_Zll_miss_  ;
        
        int nbrOfPermut_per_jet_ ;

        double CI_TFJet_ ;
        bool use_pT_TFJet_ ;
        bool use_top_compatibility_check_ ;

        bool include_hadrecoil_ ;
        bool force_nonzero_integral_ ;

        double eta_acceptance_ ;
        double jet_radius_ ;
        double dR_veto_jet_lep_ ;
        double rel_iso_lep_ ;
        double pT_cut_ ;

        //
        // Verbose mode
        int verbose_;
        //
        // Madgraph computation ()
        int MEversion_;/**/
        
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      // Others 
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
      int nbOfTriplets_ = 0; // nb for all Events
      int nbOfNplets_ = 0; // nb for all Events
      int ttlnbOfElectrons_ = 0; // nb for all Events
      int ttlnbOfMuons_ = 0; // nb for all Events
      int ttlnbOfTaus_ = 0; // nb by Event
      int nbOfElectrons_ = 0; // nb by Event
      int nbOfMuons_ = 0; // nb by Event
      int nbOfTaus_ = 0; // nb by Event
      
};

class ConfigException: public exception {
  virtual const char* what() const throw()
  {
    return "My exception happened here";
  }
};

#endif	/* RUNCONFIG_H */

