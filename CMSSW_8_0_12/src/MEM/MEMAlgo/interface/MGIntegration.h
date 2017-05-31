/* 
 * File:   GMIntegrator.h
 * Author: grasseau
 *
 * Created on 21 novembre 2014, 15:55
 */

// GG TODO : put the message structure inside the class to avoid 2 separate data types
// AC : most of those variables are integrated into runConfig

#ifndef MGINTEGRATOR_H
#define	MGINTEGRATOR_H

# include <stdint.h>
# include <vector>

# include "TVector3.h"
# include "TLorentzVector.h"

#include <stdlib.h>
#include "gsl/gsl_math.h"
#include "gsl/gsl_monte_vegas.h"
#include "gsl/gsl_pow_int.h"

// # include "IOLib/FctValuesRootOutputs.h"
// # include "IOLib/FctValuesRootOutputs_ttW.h"
// # include "IOLib/FctValuesRootOutputs_ttjets.h"
// # include "IOLib/FctValuesRootOutputs_ttbar_SL.h"
// # include "IOLib/FctValuesRootOutputs_ttbar_DL.h"
# include "MEMDataFormats/MEMDataFormats/interface/RunConfig.h"

# define PI TMath::Pi()

////////////////////////////////////////////////////////////////////////////////
// Notations / Abreviations used
// -----------------------------
//   + 'Lep' Lepton or leptonic
//   + 'Had' Hadron or hadronic
//   + 'ATau' anti-Tau
//   + 't' top
// Prefixes:
//   + 'ev' Event observable (hat variable in Thomas et al document)
//     example: evLep4P_
//   + 'P_' Module of the momentum
//     example: P_Tau
//   + 'm_' mass
//     example: m_TauTau_2
// Suffixes:
//   + '_' Class private attribute 
//     example: nbrOfDim_
//   + '4P' Momentum quadri-vector 
//     example: evJet1_4P_
//
////////////////////////////////////////////////////////////////////////////////

using namespace std;


// Argument index of the VEGAS integration 
# define PTauLep_id        0 //module of the momentum of the tau
# define cosThetaTauLepTauHad_id 1 //cos of the TauLep/TauHad angle
# define EQuark1_id     2 //module of the momentum of the light Quark1 from hadronic top
# define cosThetaNu_tlep_id     3 //costheta of the neutrino from the leptonic top
# define phiNu_tlep_id          4 //phi of the neutrino from the leptonic top
# define cosTheta_missing_jet_id     5 //costheta of the missing jet
# define phi_missing_jet_id          6 //phi of the missing jet


//ttjets
# define EQuark1_ttjets_id     0 //module of the momentum of the light Quark1 from hadronic top
# define cosThetaNu_tlep_ttjets_id     1 //costheta of the neutrino from the leptonic top
# define phiNu_tlep_ttjets_id          2 //phi of the neutrino from the leptonic top
#define POutFake_ttjets_id        3  //module of the momentum of the final state parton which fakes a tauh
# define cosTheta_missing_jet_ttjets_id     4 //costheta of the missing jet
# define phi_missing_jet_ttjets_id          5 //phi of the missing jet

//ttbar_SL
# define EQuark1_ttbar_SL_ttZ_Zll_id     0 //module of the momentum of the light Quark1 from hadronic top
# define cosThetaNu_tlep_ttbar_SL_ttZ_Zll_id     1 //costheta of the neutrino from the leptonic top
# define phiNu_tlep_ttbar_SL_ttZ_Zll_id          2 //phi of the neutrino from the leptonic top
# define cosTheta_missing_jet_ttbar_SL_ttZ_Zll_id     3 //costheta of the missing jet
# define phi_missing_jet_ttbar_SL_ttZ_Zll_id          4 //phi of the missing jet

//ttbar_DL
# define cosThetaNu_tlep_ttbar_DL_ttW_id     0 //costheta of the neutrino from the leptonic top
# define phiNu_tlep_ttbar_DL_ttW_id          1 //phi of the neutrino from the leptonic top
# define cosThetaNu_ttau_ttbar_DL_ttW_id     2 //costheta of the neutrino from the tauonic top
# define phiNu_ttau_ttbar_DL_ttW_id          3 //phi of the neutrino from the tauonic to
# define PTauHad_ttbar_DL_ttW_id                     4 //module of the momentum of the tau


//ttW
#define cosThetaNu_W_ttW_id      5  //cosTheta of the nu from W
#define phiNu_W_ttW_id           6  //phi of the nu from W


// Verbose values
# define ResultLevel      1  // Display integration results
# define IntegrationLevel 2  // Display integration parameters
# define IntegrandLevel   3  // Display Integrand information

// Dimension max of the integration
# define DimensionMax 8
# define nbrOfPermutMax 40
# define nbrOfJetMax 10

# define integration_type_wo_miss 0
# define integration_type_w_miss  1


typedef struct MGIntegration_s {
  
  // Computation mode
  int            verbose_;

  int event_type_; //0 = loose event (1 MVA + 1 fakeable lepton), 1 = tight event (2 MVA leptons)
  int integration_type_ ; //0 = no missing jet, 1 = missing jet
 
  // MC integration
  int    nbrOfDim_ttH_;
  int    nbrOfDim_ttZ_;
  int    nbrOfDim_ttW_;
  int    nbrOfDim_ttjets_;
  int    nbrOfDim_ttbar_SL_;
  int    nbrOfDim_ttbar_DL_;
  int    nbrOfDim_ttZ_Zll_;


  int    nbrOfDim_ttH_miss_;
  int    nbrOfDim_ttZ_miss_;
  int    nbrOfDim_ttjets_miss_;
  int    nbrOfDim_ttbar_SL_miss_;
  int    nbrOfDim_ttZ_Zll_miss_;
  
  int64_t   nbrOfPoints_ttH_;
  int64_t   nbrOfPoints_ttZ_;
  int64_t   nbrOfPoints_ttW_;
  int64_t   nbrOfPoints_ttjets_;
  int64_t   nbrOfPoints_ttbar_SL_;
  int64_t   nbrOfPoints_ttbar_DL_;
  int64_t   nbrOfPoints_ttZ_Zll_;

  int64_t   nbrOfPoints_ttH_miss_;
  int64_t   nbrOfPoints_ttZ_miss_;
  int64_t   nbrOfPoints_ttjets_miss_;
  int64_t   nbrOfPoints_ttbar_SL_miss_;
  int64_t   nbrOfPoints_ttZ_Zll_miss_;

  int    nbrOfPermut_per_jet_;
  int    nbrOfPermut_;

  double lowerValues_[DimensionMax];
  double upperValues_[DimensionMax];

  double CI_TFJet_;
  bool use_pT_TFJet_;
  bool use_top_compatibility_check_;

  //Common boundaries
  double phiNu_tlep_Boundaries_[2];
  double cosThetaNu_tlep_Boundaries_[2];
  double phi_missing_jet_Boundaries_[2];
  double cosTheta_missing_jet_Boundaries_[2];

  // Missing jet
  double eta_acceptance_;
  double jet_radius_;
  double dR_veto_jet_lep_;
  double rel_iso_lep_;
  double pT_cut_;

  double EQuark1_Lower_[nbrOfPermutMax];
  double EQuark1_Upper_[nbrOfPermutMax];

  // ttH
  int include_perm_ttH_[nbrOfPermutMax];
  double mTauTau_ttH_[nbrOfPermutMax];
  double PTauLep_ttH_Lower_[nbrOfPermutMax]; 
  double PTauLep_ttH_Upper_[nbrOfPermutMax];  
  double cosTheta_diTau_ttH_Lower_[nbrOfPermutMax];
  double cosTheta_diTau_ttH_Upper_[nbrOfPermutMax];  

  // ttZ
  int include_perm_ttZ_[nbrOfPermutMax];
  double mTauTau_ttZ_[nbrOfPermutMax];
  double PTauLep_ttZ_Lower_[nbrOfPermutMax]; 
  double PTauLep_ttZ_Upper_[nbrOfPermutMax];  
  double cosTheta_diTau_ttZ_Lower_[nbrOfPermutMax];
  double cosTheta_diTau_ttZ_Upper_[nbrOfPermutMax];

  // ttW
  int include_perm_ttW_[nbrOfPermutMax];
  double phiNu_W_ttW_Boundaries_[2];
  double cosThetaNu_W_ttW_Boundaries_[2];  

  // ttjets
  double POutFakeBoundaries_ttjets_[2];
  int include_perm_ttjets_[nbrOfPermutMax];

  // ttbar_SL
  int include_perm_ttbar_SL_[nbrOfPermutMax];

  // ttbar_DL
  int include_perm_ttbar_DL_[nbrOfPermutMax];
  int include_perm_ttbar_DL_fakelep_tlep_[nbrOfPermutMax];
  int include_perm_ttbar_DL_fakelep_ttau_[nbrOfPermutMax];

  double phiNu_ttau_ttbar_DL_ttW_Boundaries_[2];
  double cosThetaNu_ttau_ttbar_DL_ttW_Boundaries_[2];
  double PTauHad_ttbar_DL_ttW_Boundaries_[2];

  // ttZ_Zll
  int include_perm_ttZ_Zll_[nbrOfPermutMax];

  // LHC 
  double sqrtS_; // LHC beam energy    
  
  // Particles
  double mTau_;
  double mHiggs_; 
  double mtop_; 
  double mW_;

  //
  //  EVENT
  //
  // Particles measure in the event
  // Event ID in the code
  uint64_t eventID_;
  uint64_t nRun_;
  uint64_t nLumi_;
  // Root Event ID
  int nEvent_;

  bool run_loose_event_integration_;

  double evLep1_4P_[4];
  double evLep2_4P_[4];
  double evBJet1_4P_[4];
  double evBJet2_4P_[4];
  //All untagged jets
  double evJets_4P_[nbrOfJetMax][4];
  int n_lightJets_;
  
  double evLep_Tau_4P_[4];
  double evLep_top_4P_[4];
  double evHadSys_Tau_4P_[4];
  double evJet1_4P_[4];
  double evJet2_4P_[4];
  double evBJet_leptop_4P_[4];
  double evBJet_hadtop_4P_[4];

  //
  // Tau Decay mode + lep type
  int HadtauDecayMode_;
  int lepton_Tau_Type_;
  int lepton_top_Type_;
  int lepton_Tau_charge;
  int lepton_top_charge;

  int lepton1_Type_;
  int lepton2_Type_;
  int lepton1_charge_;
  int lepton2_charge_;
 
  //
  // Reconstructed (Reco) MET
  double evRecoMET4P_[4];
  // MET covariance matrix
  // Order : (0,0), (0,1), (1,0), (1,1) 
  double evV_[4];
  
  // Computation / run
  bool flagTFLepTau_; 
  bool flagTFHadTau_;
  bool flagTFFake_;
  bool flagTFMET_;   
  bool flagTFJet1_;  
  bool flagTFJet2_;
  bool flagTFBJet_leptop_;  
  bool flagTFBJet_hadtop_;
  bool flagTF_fakelep_;  
  bool flagTF_fakeleptau_;  
  bool flagTFTop_;
  bool flagJac_; 
  bool flagWME_;

  // MET TF
  bool include_hadrecoil_; 
  bool force_nonzero_integral_;

  // Results  
  double integralttH_[nbrOfPermutMax]; // Chichi = {0.}
  double stderrttH_[nbrOfPermutMax];
  double chiSquarettH_[nbrOfPermutMax];
  int64_t  integrationEfficiencyttH_[nbrOfPermutMax];
  int64_t  totalDrawsttH_[nbrOfPermutMax];
  double compTimettH_[nbrOfPermutMax];
  //
  double integralttZ_[nbrOfPermutMax];
  double stderrttZ_[nbrOfPermutMax];
  double chiSquarettZ_[nbrOfPermutMax];
  int64_t  integrationEfficiencyttZ_[nbrOfPermutMax];
  int64_t  totalDrawsttZ_[nbrOfPermutMax];
  double compTimettZ_[nbrOfPermutMax];
  //
  double integralttW_[nbrOfPermutMax];
  double stderrttW_[nbrOfPermutMax];
  double chiSquarettW_[nbrOfPermutMax];
  int64_t  integrationEfficiencyttW_[nbrOfPermutMax];
  int64_t  totalDrawsttW_[nbrOfPermutMax];
  double compTimettW_[nbrOfPermutMax];
  //
  double integralttjets_[nbrOfPermutMax];
  double stderrttjets_[nbrOfPermutMax];
  double chiSquarettjets_[nbrOfPermutMax];
  int64_t  integrationEfficiencyttjets_[nbrOfPermutMax];
  int64_t  totalDrawsttjets_[nbrOfPermutMax];
  double compTimettjets_[nbrOfPermutMax];
  //
  double integralttbar_SL_[nbrOfPermutMax];
  double stderrttbar_SL_[nbrOfPermutMax];
  double chiSquarettbar_SL_[nbrOfPermutMax];
  int64_t  integrationEfficiencyttbar_SL_[nbrOfPermutMax];
  int64_t  totalDrawsttbar_SL_[nbrOfPermutMax];
  double compTimettbar_SL_[nbrOfPermutMax];
  //
  double integralttbar_DL_[nbrOfPermutMax];
  double stderrttbar_DL_[nbrOfPermutMax];
  double chiSquarettbar_DL_[nbrOfPermutMax];
  int64_t  integrationEfficiencyttbar_DL_[nbrOfPermutMax];
  int64_t  totalDrawsttbar_DL_[nbrOfPermutMax];
  double compTimettbar_DL_[nbrOfPermutMax];
  //
  double integralttbar_DL_fakelep_tlep_[nbrOfPermutMax];
  double stderrttbar_DL_fakelep_tlep_[nbrOfPermutMax];
  double chiSquarettbar_DL_fakelep_tlep_[nbrOfPermutMax];
  int64_t  integrationEfficiencyttbar_DL_fakelep_tlep_[nbrOfPermutMax];
  int64_t  totalDrawsttbar_DL_fakelep_tlep_[nbrOfPermutMax];
  double compTimettbar_DL_fakelep_tlep_[nbrOfPermutMax];
  //
  double integralttbar_DL_fakelep_ttau_[nbrOfPermutMax];
  double stderrttbar_DL_fakelep_ttau_[nbrOfPermutMax];
  double chiSquarettbar_DL_fakelep_ttau_[nbrOfPermutMax];
  int64_t  integrationEfficiencyttbar_DL_fakelep_ttau_[nbrOfPermutMax];
  int64_t  totalDrawsttbar_DL_fakelep_ttau_[nbrOfPermutMax];
  double compTimettbar_DL_fakelep_ttau_[nbrOfPermutMax];
  //
  double integralttZ_Zll_[nbrOfPermutMax];
  double stderrttZ_Zll_[nbrOfPermutMax];
  double chiSquarettZ_Zll_[nbrOfPermutMax];
  int64_t  integrationEfficiencyttZ_Zll_[nbrOfPermutMax];
  int64_t  totalDrawsttZ_Zll_[nbrOfPermutMax];
  double compTimettZ_Zll_[nbrOfPermutMax];

  
  //Efficiency filled during the integration
  int64_t integr_EfficiencyttH_;
  int64_t integr_EfficiencyttZ_;
  int64_t integr_EfficiencyttW_;
  int64_t integr_Efficiencyttjets_;
  int64_t integr_Efficiencyttbar_SL_;
  int64_t integr_Efficiencyttbar_DL_;
  int64_t integr_Efficiencyttbar_DL_fakelep_tlep_;
  int64_t integr_Efficiencyttbar_DL_fakelep_ttau_;
  int64_t integr_EfficiencyttZ_Zll_;

  int64_t  tot_DrawsttH_;
  int64_t  tot_DrawsttZ_;
  int64_t  tot_DrawsttW_;
  int64_t  tot_Drawsttjets_;
  int64_t  tot_Drawsttbar_SL_;
  int64_t  tot_Drawsttbar_DL_;
  int64_t  tot_Drawsttbar_DL_fakelep_tlep_;
  int64_t  tot_Drawsttbar_DL_fakelep_ttau_;
  int64_t  tot_DrawsttZ_Zll_;


  // MPI
  int MPIInfo_;
} IntegrationMsg_t;

class MGIntegration : public IntegrationMsg_t {
  
  
  // Integration class
  /*enum IntegrationType {
    // Hyp signal, background ???
    
    // Signal
    SignalLeptonLepton=0,
    SignalLeptonHadron, // Semi-leptonic
    SignalHadronHadron,
    Unknown
    } IntegrationType_t;*/
  
  // Jacobian term
  typedef enum JacobianKind {
    NoJacobian = 0,
    MassTauTauBarJacobian,
    UnknownJacobian
  } JacobianKind_t;
  
 // GG To fuse Msg struct and class
 // GG ??? 
 // private :
public :
  
  // Computation
  int            usePDF_;
  JacobianKind_t jacobianKind_;
  int            verbose_;
  
  //
  // Ouptuts
  //
  ofstream  *stream_;
  // Function values
  /*
  string            valuesOutFileName_; 
  ofstream         *valuesOut_;
  FctValuesOutRoot *valuesOutRoot_;
  //
  ofstream      *outgoingJetsOut_;
  //
  string            valuesOutFileName_ttW_; 
  ofstream         *valuesOut_ttW_;
  FctValuesOutRoot_ttW *valuesOutRoot_ttW_;
  //
  ofstream      *outgoingJetsOut_ttW_;
  //
  string            valuesOutFileName_ttjets_; 
  ofstream         *valuesOut_ttjets_;
  FctValuesOutRoot_ttjets *valuesOutRoot_ttjets_;
  //
  ofstream      *outgoingJetsOut_ttjets_;
  //
  string            valuesOutFileName_ttbar_SL_; 
  ofstream         *valuesOut_ttbar_SL_;
  FctValuesOutRoot_ttbar_SL *valuesOutRoot_ttbar_SL_;
  //
  ofstream      *outgoingJetsOut_ttbar_SL_;
  //
  string            valuesOutFileName_ttbar_DL_; 
  ofstream         *valuesOut_ttbar_DL_;
  FctValuesOutRoot_ttbar_DL *valuesOutRoot_ttbar_DL_;
  //
  ofstream      *outgoingJetsOut_ttbar_DL_;
*/

  // LHC 
  double sqrtS_; // LHC beam energy   
  // PDF  
  double Q_;      // PDF calibration (Gev)
  // Particles
  double mTau_;
  double mHiggs_; 
  double mtop_; 
  double mW_; 

  //
  //  EVENT
  //
  // Particles measure in the event
  int64_t eventID_;
  uint64_t nRun_;
  uint64_t nLumi_;
  uint64_t nEvent_;

  TLorentzVector evLep_Tau_4P_;
  TLorentzVector evLep_top_4P_;
  TLorentzVector evHadSys_Tau_4P_;
  TLorentzVector evJet1_4P_;
  TLorentzVector evJet2_4P_;
  TLorentzVector evBJet_leptop_4P_;
  TLorentzVector evBJet_hadtop_4P_;

  TLorentzVector evLep1_4P_;
  TLorentzVector evLep2_4P_;
  TLorentzVector evBJet1_4P_;
  TLorentzVector evBJet2_4P_;
  TLorentzVector evJets_4P_[nbrOfJetMax];
  int n_lightJets_;
  
  // Tau Decay mode + lep type
  int HadtauDecayMode_;
  int lepton_Tau_Type_;
  int lepton_top_Type_;
  int lepton_Tau_charge_;
  int lepton_top_charge_;
  
  int lepton1_Type_;
  int lepton2_Type_;
  int lepton1_charge_;
  int lepton2_charge_;

  
  // Reconstructed (Reco) MET
  TLorentzVector evRecoMET4P_;
  // MET covariance matrix
  double evV_[4];
  
  // Run ttZ/ttW integration
  bool runttZ_integration_;
  bool runttW_integration_;
  bool runttjets_integration_;
  bool runttbar_SL_integration_;
  bool runttbar_DL_integration_;
  bool runttbar_DL_fakelep_integration_;
  bool runttZ_Zll_integration_;

  bool run_missing_jet_integration_;
  bool force_missing_jet_integration_;
  bool force_missing_jet_integration_ifnoperm_;

  // RNG
  bool flagSameRNG_;
  
  // Computation / run
  bool flagTFLepTau_; 
  bool flagTFHadTau_;
  bool flagTFFake_;
  bool flagTFMET_;   
  bool flagTFJet1_;  
  bool flagTFJet2_;
  bool flagTFBJet_leptop_;  
  bool flagTFBJet_hadtop_;
  bool flagTF_fakelep_;   
  bool flagTF_fakeleptau_; 
  bool flagTFTop_;     
  bool flagJac_; 
  bool flagWME_; 

  // MET TF
  bool include_hadrecoil_;
  bool force_nonzero_integral_;  

  // Jets TF
  double CI_TFJet_;
  bool use_pT_TFJet_;

  // Selected processes 
  bool signalME_;
  double m_TauTau_2_;

  int MEversion_; //1 full version, 2 v2, 3 v3
  
  // Storage for useful quantities
  double   evTMinus_, evTPlus_;
  double   evSMinus_, evSPlus_;

  
 public :
  MGIntegration( );
  MGIntegration( RunConfig runConfig, int processID = -1);
  //MGIntegration( const char *configname, int processID = -1);
  //MGIntegration( const char *configname, RunConfig &runConfig, int processID = -1);
  ~MGIntegration( );
  
  // Config
  void setDefault( );
  void logConfig( ofstream & logStream);
  int  getVerbose( ) { return verbose_; };
  void setLogStream( ofstream &stream ) { stream_ = &stream; };
  
  const char*   getJacobianKindAsString( );
  
  // Setters of Event measured quantities
  void setLep_4Ps( TLorentzVector lep1, int lep1_Type, TLorentzVector lep2, int lep2_Type ) {

    if(nbrOfPermut_per_jet_==1){
      evLep1_4P_ = lep1; lepton1_Type_ = lep1_Type; lepton1_charge_ = lep1_Type>0 ? -1:1;
      evLep2_4P_ = lep2; lepton2_Type_ = lep2_Type; lepton2_charge_ = lep2_Type>0 ? -1:1;  
    }

    else{
      if(lep1.Pt()>lep2.Pt()){
	evLep1_4P_ = lep1; lepton1_Type_ = lep1_Type; lepton1_charge_ = lep1_Type>0 ? -1:1;
	evLep2_4P_ = lep2; lepton2_Type_ = lep2_Type; lepton2_charge_ = lep2_Type>0 ? -1:1;     
      }
      else{
	evLep1_4P_ = lep2; lepton1_Type_ = lep2_Type;  lepton1_charge_ = lep2_Type>0 ? -1:1;
	evLep2_4P_ = lep1; lepton2_Type_ = lep1_Type;  lepton2_charge_ = lep1_Type>0 ? -1:1;
      }
    }
    
  }

  void setLep_Tau4P( TLorentzVector tlv, int leptonType ) {
    evLep_Tau_4P_ = tlv; lepton_Tau_Type_ = leptonType; lepton_Tau_charge_ = leptonType>0 ? -1:1;
    evTMinus_ = TMath::Abs( t_minus ( evLep_Tau_4P_ ) );
    evTPlus_  = t_plus( evLep_Tau_4P_ );         
  }

  void setLep_top4P( TLorentzVector tlv, int leptonType ) {
    evLep_top_4P_ = tlv; lepton_top_Type_ = leptonType;  lepton_top_charge_ = leptonType>0 ? -1:1; 
  }

 
  void setHadronicTau4P( TLorentzVector tlv, int tauDecayMode ) {
    evHadSys_Tau_4P_ = tlv; HadtauDecayMode_ = tauDecayMode;
    evSMinus_ = t_minus ( evHadSys_Tau_4P_ ) ;
    evSPlus_  = t_plus( evHadSys_Tau_4P_ );     
  }
  
  void setlightJet_4Ps( TLorentzVector jet1,  TLorentzVector jet2 ) {

    if(jet1.Pt()>jet2.Pt()){
      evJet1_4P_ = jet1; evJet2_4P_ = jet2;
    }
    else{
      evJet1_4P_ = jet2; evJet2_4P_ = jet1;
    }
    n_lightJets_ = 2;
    nbrOfPermut_ = nbrOfPermut_per_jet_;
    //std::cout << "----- nbrOfPermut_per_jet_ inside setlightJet_4Ps : " << nbrOfPermut_per_jet_ << std::endl; // TEMPORAIRE //
    //std::cout << "----- nbrOfPermut_         inside setlightJet_4Ps : " << nbrOfPermut_ << std::endl; // TEMPORAIRE

  }

  void setlightJets_list_4P( const double jets[10][4], int nJets ) {
    
    n_lightJets_ = min(10,nJets);
    nbrOfPermut_ = nbrOfPermut_per_jet_;
    //std::cout << "----- nbrOfPermut_         inside setlightJets_list_4P : " << nbrOfPermut_ << std::endl; // TEMPORAIRE
    if(integration_type_ == integration_type_w_miss || force_missing_jet_integration_)
      nbrOfPermut_ *= nJets;    
    //std::cout << "----- nbrOfPermut_         inside setlightJets_list_4P : " << nbrOfPermut_ << std::endl; // TEMPORAIRE

    for(int i=0; i<nJets; i++){
      for(int j=0; j<4; j++)
	evJets_4P_[i][j] = jets[i][j];
    }
    if(n_lightJets_<10){
      for(int i=n_lightJets_; i<10; i++){
	for(int j=0; j<4; j++)
	  evJets_4P_[i][j] = 0;
      }
    }
    //std::cout << "----- nbrOfPermut_per_jet_ inside setlightJets_list_4P : " << nbrOfPermut_per_jet_ << std::endl; // TEMPORAIRE
    //std::cout << "----- nbrOfPermut_         inside setlightJets_list_4P : " << nbrOfPermut_ << std::endl; // TEMPORAIRE
      
  }
  
  void setBJet_4Ps( TLorentzVector jet1,  TLorentzVector jet2 ) {
    
    if(nbrOfPermut_per_jet_==1){
      evBJet1_4P_ = jet1; evBJet2_4P_ = jet2;
    }
    else{
      if(jet1.Pt()>jet2.Pt()){
	evBJet1_4P_ = jet1; evBJet2_4P_ = jet2;
      }
      else{
	evBJet1_4P_ = jet2; evBJet2_4P_ = jet1;
      }
    }
  }

  void setBJet_leptop_4P( TLorentzVector tlv ) {
    evBJet_leptop_4P_ = tlv;
  }
  
  void setBJet_hadtop_4P( TLorentzVector tlv ) {
    evBJet_hadtop_4P_ = tlv;
  }

 
  void setMETCov( const double V[4] ) {
    // Order : (0,0), (0,1), (1,0), (1,1)         
    evV_[0] = V[0]; evV_[1] = V[1]; 
    evV_[2] = V[2]; evV_[3] = V[3]; 
  }
  
  void setMET_4P( TLorentzVector tlv ) {
    evRecoMET4P_ = tlv;
  }


  // For permutations
  //Perm[0]: lep_top=lep1, lep_tau=lep2, b_leptop=BJet1, b_hadtop=BJet2
  //Perm[1]: lep_top=lep2, lep_tau=lep1, b_leptop=BJet1, b_hadtop=BJet2
  //Perm[2]: lep_top=lep1, lep_tau=lep2, b_leptop=BJet2, b_hadtop=BJet1
  //Perm[3]: lep_top=lep2, lep_tau=lep1, b_leptop=BJet2, b_hadtop=BJet1
  
  void initVersors(int perm){

    switch(perm){
    case 0:
      evLep_top_4P_ = evLep1_4P_;
      lepton_top_Type_ = lepton1_Type_;
      lepton_top_charge_ = lepton1_charge_;
      evLep_Tau_4P_ = evLep2_4P_;
      lepton_Tau_Type_ = lepton2_Type_;
      lepton_Tau_charge_ = lepton2_charge_;
      evBJet_leptop_4P_ = evBJet1_4P_;
      evBJet_hadtop_4P_ = evBJet2_4P_;
      break;

    case 1:
      evLep_top_4P_ = evLep2_4P_;
      lepton_top_Type_ = lepton2_Type_;
      lepton_top_charge_ = lepton2_charge_;
      evLep_Tau_4P_ = evLep1_4P_;
      lepton_Tau_Type_ = lepton1_Type_;
      lepton_Tau_charge_ = lepton1_charge_;
      evBJet_leptop_4P_ = evBJet1_4P_;
      evBJet_hadtop_4P_ = evBJet2_4P_;
      break;

    case 2:
      evLep_top_4P_ = evLep1_4P_;
      lepton_top_Type_ = lepton1_Type_;
      lepton_top_charge_ = lepton1_charge_;
      evLep_Tau_4P_ = evLep2_4P_;
      lepton_Tau_Type_ = lepton2_Type_;
      lepton_Tau_charge_ = lepton2_charge_;
      evBJet_leptop_4P_ = evBJet2_4P_;
      evBJet_hadtop_4P_ = evBJet1_4P_;
      break;
      
    case 3:
      evLep_top_4P_ = evLep2_4P_;
      lepton_top_Type_ = lepton2_Type_;
      lepton_top_charge_ = lepton2_charge_;
      evLep_Tau_4P_ = evLep1_4P_;
      lepton_Tau_Type_ = lepton1_Type_;
      lepton_Tau_charge_ = lepton1_charge_;
      evBJet_leptop_4P_ = evBJet2_4P_;
      evBJet_hadtop_4P_ = evBJet1_4P_;
      break;
    
    default:
      TLorentzVector dummy;
      evLep_top_4P_ = dummy;
      lepton_top_Type_ = 0;
      lepton_top_charge_ = 0;
      evLep_Tau_4P_ = dummy;
      lepton_Tau_Type_ = 0;
      lepton_Tau_charge_ = 0;
      evBJet_leptop_4P_ = dummy;
      evBJet_hadtop_4P_ = dummy;
      break;

    }
  }



  void initVersors_miss(int perm){

    int ijet=perm/4;
    if(ijet>n_lightJets_){
      cout<<"Less than "<<ijet+1<<" jets in the MEM event"<<endl;
      cout<<"EXITING"<<endl;
      return;
    }
    else
      evJet1_4P_ = evJets_4P_[ijet];

    switch(perm%4){
    case 0:
      evLep_top_4P_ = evLep1_4P_;
      lepton_top_Type_ = lepton1_Type_;
      lepton_top_charge_ = lepton1_charge_;
      evLep_Tau_4P_ = evLep2_4P_;
      lepton_Tau_Type_ = lepton2_Type_;
      lepton_Tau_charge_ = lepton2_charge_;
      evBJet_leptop_4P_ = evBJet1_4P_;
      evBJet_hadtop_4P_ = evBJet2_4P_;
      break;

    case 1:
      evLep_top_4P_ = evLep2_4P_;
      lepton_top_Type_ = lepton2_Type_;
      lepton_top_charge_ = lepton2_charge_;
      evLep_Tau_4P_ = evLep1_4P_;
      lepton_Tau_Type_ = lepton1_Type_;
      lepton_Tau_charge_ = lepton1_charge_;
      evBJet_leptop_4P_ = evBJet1_4P_;
      evBJet_hadtop_4P_ = evBJet2_4P_;
      break;

    case 2:
      evLep_top_4P_ = evLep1_4P_;
      lepton_top_Type_ = lepton1_Type_;
      lepton_top_charge_ = lepton1_charge_;
      evLep_Tau_4P_ = evLep2_4P_;
      lepton_Tau_Type_ = lepton2_Type_;
      lepton_Tau_charge_ = lepton2_charge_;
      evBJet_leptop_4P_ = evBJet2_4P_;
      evBJet_hadtop_4P_ = evBJet1_4P_;
      break;
      
    case 3:
      evLep_top_4P_ = evLep2_4P_;
      lepton_top_Type_ = lepton2_Type_;
      lepton_top_charge_ = lepton2_charge_;
      evLep_Tau_4P_ = evLep1_4P_;
      lepton_Tau_Type_ = lepton1_Type_;
      lepton_Tau_charge_ = lepton1_charge_;
      evBJet_leptop_4P_ = evBJet2_4P_;
      evBJet_hadtop_4P_ = evBJet1_4P_;
      break;
    
    default:
      TLorentzVector dummy;
      evLep_top_4P_ = dummy;
      lepton_top_Type_ = 0;
      lepton_top_charge_ = 0;
      evLep_Tau_4P_ = dummy;
      lepton_Tau_Type_ = 0;
      lepton_Tau_charge_ = 0;
      evBJet_leptop_4P_ = dummy;
      evBJet_hadtop_4P_ = dummy;
      break;

    }


  }
  

  // Integration
  int getIntegrationType() {return integration_type_; };

  int  getNbrOfDimttH()    { return nbrOfDim_ttH_; };
  int  getNbrOfDimttZ()    { return nbrOfDim_ttZ_; };
  int  getNbrOfDimttW()    { return nbrOfDim_ttW_; };
  int  getNbrOfDimttjets()    { return nbrOfDim_ttjets_; };
  int  getNbrOfDimttbar_SL()    { return nbrOfDim_ttbar_SL_; };
  int  getNbrOfDimttbar_DL()    { return nbrOfDim_ttbar_DL_; };

  int  getNbrOfDim_ttH_miss()    { return nbrOfDim_ttH_miss_; };
  int  getNbrOfDim_ttZ_miss()    { return nbrOfDim_ttZ_miss_; };
  int  getNbrOfDim_ttjets_miss()    { return nbrOfDim_ttjets_miss_; };
  int  getNbrOfDim_ttbar_SL_miss()    { return nbrOfDim_ttbar_SL_miss_; };
 
  int64_t getNbrOfPoints_ttH() { return nbrOfPoints_ttH_; };
  int64_t getNbrOfPoints_ttZ() { return nbrOfPoints_ttZ_; };
  int64_t getNbrOfPoints_ttW() { return nbrOfPoints_ttW_; };
  int64_t getNbrOfPoints_ttjets() { return nbrOfPoints_ttjets_; };
  int64_t getNbrOfPoints_ttbar_SL() { return nbrOfPoints_ttbar_SL_; };
  int64_t getNbrOfPoints_ttbar_DL() { return nbrOfPoints_ttbar_DL_; };
  int64_t getNbrOfPoints_ttZ_Zll() { return nbrOfPoints_ttZ_Zll_; };


  int64_t getNbrOfPoints_ttH_miss() { return nbrOfPoints_ttH_miss_; };
  int64_t getNbrOfPoints_ttZ_miss() { return nbrOfPoints_ttZ_miss_; };
  int64_t getNbrOfPoints_ttjets_miss() { return nbrOfPoints_ttjets_miss_; };
  int64_t getNbrOfPoints_ttbar_SL_miss() { return nbrOfPoints_ttbar_SL_miss_; };
  int64_t getNbrOfPoints_ttZ_Zll_miss() { return nbrOfPoints_ttZ_Zll_miss_; };

  // Integral Boundaries
  const double* getLowerValues() { return lowerValues_; };
  const double* getUpperValues() { return upperValues_; };  
  static double u_plus( double mTauTau_2, double cosThetaTauLepTauHad );
  static bool setCosThetaTauLepTauHadBoundaries( int leptonType,
              const TLorentzVector &Lep4P, const TLorentzVector &HadSys4P,
              double PTauLepBoundaries[2], double cosTauLepTauHadBoundaries[2], 
              double mTauTauHZ2);
//  bool setmTauTau_CosThetaTauLepTauHadBoundaries();
  void setmTauTau_CosThetaTauLepTauHadBoundaries();
  void setTauHadMomentumBoundaries();  

  double t_plus( const TLorentzVector &HadSys4P );
  double t_minus( const TLorentzVector &HadSys4P );
  static double t_minusLepton( const TLorentzVector &vec4P, int leptonType );
  static double t_plusLepton( const TLorentzVector &vec4P, int leptonType );

  pair<double,double> getE_Quark_Boundaries(TLorentzVector jet4P, TString flavor);
  pair<double,double> getE_fakelep_Boundaries(TLorentzVector lep4P);

  void setE_lightQuark_Boundaries();
  void checkCompatibility_TopHad();
  void checkCompatibility_TopHad_missing_jet();
  void checkCompatibility_TopLep();
  void checkCompatibility_TopTau();
  void checkCompatibility_TopLep_fakelep();
  void checkCompatibility_TopTau_fakelep();
  void checkCompatibility_Zll();

  // Helper functions for di-tau system  
  static double getPTauHad( double mTauTau_2, double PTauLep, double cosThetaTauLepTauHad );
  static double getGamma2( double PTauLep, double cosThetaTauLepTauHad, 
                         double cosThetaTauLepPi, const TLorentzVector &HadSys4P, double mTauTauHZ2 );
  static double getCosThetaTauHadPi( double P_TauHad, const TLorentzVector &HadSys4P); 
  static void getAlphaBetaGamma(double CosThetaTauLepTauHad, 
        double CosThetaTauLepPi, double CosThetaTauHadPi,
        double &alpha, double &beta, double &gamma, double &gamma2, const char **error);


  // Helper functions for tops
  double getEqbar_Enu(double CosTheta_qq, double Eq);
  double getEb(TLorentzVector W, TLorentzVector bjet);


  // Transfer functions
  vector<double> getJetTFresolution( double Egen, double eta, TString flavor) const;
  vector<double> getJetTFmean( double Egen, double eta, TString flavor) const;
  double getJetTFfraction( double eta, TString flavor) const;
  double getJetTF( const TLorentzVector &quark4P, const TLorentzVector &evJet4P, TString flavor) const;  

  double getFakeLepTFresolution(double pT_gen) const;
  double getFakeLepTFmean(double pT_gen) const;
  double getFakeLepTF( const TLorentzVector &quark4P, const TLorentzVector &evLep4P) const;

  double getFakeTauLepTFresolution(double pT_gen,int leptype) const;
  double getFakeTauLepTFmean(double pT_gen,int leptype) const;
  double getFakeTauLepTFn(double pT_gen,int leptype) const;
  double getFakeTauLepTF( const TLorentzVector &lep4P, const TLorentzVector &evTau4P, int leptype) const;


  double getMissJetAcceptance( const TLorentzVector &quark4P) const;

  double getTauLeptonicTF( double P_TauLep, const TLorentzVector &Lep4P ) const;  
  double getTauHadronicTF( double P_TauHad, const TLorentzVector &HadSys4P ) const; 
  double getTopTF(const TLorentzVector &t, const TLorentzVector &b, const TLorentzVector &l, const TLorentzVector &nu);
 

  // Jacobian
  double getJacobian_Tau(const TLorentzVector &tau) const;
  double getJacobian_diTau( const TLorentzVector &TauLep4P, const TLorentzVector &TauHad4P,
			   double gamma, double cosTheta_TauLepTauHad,  double sinTheta_TauLepPi) const ;
  double getJacobian_tophad(TLorentzVector b, TLorentzVector W, double Eqbar, double CosTheta_qq);
  double getJacobian_toplep(TLorentzVector b, TLorentzVector W, double El, double Enu);


  //ttjets specific functions
  double getFakeTF( const TLorentzVector &part4P, const TLorentzVector &evHadSys4P) const;
  double getPFake_upperBoundary_CI(double quantile,
				   const TLorentzVector &evHadSys4P);
  void setPOutFakeBoundaries();

  //ttW specific functions
  double getJacobian_W( double Enu, double Elep ) const;
  

  // Integrand utilities

  void setEventID( int64_t currentEvent_ ) { eventID_ = currentEvent_;}
  void setEventParameters( const IntegrationMsg_t &data, bool force_missing_jet_integration_forevent); 
  void copyBoundaries( IntegrationMsg_t *integration );

  
  // Integrands
  double evalttH(const double*);
  double evalttW(const double*);
  double evalttjets(const double*);
  double evalttbar_SL(const double*);
  double evalttbar_DL(const double*);
  double evalttbar_DL_fakelep_tlep(const double*);
  double evalttbar_DL_fakelep_ttau(const double*);
  double evalttZ_Zll(const double*);


 // Outputs
/*
  void writeFctValues(int64_t eventID, bool flagSignal,
		      double mTauTau2,
		      double cosTheta_miss_jet, double phi_miss_jet,
		      double PTauLep,  double cosTheta_diTau, 
		      double EQuark1,
		      double cosThetaNu, double phiNu,
		      double boost,
		      TLorentzVector rho4P, TLorentzVector boost4P, 
		      double TFJet1, double TFJet2,
		      double TFBjet_hadtop, double TFBjet_leptop,
		      double TFMET, 
		      double TFLepTau, double TFHadTau,
		      double TFLepTop, double TFHadTop,
		      double Jac, double wME, 
		      double eval,
		      const char *errStr) const;

  void writeFctValues_ttW(int64_t eventID, bool lep_sign,
			  double cosThetaNu_tlep, double phiNu_tlep,
			  double cosThetaNu_ttau, double phiNu_ttau,
			  double PTauHad,
			  double cosThetaNu_W, double phiNu_W,
			  double boost, 
			  TLorentzVector rho4P, TLorentzVector boost4P, 
			  double TFBjet_tautop, double TFBjet_leptop,
			  double TFMET, double TFHadTau,
			  double TFLepTop, double TFTauTop,
			  double Jac, double wME, 
			  double eval,
			  const char *errStr) const;


  void writeFctValues_ttjets(int64_t eventID,
			     double cosTheta_miss_jet, double phi_miss_jet,
			     double POutFake,
			     double EQuark1,
			     double cosThetaNu_tlep, double phiNu_tlep,
			     double boost, 
			     TLorentzVector rho4P, TLorentzVector boost4P, 
			     double TFJet1, double TFJet2,
			     double TFBjet_hadtop, double TFBjet_leptop,
			     double TFMET, 
			     double TFFake,
			     double TFLepTop, double TFHadTop,
			     double Jac, double wME, 
			     double eval,
			     const char *errStr) const;

  void writeFctValues_ttbar_SL(int64_t eventID,
			       double cosTheta_miss_jet, double phi_miss_jet,
			       double EQuark1,
			       double cosThetaNu_tlep, double phiNu_tlep,
			       double boost, 
			       TLorentzVector rho4P, TLorentzVector boost4P, 
			       double TFJet1, double TFJet2,
			       double TFBjet_hadtop, double TFBjet_leptop,
			       double TFMET, 
			       double TFLepTop, double TFHadTop,
			       double Jac, double wME, 
			       double eval,
			       const char *errStr) const;

  void writeFctValues_ttbar_DL(int64_t eventID,
			       double cosThetaNu_tlep, double phiNu_tlep,
			       double cosThetaNu_ttau, double phiNu_ttau,
			       double PTauHad,
			       double boost, 
			       TLorentzVector rho4P, TLorentzVector boost4P, 
			       double TFBjet_tautop, double TFBjet_leptop,
			       double TFMET, double TFHadTau,
			       double TFLepTop, double TFTauTop,
			       double Jac, double wME, 
			       double eval,
			       const char *errStr) const;

*/

  

};

extern "C" gsl_rng *gslRNGInit( );

extern "C" void gslIntegrate( gsl_monte_function *fct_descr, 
			      IntegrationMsg_t* integration, int nbrOfDim, int nbrOfPoints,
			      gsl_rng *rng, gsl_monte_vegas_state *state, bool newGrid,
			      double *res, double *err, double *chisq);

extern "C" double wrapper_evalttH(double *x, size_t dim, void* param);
extern "C" double wrapper_evalttW(double *x, size_t dim, void* param);
extern "C" double wrapper_evalttjets(double *x, size_t dim, void* param);
extern "C" double wrapper_evalttbar_SL(double *x, size_t dim, void* param);
extern "C" double wrapper_evalttbar_DL(double *x, size_t dim, void* param);
extern "C" double wrapper_evalttbar_DL_fakelep_tlep(double *x, size_t dim, void* param);
extern "C" double wrapper_evalttbar_DL_fakelep_ttau(double *x, size_t dim, void* param);
extern "C" double wrapper_evalttZ_Zll(double *x, size_t dim, void* param);


#endif	/* MGINTEGRATOR_H */

