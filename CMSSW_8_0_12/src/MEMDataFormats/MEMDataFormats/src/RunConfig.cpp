# include "MEMDataFormats/MEMDataFormats/interface/RunConfig.h"

using namespace std;

RunConfig::~RunConfig()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

  //std::cout << "Destructor RunConfig" << std::endl;

}

RunConfig::RunConfig(  ) { // const std::string configname

}

  
void RunConfig::print() const {

  cout << "\nConfig File Name : " << configName_ << endl;
  // 
  cout << "Event Type : " << eventType_<< endl;
  cout << "Max. number of events to process : " << maxNbrOfEventsToRead_ << endl;
  cout << "Input file names : " << endl;
  for ( uint i=0; i < inputFileNames_.size(); i++) {
      cout << "  " << i << " : " << inputFileNames_[i] << endl;
  }
  //
  cout << "File name of integral outputs           : " << integralsOutFileName_<< endl;
  //cout << "File name of function evaluation points : " << fctValuesOutFileName_<< endl;
  cout << "treeName: " << treeName_ << endl;
  //
  cout << "MadGraph card file name" << MGParamCardFileName_ << endl;
/*  cout << "File name of the One Prong LUT  : " << LUTOneProngPi0FileName_ << endl;
  cout << "File name of the Three Prong LUT: " << LUTThreeProngsFileName_ << endl;*/
  cout << "LHAPDF file name: " << LHAPDFFileName_ << endl;
  //
  cout << "sqrt(S): " << sqrtS_ << endl;
  cout << "Q_ (LHAPDF): " << Q_ << endl;
  //
  cout << "VEGAS" << endl;
  cout << "  # Of DimttH      : " << nbrOfDimttH_ << endl;
  cout << "  # Of DimttZ      : " << nbrOfDimttZ_ << endl;
  cout << "  # Of DimttW      : " << nbrOfDimttW_ << endl;
  cout << "  # Of Dimttjets   : " << nbrOfDimttjets_ << endl;
  cout << "  # Of Dimttbar_SL : " << nbrOfDimttbar_SL_ << endl;
  cout << "  # Of Dimttbar_DL : " << nbrOfDimttbar_DL_ << endl;
  cout << "  # Of DimttZ_Zll  : " << nbrOfDimttZ_Zll_ << endl;
  cout << "  # Of DimttH_miss       : " << nbrOfDimttH_miss_ << endl;
  cout << "  # Of DimttZ_miss       : " << nbrOfDimttZ_miss_ << endl;
  cout << "  # Of Dimttjets_miss    : " << nbrOfDimttjets_miss_ << endl;
  cout << "  # Of Dimttbar_SL_miss  : " << nbrOfDimttbar_SL_miss_ << endl;
  cout << "  # Of DimttZ_Zll_miss   : " << nbrOfDimttZ_Zll_miss_ << endl;

  cout << "  # Of runttZ_integration   : " << runttZ_integration_ << endl;
  cout << "  # Of runttW_integration   : " << runttW_integration_ << endl;
  cout << "  # Of runttjets_integration     : " << runttjets_integration_ << endl;
  cout << "  # Of runttbar_SL_integration   : " << runttbar_SL_integration_ << endl;
  cout << "  # Of runttbar_DL_integration   : " << runttbar_DL_integration_ << endl;
  cout << "  # Of runttbar_DL_fakelep_integration   : " << runttbar_DL_fakelep_integration_ << endl;
  cout << "  # Of runttZ_Zll_integration            : " << runttZ_Zll_integration_ << endl;

  cout << "  # Of run_missing_jet_integration              : " << run_missing_jet_integration_ << endl;
  cout << "  # Of force_missing_jet_integration            : " << force_missing_jet_integration_ << endl;
  cout << "  # Of force_missing_jet_integration_ifnoperm   : " << force_missing_jet_integration_ifnoperm_ << endl;

  cout << "  # Of Points_ttZ   : " << nbrOfPoints_ttZ_ << endl;
  cout << "  # Of Points_ttH   : " << nbrOfPoints_ttH_ << endl;
  cout << "  # Of Points_ttW   : " << nbrOfPoints_ttW_ << endl;
  cout << "  # Of Points_ttjets     : " << nbrOfPoints_ttjets_ << endl;
  cout << "  # Of Points_ttbar_SL   : " << nbrOfPoints_ttbar_SL_ << endl;
  cout << "  # Of Points_ttbar_DL   : " << nbrOfPoints_ttbar_DL_ << endl;
  cout << "  # Of Points_ttZ_Zll    : " << nbrOfPoints_ttZ_Zll_ << endl;

  cout << "  # Of Points_ttZ_miss        : " << nbrOfPoints_ttZ_miss_ << endl;
  cout << "  # Of Points_ttH_miss        : " << nbrOfPoints_ttH_miss_ << endl;
  cout << "  # Of Points_ttjets_miss     : " << nbrOfPoints_ttjets_miss_ << endl;
  cout << "  # Of Points_ttbar_SL_miss   : " << nbrOfPoints_ttbar_SL_miss_ << endl;
  cout << "  # Of Points_ttZ_Zll_miss    : " << nbrOfPoints_ttZ_Zll_miss_ << endl;

  cout << "  # Of Permut_per_jet      : " << nbrOfPermut_per_jet_ << endl;

  cout << "  Set the same RNG seed for all integrations: " << flagSameRNG_ << endl;
  cout << "  MadGraph reduced scheme: " << MEversion_ << endl; /**/ 

  cout << "  phiNu boundaries [ttH]           : " << phiNu_tlep_Boundaries_[0] << ", " << phiNu_tlep_Boundaries_[1] << endl;
  cout << "  cosThetaNu boundaries [ttH]      : " << cosThetaNu_tlep_Boundaries_[0] << ", " << cosThetaNu_tlep_Boundaries_[1] << endl;
  cout << "  phiNu_ttau boundaries [ttH]      : " << phiNu_ttau_ttbar_DL_ttW_Boundaries_[0] << ", " << phiNu_ttau_ttbar_DL_ttW_Boundaries_[1] << endl;
  cout << "  cosThetaNu_ttau boundaries [ttH] : " << cosThetaNu_ttau_ttbar_DL_ttW_Boundaries_[0] << ", " << cosThetaNu_ttau_ttbar_DL_ttW_Boundaries_[1] << endl;
  cout << "  phiNu_W_ttW boundaries [ttH]     : " << phiNu_W_ttW_Boundaries_[0] << ", " << phiNu_W_ttW_Boundaries_[1] << endl;
  cout << "  mTauTau cosThetaNu_W_ttW [ttH]   : " << cosThetaNu_W_ttW_Boundaries_[0] << ", " << cosThetaNu_W_ttW_Boundaries_[1] << endl;
 
  cout << "  phi_missing_jetBoundaries_       : " << phi_missing_jet_Boundaries_[0] << ", " << phi_missing_jet_Boundaries_[1] << endl;
  cout << "  cosTheta_missing_jetBoundaries_  : " << cosTheta_missing_jet_Boundaries_[0] << ", " << cosTheta_missing_jet_Boundaries_[1] << endl;

  cout << "  CI_TFJet         : " << CI_TFJet_ << endl;
  cout << "  use_pT_TFJet     : " << use_pT_TFJet_ << endl;
  cout << "  use_top_compatibility_check     : " << use_top_compatibility_check_ << endl;

  cout << "  include_hadrecoil          : " << include_hadrecoil_ << endl;
  cout << "  force_nonzero_integral     : " << force_nonzero_integral_ << endl;

  cout << "  eta_acceptance          : " << eta_acceptance_ << endl;
  cout << "  jet_radius              : " << jet_radius_ << endl;
  cout << "  dR_veto_jet_lep         : " << dR_veto_jet_lep_ << endl;
  cout << "  rel_iso_lep             : " << rel_iso_lep_ << endl;
  cout << "  pT_cut                  : " << pT_cut_ << endl;

  cout << "Function value flags" << endl;
  cout << "  flagTFLepTau     : " << flagTFLepTau_ << endl;  
  cout << "  flagTFHadTau     : " << flagTFHadTau_ << endl;
  cout << "  flagTFFake       : " << flagTFFake_ << endl;
  cout << "  flagTFMET        : " << flagTFMET_ << endl;
  cout << "  flagTFJet1       : " << flagTFJet1_ << endl;
  cout << "  flagTFJet2       : " << flagTFJet2_ << endl;
  cout << "  flagTFBJet_leptop       : " << flagTFBJet_leptop_ << endl;
  cout << "  flagTFBJet_hadtop       : " << flagTFBJet_hadtop_ << endl;
  cout << "  flagTF_fakelep          : " << flagTF_fakelep_ << endl;
  cout << "  flagTF_fakeleptau       : " << flagTF_fakeleptau_ << endl;
  cout << "  flagTFTop        : " << flagTFTop_ << endl;
  cout << "  flagJac          : " << flagJac_ << endl;
  cout << "  flagWME          : " << flagWME_ << endl;

  cout << "Verbose mode : " << verbose_ << endl;
  cout << "eventType_ : " << eventType_ << endl;
  
}
