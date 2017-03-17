// -*- C++ -*-
//
// Package:    MEM/MEMProducer
// Class:      MEMProducer
// 
/**\class MEMProducer MEMProducer.cc MEM/MEMProducer/plugins/MEMProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation] 
*/
//
// Original Author:  Arnaud Chiron
//         Created:  Mon, 15 Feb 2016 10:41:27 GMT
//     ttH Version:  Tue, 29 Nov 2016 17:32:26 GMT
//
//

//# include "mpi.h"

// system include files
#include "MEM/MEMAlgo/interface/ThreadScheduler.h"
#include "MEM/MEMProducer/plugins/MEMProducer.h"
#include "MEMDataFormats/MEMDataFormats/interface/IntegralResult.h"
#include "MEMDataFormats/MEMDataFormats/interface/MEMResult.h"

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "CondFormats/DataRecord/interface/JetResolutionRcd.h" // not used ??
#include "FWCore/Framework/interface/eventsetuprecord_registration_macro.h" // not used ??

using namespace reco;
RunConfig runConfig;
NodeScheduler *scheduler; 

//
// constants, enums and typedefs
//
# define debugLEVEL 2 // from MPIScheduler.h

//static int nbrOfProcess_;
//static int processID_;
//static int namelen_ ;
//static char processor_name[MPI_MAX_PROCESSOR_NAME];

//
// static data member definitions
//

bool MEMProducer::tauFilter(const pat::Tau& tauFiltre)
{
    bool b (true);
    bool passPtSelection(true), passEtaSelection(true) ;
    pat::Tau f_tau(tauFiltre);
    passPtSelection = tauFiltre.pt() > 20.;
    passEtaSelection = fabs(tauFiltre.eta()) < 2.3;
    
    /*if ( !passPtSelection ) { // Pt <= 18. 
        b = false;
	}*/
    
    if ( !(passPtSelection && passEtaSelection) ) {
        b = false;
	}

    /*int tauid1 = (f_tau.isTauIDAvailable("byTightCombinedIsolationDeltaBetaCorr3Hits") ? f_tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits") : -999);
    int tauid2 = (f_tau.isTauIDAvailable("decayModeFindingNewDMs") ? f_tau.tauID("decayModeFindingNewDMs") : -999);
    int tauid3 = (f_tau.isTauIDAvailable("againstMuonTight3") ? f_tau.tauID("againstMuonTight3") : -999);
    int tauid4 = (f_tau.isTauIDAvailable("againstElectronVLooseMVA5") ? f_tau.tauID("againstElectronVLooseMVA5") : -999);
    if (tauid1 != 1 || tauid1 == -999) { 
        b = false;
    }
    if (tauid2 <= 0 || tauid2 == -999) { 
        b = false;
    }
    if (tauid3 <= 0 || tauid3 == -999) { 
        b = false;
    }
    if (tauid4 <= 0 || tauid4 == -999) { 
        b = false;
	}*/

    int tauid1 = (f_tau.isTauIDAvailable("decayModeFinding") ? f_tau.tauID("decayModeFinding") : -999);
    int tauid2 = (f_tau.isTauIDAvailable("byLooseIsolationMVArun2v1DBdR03oldDMwLT") ? f_tau.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT") : -999);
    if (tauid1 != 1 || tauid1 == -999) { 
      b = false;
    }/**/
    if (tauid2 <= 0 || tauid2 == -999) { 
      b = false;
    }
    
    return b;
}

bool MEMProducer::muonFilter(const pat::Muon& muonFiltre)
{
    bool b (false);
    bool passPtSelection(true), passEtaSelection(true) ;
    passPtSelection = muonFiltre.pt() > 8.;
    passEtaSelection = fabs(muonFiltre.eta())<2.4 ;
    
    /*if ( passPtSelection && muonFiltre.isLooseMuon() ) { // Pt <= 8. && isLooseMuon()
        b = true;
    }*/
    /*if ( muonFiltre.chargedHadronIso() >= 0.1 ) { // Iso <= 0.1
    //if ( muonFiltre.chargedHadronIso()/muonFiltre.pt() >= 0.1 ) { // Iso <= 0.1
        b = false;
	}*/

    if ( passPtSelection && passEtaSelection && muonFiltre.isLooseMuon() ) {
        b = true;
    }

    return b;
}

bool MEMProducer::electronFilter(const pat::Electron& electronFiltre)
{
    bool b (false);
    bool passPtSelection(true), passEtaSelection(true) ;
    passPtSelection = electronFiltre.pt() > 8.;
    passEtaSelection = fabs(electronFiltre.eta()) < 2.5 ;
    
    /*if ( passPtSelection && muonFiltre.isLooseMuon() ) { // Pt <= 8. && isLooseMuon()
        b = true;
    }*/
    /*if ( muonFiltre.chargedHadronIso() >= 0.1 ) { // Iso <= 0.1
    //if ( muonFiltre.chargedHadronIso()/muonFiltre.pt() >= 0.1 ) { // Iso <= 0.1
        b = false;
	}*/

    if ( passPtSelection && passEtaSelection ) {
        b = true;
    }

    return b;
}

bool MEMProducer::jetFilter(const pat::Jet& jetFiltre)
{
    bool b (false);   
    bool passPtSelection(true), passEtaSelection(true) ;
    passPtSelection = jetFiltre.pt() > 25.;
    passEtaSelection = fabs(jetFiltre.eta())<2.4;
    float NHF = jetFiltre.neutralHadronEnergyFraction();
    float NEMF = jetFiltre.neutralEmEnergyFraction();
    float CHF = jetFiltre.chargedHadronEnergyFraction();
    float CEMF = jetFiltre.chargedEmEnergyFraction();
    int NumNeutralParticles =jetFiltre.neutralMultiplicity();
    int chargedMult = jetFiltre.chargedMultiplicity();
    int NumConst = chargedMult + NumNeutralParticles;
    float CHM = jetFiltre.chargedMultiplicity();
    bool looseJetID = ( (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((passEtaSelection && CHF>0 && CHM>0 && CEMF<0.99) || !passEtaSelection ) ); //   absjeta>2.4) );
  
   if ( passPtSelection && passEtaSelection && looseJetID ) {
        b = true;
    }
    
    
    return b;
}

//Returns a vector of cleaned jets
//Empty collection if event does not pass basic filters
//May probably have to be re-optimized in the future

//Implementation for ele+ele
//Similar implementation to be used for ele+mu and mu+mu
/*vector<pat::Jet> MEMProducer::npletFilter_ee(const pat::Electron& ele1, const pat::Electron& ele2, const pat::Tau& tau, edm::Handle<pat::JetCollection> jets){

  vector<pat::Jet> cleaned_jets;
  
  std::cout << "inside npletFilter_ee" << std::endl;

  math::XYZTLorentzVector ele1_P4 = ele1.p4();
  math::XYZTLorentzVector ele2_P4 = ele2.p4();
  math::XYZTLorentzVector tau_P4 = tau.p4();

  if ( !(deltaR(ele1_P4,ele2_P4)>0.3 && deltaR(ele1_P4,tau_P4)>0.3 && deltaR(ele2_P4,tau_P4)>0.3) )
    return cleaned_jets ; // null cleaned_jets_ordered

  for (const pat::Jet &jet : *jets ) {
    //std::cout << "\t\t jet eta : " << jet.eta() << std::endl;
    math::XYZTLorentzVector jet_P4 = jet.p4();
    if(jetFilter(jet) && deltaR(ele1_P4,jet_P4)>0.4 && deltaR(ele2_P4,jet_P4)>0.4 && deltaR(tau_P4,jet_P4)>0.4)
      cleaned_jets.push_back(jet);

  }

  return cleaned_jets;

}*/

/*vector<pat::Jet> MEMProducer::npletFilter_mm(const pat::Muon& muon1, const pat::Muon& muon2, const pat::Tau& tau, edm::Handle<pat::JetCollection> jets){

  vector<pat::Jet> cleaned_jets;
  
  std::cout << "inside npletFilter_mm" << std::endl;

  math::XYZTLorentzVector muon1_P4 = muon1.p4();
  math::XYZTLorentzVector muon2_P4 = muon2.p4();
  math::XYZTLorentzVector tau_P4 = tau.p4();

  if ( !(deltaR(muon1_P4,muon2_P4)>0.3 && deltaR(muon1_P4,tau_P4)>0.3 && deltaR(muon2_P4,tau_P4)>0.3) )
    return cleaned_jets ; // null cleaned_jets_ordered

  for (const pat::Jet &jet : *jets ) {
    math::XYZTLorentzVector jet_P4 = jet.p4();
    if(jetFilter(jet) && deltaR(muon1_P4,jet_P4)>0.4 && deltaR(muon2_P4,jet_P4)>0.4 && deltaR(tau_P4,jet_P4)>0.4)
      cleaned_jets.push_back(jet);
  }

  return cleaned_jets;

}*/

/*vector<pat::Jet> MEMProducer::npletFilter_em(const pat::Electron& ele1, const pat::Muon& muon2, const pat::Tau& tau, edm::Handle<pat::JetCollection> jets){

  vector<pat::Jet> cleaned_jets;
  
  std::cout << "inside npletFilter_em" << std::endl;

  math::XYZTLorentzVector ele1_P4 = ele1.p4();
  math::XYZTLorentzVector muon2_P4 = muon2.p4();
  math::XYZTLorentzVector tau_P4 = tau.p4();

  if ( !(deltaR(ele1_P4,muon2_P4)>0.3 && deltaR(ele1_P4,tau_P4)>0.3 && deltaR(muon2_P4,tau_P4)>0.3) )
    return cleaned_jets ; // null cleaned_jets_ordered

  for (const pat::Jet &jet : *jets ) {
    math::XYZTLorentzVector jet_P4 = jet.p4();
    if(jetFilter(jet) && deltaR(ele1_P4,jet_P4)>0.4 && deltaR(muon2_P4,jet_P4)>0.4 && deltaR(tau_P4,jet_P4)>0.4)
      cleaned_jets.push_back(jet);
  }

  return cleaned_jets;

}*/

//
// constructors and destructor
//
MEMProducer::MEMProducer(const edm::ParameterSet& iConfig)
{
   //now do what ever other initialization is needed
   
//   MPI_Init(0, 0);
//   MPI_Comm_size(MPI_COMM_WORLD, &nbrOfProcess_);
//   MPI_Comm_rank(MPI_COMM_WORLD, &processID_ );
//   MPI_Get_processor_name(processor_name,&namelen_ );

/*  std::cout 
    << "cout:Constructor nbrOfProcess = " << nbrOfProcess_ 
    << "  processID = " << processID_ 
    << std::endl;*/
  
    electronToken_ = consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons")); // OK
    muonToken_     = consumes<pat::MuonCollection>    (iConfig.getParameter<edm::InputTag>("muons"));
    tauToken_      = consumes<pat::TauCollection>     (iConfig.getParameter<edm::InputTag>("taus"));
    metToken_      = consumes<pat::METCollection>     (iConfig.getParameter<edm::InputTag>("mets"));
    jetToken_      = consumes<pat::JetCollection>     (iConfig.getParameter<edm::InputTag>("jets"));/**/
    
    theSigTag      = consumes<double>                 (iConfig.getParameter<edm::InputTag>("srcSig"));
    theCovTag      = consumes<math::Error<2>::type>   (iConfig.getParameter<edm::InputTag>("srcCov"));
      
    produces<pat::ElectronCollection>();
    produces<pat::MuonCollection>();
    produces<pat::METCollection>();
    produces<pat::TauCollection>();
    produces<pat::JetCollection>();

    produces<reco::MEMResultCollection>();
    
    std::cout << "debut appel RunConfig" << std::endl;                               
    std::string str;
    
    edm::ParameterSet histosSet = iConfig.getParameter<edm::ParameterSet>("parameters") ;
    
    runConfig.configName_ = iConfig.getParameter<std::string>("configName") ;
    runConfig.maxNbrOfEventsToRead_ = histosSet.getParameter<int>("maxNbrOfEventsToRead") ; // 
    if  (runConfig.maxNbrOfEventsToRead_ < 0) 
        runConfig.maxNbrOfEventsToRead_ = LLONG_MAX;

    runConfig.sqrtS_ = histosSet.getParameter<double>("sqrtS") ; // 
    
    str = histosSet.getParameter<std::string>("InputType") ; //
    std::cout << "str : " << str << std::endl;
    if(str == "MonteCarlo") {
        runConfig.eventType_ = MCEventType;
    } else if (str == "Run1") {
        runConfig.eventType_ = Run1EventType;
    } else if (str == "Run1_etau") {
        runConfig.eventType_ = Run1_etau_EventType;
    } else if (str == "PyRun1") {
        runConfig.eventType_ = PyRun1EventType;
    } else if (str == "PyRun2") {
        runConfig.eventType_ = PyRun2EventType;
    } else {
        cerr << "PyConfig: Bad input filename type (MonteCarlo/Run1/PyRun1)" << endl; 
        runConfig.eventType_ = BadEventType;
    }/**/
    
    runConfig.Q_ = histosSet.getParameter<double>("Q"); //
    runConfig.flagSameRNG_ = histosSet.getParameter<bool>("flagSameRNG"); //
    runConfig.flagTFLepTau_ = histosSet.getParameter<bool>("flagTFLepTau"); //
    runConfig.flagTFHadTau_ = histosSet.getParameter<bool>("flagTFHadTau"); //
    runConfig.flagTFFake_ =  histosSet.getParameter<bool>("flagTFFake"); //
    runConfig.flagTFMET_ = histosSet.getParameter<bool>("flagTFMET"); //
    runConfig.flagTFJet1_ = histosSet.getParameter<bool>("flagTFJet1"); //
    runConfig.flagTFJet2_ = histosSet.getParameter<bool>("flagTFJet2"); //
    runConfig.flagTFBJet_leptop_ = histosSet.getParameter<bool>("flagTFBJet_leptop"); //
    runConfig.flagTFBJet_hadtop_ = histosSet.getParameter<bool>("flagTFBJet_hadtop"); //
    runConfig.flagTF_fakelep_ = histosSet.getParameter<bool>("flagTF_fakelep"); //
    runConfig.flagTF_fakeleptau_ = histosSet.getParameter<bool>("flagTF_fakeleptau"); //
    runConfig.flagTFTop_ = histosSet.getParameter<bool>("flagTFTop"); //
    runConfig.flagJac_ = histosSet.getParameter<bool>("flagJac"); //
    runConfig.flagWME_ = histosSet.getParameter<bool>("flagWME"); //
    runConfig.inputFileNames_ = histosSet.getParameter<std::vector<std::string> >("InputFileList"); //
    runConfig.verbose_ = histosSet.getParameter<int>("verbose") ; //
    runConfig.mTau_ = histosSet.getParameter<double>("mTau"); //
    runConfig.mHiggs_ = histosSet.getParameter<double>("mHiggs"); //
    runConfig.mZ_ = histosSet.getParameter<double>("mZ"); //
    runConfig.MEversion_ = histosSet.getParameter<int>("MEversion"); //
    runConfig.integralsOutFileName_ = histosSet.getParameter<std::string>("FileOfIntegrals") ; //
//    runConfig.fctValuesOutFileName_ = histosSet.getParameter<std::string>("FileOfFunctionValues") ; //
    runConfig.nbrOfDimttH_ = histosSet.getParameter<int>("nbrOfDimttH");
    runConfig.nbrOfDimttZ_ = histosSet.getParameter<int>("nbrOfDimttZ");
    runConfig.nbrOfDimttW_ = histosSet.getParameter<int>("nbrOfDimttW");
    runConfig.nbrOfDimttjets_ = histosSet.getParameter<int>("nbrOfDimttjets");
    runConfig.nbrOfDimttbar_SL_ = histosSet.getParameter<int>("nbrOfDimttbar_SL");
    runConfig.nbrOfDimttbar_DL_ = histosSet.getParameter<int>("nbrOfDimttbar_DL");
    runConfig.nbrOfDimttZ_Zll_ = histosSet.getParameter<int>("nbrOfDimttZ_Zll");
    runConfig.nbrOfDimttH_miss_ = histosSet.getParameter<int>("nbrOfDimttH_miss");
    runConfig.nbrOfDimttZ_miss_ = histosSet.getParameter<int>("nbrOfDimttZ_miss");
    runConfig.nbrOfDimttjets_miss_ = histosSet.getParameter<int>("nbrOfDimttjets_miss");
    runConfig.nbrOfDimttbar_SL_miss_ = histosSet.getParameter<int>("nbrOfDimttbar_SL_miss");
    runConfig.nbrOfDimttZ_Zll_miss_ = histosSet.getParameter<int>("nbrOfDimttZ_Zll_miss");
    runConfig.runttZ_integration_ = histosSet.getParameter<bool>("runttZ_integration");
    runConfig.runttW_integration_ = histosSet.getParameter<bool>("runttW_integration");
    runConfig.runttjets_integration_ = histosSet.getParameter<bool>("runttjets_integration");
    runConfig.runttbar_SL_integration_ = histosSet.getParameter<bool>("runttbar_SL_integration");
    runConfig.runttbar_DL_integration_ = histosSet.getParameter<bool>("runttbar_DL_integration");
    runConfig.runttbar_DL_fakelep_integration_ = histosSet.getParameter<bool>("runttbar_DL_fakelep_integration");
    runConfig.runttZ_Zll_integration_ = histosSet.getParameter<bool>("runttZ_Zll_integration");

    runConfig.run_missing_jet_integration_ = histosSet.getParameter<bool>("run_missing_jet_integration");
    runConfig.force_missing_jet_integration_ = histosSet.getParameter<bool>("force_missing_jet_integration");
    runConfig.force_missing_jet_integration_ifnoperm_ = histosSet.getParameter<bool>("force_missing_jet_integration_ifnoperm");

    runConfig.nbrOfPoints_ttZ_ = histosSet.getParameter<int>("nbrOfPoints_ttZ");
    runConfig.nbrOfPoints_ttH_ = histosSet.getParameter<int>("nbrOfPoints_ttH");
    runConfig.nbrOfPoints_ttW_ = histosSet.getParameter<int>("nbrOfPoints_ttW");
    runConfig.nbrOfPoints_ttjets_ = histosSet.getParameter<int>("nbrOfPoints_ttjets");
    runConfig.nbrOfPoints_ttbar_SL_ = histosSet.getParameter<int>("nbrOfPoints_ttbar_SL");
    runConfig.nbrOfPoints_ttbar_DL_ = histosSet.getParameter<int>("nbrOfPoints_ttbar_DL");
    runConfig.nbrOfPoints_ttZ_Zll_ = histosSet.getParameter<int>("nbrOfPoints_ttZ_Zll");

    runConfig.nbrOfPoints_ttZ_miss_ = histosSet.getParameter<int>("nbrOfPoints_ttZ_miss");
    runConfig.nbrOfPoints_ttH_miss_ = histosSet.getParameter<int>("nbrOfPoints_ttH_miss");
    runConfig.nbrOfPoints_ttjets_miss_ = histosSet.getParameter<int>("nbrOfPoints_ttjets_miss");
    runConfig.nbrOfPoints_ttbar_SL_miss_ = histosSet.getParameter<int>("nbrOfPoints_ttbar_SL_miss");
    runConfig.nbrOfPoints_ttZ_Zll_miss_ = histosSet.getParameter<int>("nbrOfPoints_ttZ_Zll_miss");

    runConfig.nbrOfPermut_per_jet_ = histosSet.getParameter<int>("nbrOfPermut_per_jet");

    runConfig.phiNu_tlep_Boundaries_ = histosSet.getParameter<std::vector<double> >("phiNu_tlep");
    runConfig.cosThetaNu_tlep_Boundaries_ = histosSet.getParameter<std::vector<double> >("cosThetaNu_tlep");
    runConfig.phiNu_ttau_ttbar_DL_ttW_Boundaries_ = histosSet.getParameter<std::vector<double> >("phiNu_ttau");
    runConfig.cosThetaNu_ttau_ttbar_DL_ttW_Boundaries_ = histosSet.getParameter<std::vector<double> >("cosThetaNu_ttau");
    runConfig.phiNu_W_ttW_Boundaries_ = histosSet.getParameter<std::vector<double> >("phiNu_W_ttW");
    runConfig.cosThetaNu_W_ttW_Boundaries_ = histosSet.getParameter<std::vector<double> >("cosThetaNu_W_ttW");
 
    runConfig.CI_TFJet_ = histosSet.getParameter<double>("CI_TFJet");
    runConfig.use_pT_TFJet_ = histosSet.getParameter<bool>("use_pT_TFJet");
    runConfig.use_top_compatibility_check_ = histosSet.getParameter<bool>("use_top_compatibility_check");

    runConfig.include_hadrecoil_ = histosSet.getParameter<bool>("include_hadrecoil");
    runConfig.force_nonzero_integral_ = histosSet.getParameter<bool>("force_nonzero_integral");

    runConfig.eta_acceptance_ = histosSet.getParameter<double>("eta_acceptance");
    runConfig.jet_radius_ = histosSet.getParameter<double>("jet_radius");
    runConfig.dR_veto_jet_lep_ = histosSet.getParameter<double>("dR_veto_jet_lep");
    runConfig.rel_iso_lep_ = histosSet.getParameter<double>("rel_iso_lep");
    runConfig.pT_cut_ = histosSet.getParameter<double>("pT_cut");

    runConfig.phi_missing_jet_Boundaries_ = histosSet.getParameter<std::vector<double> >("phi_missing_jet");
    runConfig.cosTheta_missing_jet_Boundaries_ = histosSet.getParameter<std::vector<double> >("cosTheta_missing_jet");

    runConfig.treeName_ = histosSet.getParameter<std::string>("treeName"); // VBF : "No Tree Name"; //
    runConfig.LHAPDFFileName_ = histosSet.getParameter<std::string>("LHAPDFFile") ; //
    runConfig.MGParamCardFileName_ = histosSet.getParameter<std::string>("MGParamCardFile") ; //

    runConfig.print();
    std::cout << "fin appel RunConfig" << std::endl;

    /* scheduler initialization */
    scheduler = new ThreadScheduler();
    std::cout << "fin creation ThreadScheduler" << std::endl;
        
}

MEMProducer::~MEMProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

//  MPI_Finalize();
  std::cout 
//    << "Destructor nbrOfProcess = " << nbrOfProcess_ 
//    << "  processID = " << processID_ 
    << std::endl;

}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
MEMProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm; 
    using namespace reco; 
    using namespace std;
//    char greeting[100];
    char tr_1[100];
//    MPI_Status stat;
    
    // get collections
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);
    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);
    edm::Handle<pat::TauCollection> taus;
    iEvent.getByToken(tauToken_, taus);
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);/**/

    std::auto_ptr<reco::MEMResultCollection> MEMResultColl( new reco::MEMResultCollection );

//    sprintf(greeting, "Hello world: processor %d of %d\n", processID_, nbrOfProcess_);

    //    the lines below are correct
//    if ( processID_ == 0 ) { 
        //
        //  Master
        //

//        std::cout << "\nTreating event "<< iEvent.id() << std::endl; 
        
        /*for (int partner = 1; partner < nbrOfProcess_; partner++){
            MPI_Recv(greeting, sizeof(greeting), MPI_BYTE, partner, 1, MPI_COMM_WORLD, &stat);
            fputs (greeting, stdout); // necessaire pour envoyer sur les noeuds
        } */

        // Output collection
        auto_ptr<pat::ElectronCollection> result_elec( new pat::ElectronCollection() );
        auto_ptr<pat::MuonCollection>     result_muon( new pat::MuonCollection() );
        auto_ptr<pat::METCollection>      result_met ( new pat::METCollection() );
        auto_ptr<pat::TauCollection>      result_tau ( new pat::TauCollection() );
        auto_ptr<pat::JetCollection>      result_jet ( new pat::JetCollection() );
      
        /* boucle générale */
        sprintf(tr_1, "----- Event: %lld \n", iEvent.id().event());
        sprintf(tr_1, "----- Lumi : %i   \n", iEvent.id().luminosityBlock());
        sprintf(tr_1, "----- Run  : %i   \n", iEvent.id().run());
        
        runConfig.nbOfElectrons_ = 0;
        runConfig.nbOfMuons_ = 0;
        runConfig.nbOfTaus_ = 0;

        /* scheduler initialization */
        scheduler->initNodeScheduler( runConfig, 0 );
        
        /* BOUCLE GENERALE */
        
        sprintf(tr_1, "----- Event: %lld \n", iEvent.id().event());
        sprintf(tr_1, "----- Lumi : %i   \n", iEvent.id().luminosityBlock());
        sprintf(tr_1, "----- Run  : %i   \n", iEvent.id().run());

        /* Muons - First Loop */
        for (const pat::Muon &muon : *muons ) {
            if ( muonFilter(muon) ){
                pat::Muon l_muon(muon);
                //std::cout << "1st loop : Muon ID : " << &muon - &(*muons->begin()) << " - eta : " << muon.eta() << std::endl ;
                edm::LogInfo("MEMProducer::produce") << "2nd loop : Muon ID : " << &muon - &(*muons->begin()) ;
                result_muon->push_back(l_muon);
                runConfig.nbOfMuons_ += 1;
                runConfig.ttlnbOfMuons_ += 1;
            }/**/ // muons filter
        } // muon loop
        
        /* Electrons ^- Second Loop */
        for (const pat::Electron &elec : *electrons ) {
            if ( electronFilter(elec) ){
                pat::Electron l_elec(elec);
                //std::cout << "2nd loop : Electron ID : " << &elec - &(*electrons->begin()) << " - eta : " << l_elec.eta() << std::endl ;
                edm::LogInfo("MEMProducer::produce") << "2nd loop : Electron ID : " << &elec - &(*electrons->begin()) ;
                l_elec.setCharge(-7);
                result_elec->push_back(l_elec);
                runConfig.nbOfElectrons_ += 1;
                runConfig.ttlnbOfElectrons_ += 1;
            } // electron filter
        }/**/ // electrons loop
        
        /* Taus third loop */
        for (const pat::Tau &tau : *taus ) {
            if ( tauFilter(tau) ){
                pat::Tau l_tau(tau);
                //std::cout << "3rd loop : Tau ID : " << &tau - &(*taus->begin()) << " - eta : " << l_tau.eta() << std::endl ;
                edm::LogInfo("MEMProducer::produce") << "2nd loop : Tau ID : " << &tau - &(*taus->begin()) ;
                result_tau->push_back(l_tau);
                runConfig.nbOfTaus_ += 1;
                runConfig.ttlnbOfTaus_ += 1;
            } // tau filter
        } // tau loop /**/
                
        std::cout << "End of particle search in Event : " << iEvent.id().event() << std::endl ;
        
        std::cout << "\t Working on particles" << std::endl;
        std::cout << "nb of Electrons : " << runConfig.nbOfElectrons_ << std::endl;
        std::cout << "nb of Muons : " << runConfig.nbOfMuons_ << std::endl;
        std::cout << "nb of Taus  : " << runConfig.nbOfTaus_ << std::endl;      
        
        bool pairs_ee = ( runConfig.nbOfElectrons_ > 1 ? true : false );
        bool pairs_mm = ( runConfig.nbOfMuons_ > 1 ? true : false );
        bool pairs_em = ( ( ( runConfig.nbOfElectrons_ > 0 ) && (runConfig.nbOfMuons_ > 0) ) ? true : false );
        std::cout << "pairs_ee : " << pairs_ee << std::endl;
        std::cout << "pairs_mm : " << pairs_mm << std::endl;
        std::cout << "pairs_em : " << pairs_em << std::endl;
        
        /* searching for pairs */
        if (runConfig.nbOfTaus_ > 1) { // if = 1, only one tau, no pair
            std::cout << "you will have " << runConfig.nbOfTaus_  << " possibilities for taus" << std::endl;
        }
        
        /* muon-muon case */
        if ( pairs_mm ) { // 
            std::cout << "you will have " << runConfig.nbOfMuons_ * (runConfig.nbOfMuons_ - 1) / 2 << " for muon-muon pairs" << std::endl;
            for ( pat::MuonCollection::const_iterator iter_3 = result_muon->begin(); iter_3 < result_muon->end() - 1 ; ++iter_3) {
                for ( pat::MuonCollection::const_iterator iter_4 = iter_3 + 1; iter_4 != result_muon->end(); ++iter_4) {
                    std::cout << "(" << (iter_3-result_muon->begin())  << "," << (iter_4-result_muon->begin()) << ")" << std::endl;
                    for (const pat::Tau &tau2 : *result_tau ) {
                        runConfig.nbOfTriplets_ += 1;
                        MEM_nplet nplet(*iter_3, *iter_4, tau2); //
                        vector<pat::Jet> cleanedJets = npletFilter<pat::Muon, pat::Muon>( *iter_3, *iter_4, tau2, jets );
                        JetFilling<pat::Muon, pat::Muon>( *iter_3, *iter_4, tau2, cleanedJets, nplet, "mm" );
                        runConfig.nbOfNplets_ += 1;
                        std::cout << "classical : nbOfNplets in mm : " << runConfig.nbOfNplets_ << std::endl;
                    }
                }
            }
            
        }
        
        /* electron-electron case */
        if ( pairs_ee ) { // 
            std::cout << "you will have " << runConfig.nbOfElectrons_ * (runConfig.nbOfElectrons_ - 1) / 2 << " for electron-electron pairs" << std::endl;
            for ( pat::ElectronCollection::const_iterator iter_1 = result_elec->begin(); iter_1 < result_elec->end() - 1 ; ++iter_1) {
                for ( pat::ElectronCollection::const_iterator iter_2 = iter_1 + 1; iter_2 != result_elec->end(); ++iter_2) {
                    std::cout << "(" << (iter_1-result_elec->begin())  << "," << (iter_2-result_elec->begin()) << ")" << std::endl;
                    for (const pat::Tau &tau2 : *result_tau ) {
                        runConfig.nbOfTriplets_ += 1;
                        MEM_nplet nplet(*iter_1, *iter_2, tau2); //
                        vector<pat::Jet> cleanedJets = npletFilter<pat::Electron, pat::Electron>( *iter_1, *iter_2, tau2, jets );
                        JetFilling<pat::Electron, pat::Electron>( *iter_1, *iter_2, tau2, cleanedJets, nplet, "ee" );
                        runConfig.nbOfNplets_ += 1;
                        std::cout << "classical : nbOfNplets in ee : " << runConfig.nbOfNplets_ << std::endl;
                    }
                }
            }
            
        }
        
        /* electron-muon case */
        if ( pairs_em ) { // 
            std::cout << "you will have " << runConfig.nbOfElectrons_ * runConfig.nbOfMuons_ << " for electron-muon pairs" << std::endl;
            for ( pat::ElectronCollection::const_iterator iter_5 = result_elec->begin(); iter_5 != result_elec->end(); ++iter_5) {
                for ( pat::MuonCollection::const_iterator iter_6 = result_muon->begin(); iter_6 < result_muon->end() ; ++iter_6) {
                    std::cout << "(" << (iter_5-result_elec->begin())  << "," << (iter_6-result_muon->begin()) << ")" << std::endl;
                    for (const pat::Tau &tau2 : *result_tau ) {
                        runConfig.nbOfTriplets_ += 1;
                               //  Mise en forme des valeurs pour calculer avec MEM 
                        MEM_nplet nplet(*iter_5, *iter_6, tau2); // 
/*                        std::cout << "VVVVVVVVVVVVVVV - AVANT - VVVVVVVVVVVVVVV" << std::endl;
                        std::cout << "Lep1_4P (" << iter_5->px() << ", " << iter_5->py() << ", " << iter_5->pz() << ", " << iter_5->energy() << ") " << std::endl;
                        std::cout << "Lep1_4P (" << iter_6->px() << ", " << iter_6->py() << ", " << iter_6->pz() << ", " << iter_6->energy() << ") " << std::endl;
                        std::cout << "Had1_4P (" << tau2.px()    << ", " << tau2.py()    << ", " << tau2.pz()    << ", " << tau2.energy() << ")" <<std::endl;
                        std::cout << "VVVVVVVVVVVVVVV - AVANT - VVVVVVVVVVVVVVV" << std::endl;*/
                        vector<pat::Jet> cleanedJets = npletFilter<pat::Electron, pat::Muon>( *iter_5, *iter_6, tau2, jets );
                        JetFilling<pat::Electron, pat::Muon>( *iter_5, *iter_6, tau2, cleanedJets, nplet, "em" );
                        std::cout << "Jets_4P size = " << nplet.Jets_4P.size() << std::endl;
/*                        std::cout << "VVVVVVVVVVVVVVV - APRES - VVVVVVVVVVVVVVV" << std::endl;
                        std::cout << "Lep1_4P (" << nplet.Lep1_4P.Px() << ", " << nplet.Lep1_4P.Py() << ", " << nplet.Lep1_4P.Pz() << ", " << nplet.Lep1_4P.E() << ") " << std::endl;
                        std::cout << "Lep1_4P (" << nplet.Lep2_4P.Px() << ", " << nplet.Lep2_4P.Py() << ", " << nplet.Lep2_4P.Pz() << ", " << nplet.Lep2_4P.E() << ") " << std::endl;
                        std::cout << "Had1_4P (" << nplet.HadSys_4P.Px() << ", " << nplet.HadSys_4P.Py() << ", " << nplet.HadSys_4P.Pz() << ", " << nplet.HadSys_4P.E() << ")" <<std::endl;
                        std::cout << "VVVVVVVVVVVVVVV - APRES - VVVVVVVVVVVVVVV" << std::endl;*/
                        runConfig.nbOfNplets_ += 1;
                        std::cout << "classical : nbOfNplets in em : " << runConfig.nbOfNplets_ << std::endl;
                        std::cout << "integration_type : " << nplet.integration_type   << std::endl;
                        std::cout << "lep1 type        : " << nplet.lep1_type << std::endl;
                        std::cout << "lep2 type        : " << nplet.lep2_type << std::endl;
                        std::cout << "tau decayMode    : " << nplet.decayMode << std::endl;
                        
                        /* Mets */
                        if ( mets->size() > 1 ) {
                            std::cout << "!!!!! WARNING Met size : " << mets->size() << std::endl;
                        }
                        else {
                            std::cout << "Met size : " << mets->size() << std::endl ;
                        }
                                
                                // initialize MET 
                        /*np.EventID = iEvent.id().event();
                        np.LumiID = iEvent.id().luminosityBlock();
                        np.RunID = iEvent.id().run();

                        TMatrixD covMET(2, 2); // ne peut pas etre initialise dans le .h
                        Handle<double> significanceHandle;
                        iEvent.getByToken (theSigTag, significanceHandle);*/
                    const pat::MET& srcMET = (*mets)[0];
                    Handle<math::Error<2>::type> covHandle;
                    iEvent.getByToken (theCovTag, covHandle);/**/
                    nplet.recoMETCov[0] = (*covHandle)(0,0);
                    nplet.recoMETCov[1] = (*covHandle)(1,0);
                    nplet.recoMETCov[2] = nplet.recoMETCov[1]; // (1,0) is the only one saved
                    nplet.recoMETCov[3] = (*covHandle)(1,1);/**/
/*                    std::cout << "Met cov matrix : " << srcMET.getSignificanceMatrix() << std::endl;
                    std::cout << "Met cov matrix : " << (*covHandle) << std::endl;
                    std::cout << "Met cov matrix : " << nplet.recoMETCov[0] << std::endl;*/
                    nplet.covarMET_display(); // to be removed
                    
                    std::cout << "Met pfMET     : " << srcMET.et() << std::endl;
                    std::cout << "Met pfMET phi : " << srcMET.phi() << std::endl;/**/
                    nplet.fill_recoMET_4P(srcMET.et(), srcMET.phi());
                    std::cout << "recoMET_4P loop (" << nplet.recoMET_4P.Pt() << ", " << nplet.recoMET_4P.Eta() << ", " << nplet.recoMET_4P.Phi() << ", " << nplet.recoMET_4P.M() << ") " << std::endl;
                        
                        // Appel MEM 
                        // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 

                    //std::cout << "one process, doing job" << std::endl;
                    oneMPIProcess(nplet); // EventReader<Run1EventData_t> &eventReader, IntegralsOutputs<T> *integralsOutput
                                    
                    //std::cout << "one process, ending job" << std::endl;
                        // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
                                
                                // Mise en forme des valeurs avec IntegralResult 
                                // Mise en forme des valeurs avec MEMResult                               
                                
                    }
                }
            }
            
        }
        
        std::cout << "That is all for Event : " << iEvent.id() << "\n" << std::endl ;

        // and save the vectors
        iEvent.put(result_elec);
        iEvent.put(result_muon);
        iEvent.put(result_met);
        iEvent.put(result_tau);
        iEvent.put(result_jet);
        
        iEvent.put(MEMResultColl);

//    }
/*    else { // processID_ != 0
        // 
        //  Workers
        //
        MPI_Send(greeting, strlen(greeting)+1, MPI_BYTE, 0,1,MPI_COMM_WORLD);

    }*/

}

void MEMProducer::oneMPIProcess(MEM_nplet &np) // EventReader<Run1EventData_t> &eventReader, IntegralsOutputs<T> *integralsOutput
{

    std::cout << "Je suis dans MEMProducer::oneMPIProcess" << std::endl;
    std::cout << "np recoMET_4P.Pt()   : " << np.recoMET_4P.Pt() << std::endl;
    std::cout << "np recoMET_4P.Eta()  : " << np.recoMET_4P.Eta() << std::endl;
    std::cout << "np recoMET_4P.Phi()  : " << np.recoMET_4P.Phi() << std::endl;
    std::cout << "np recoMET_4P.M()    : " << np.recoMET_4P.M() << std::endl;
    std::cout << "np recoMET_4P.Pt()   : " << np.eventList[0].evRecoMET4P_[0] << std::endl;
    std::cout << "np recoMET_4P.Eta()  : " << np.eventList[0].evRecoMET4P_[1] << std::endl;
    std::cout << "np recoMET_4P.Phi()  : " << np.eventList[0].evRecoMET4P_[2] << std::endl;
    std::cout << "np recoMET_4P.M()    : " << np.eventList[0].evRecoMET4P_[3] << std::endl;
     
/*    np.integration.MPIInfo_ = 0;
    np.integration.integration_type_ = np.integration_type;
    np.integration.evLep1_4P_[0] = np.Lep1_4P.Px(); np.integration.evLep1_4P_[1] = np.Lep1_4P.Py(); 
    np.integration.evLep1_4P_[2] = np.Lep1_4P.Pz(); np.integration.evLep1_4P_[3] = np.Lep1_4P.E();
    np.integration.lepton1_Type_ = np.lep1_type;
    np.integration.evLep2_4P_[0] = np.Lep2_4P.Px(); np.integration.evLep2_4P_[1] = np.Lep2_4P.Py(); 
    np.integration.evLep2_4P_[2] = np.Lep2_4P.Pz(); np.integration.evLep2_4P_[3] = np.Lep2_4P.E();
    np.integration.lepton2_Type_ = np.lep1_type;
    np.integration.evHadSys_Tau_4P_[0] = np.HadSys_4P.Px(); np.integration.evHadSys_Tau_4P_[1] = np.HadSys_4P.Py();
    np.integration.evHadSys_Tau_4P_[2] = np.HadSys_4P.Pz(); np.integration.evHadSys_Tau_4P_[3] = np.HadSys_4P.E();
    np.integration.HadtauDecayMode_ =  np.decayMode;
    
    np.integration.evBJet1_4P_[0] = np.BJet1_4P.Px(); np.integration.evBJet1_4P_[1] = np.BJet1_4P.Py();
    np.integration.evBJet1_4P_[2] = np.BJet1_4P.Pz(); np.integration.evBJet1_4P_[3] = np.BJet1_4P.E();
    np.integration.evBJet2_4P_[0] = np.BJet2_4P.Px(); np.integration.evBJet2_4P_[1] = np.BJet2_4P.Py();
    np.integration.evBJet2_4P_[2] = np.BJet2_4P.Pz(); np.integration.evBJet2_4P_[3] = np.BJet2_4P.E();*/
      //
    np.eventList[0].MPIInfo_ = 0;
    np.eventList[0].integration_type_ = np.integration_type;
    np.eventList[0].evLep1_4P_[0] = np.Lep1_4P.Px(); np.eventList[0].evLep1_4P_[1] = np.Lep1_4P.Py(); 
    np.eventList[0].evLep1_4P_[2] = np.Lep1_4P.Pz(); np.eventList[0].evLep1_4P_[3] = np.Lep1_4P.E();
    np.eventList[0].lepton1_Type_ = np.lep1_type;
    np.eventList[0].evLep2_4P_[0] = np.Lep2_4P.Px(); np.eventList[0].evLep2_4P_[1] = np.Lep2_4P.Py(); 
    np.eventList[0].evLep2_4P_[2] = np.Lep2_4P.Pz(); np.eventList[0].evLep2_4P_[3] = np.Lep2_4P.E();
    np.eventList[0].lepton2_Type_ = np.lep1_type;
    np.eventList[0].evHadSys_Tau_4P_[0] = np.HadSys_4P.Px(); np.eventList[0].evHadSys_Tau_4P_[1] = np.HadSys_4P.Py();
    np.eventList[0].evHadSys_Tau_4P_[2] = np.HadSys_4P.Pz(); np.eventList[0].evHadSys_Tau_4P_[3] = np.HadSys_4P.E();
    np.eventList[0].HadtauDecayMode_ =  np.decayMode;
    
    np.eventList[0].evBJet1_4P_[0] = np.BJet1_4P.Px(); np.eventList[0].evBJet1_4P_[1] = np.BJet1_4P.Py();
    np.eventList[0].evBJet1_4P_[2] = np.BJet1_4P.Pz(); np.eventList[0].evBJet1_4P_[3] = np.BJet1_4P.E();
    np.eventList[0].evBJet2_4P_[0] = np.BJet2_4P.Px(); np.eventList[0].evBJet2_4P_[1] = np.BJet2_4P.Py();
    np.eventList[0].evBJet2_4P_[2] = np.BJet2_4P.Pz(); np.eventList[0].evBJet2_4P_[3] = np.BJet2_4P.E();
      //
    //np.integration.n_lightJets_ = min(10,int(np.Jets_4P.size()));
    np.eventList[0].n_lightJets_ = min(10,int(np.Jets_4P.size()));
//    for( int i=0; i<np.integration.n_lightJets_; i++){
    for( int i=0; i<np.eventList[0].n_lightJets_; i++){
/*        np.integration.evJets_4P_[i][0] = np.Jets_4P[i].Px(); np.integration.evJets_4P_[i][1] = np.Jets_4P[i].Py();
        np.integration.evJets_4P_[i][2] = np.Jets_4P[i].Pz(); np.integration.evJets_4P_[i][3] = np.Jets_4P[i].E();*/
        np.eventList[0].evJets_4P_[i][0] = np.Jets_4P[i].Px(); np.eventList[0].evJets_4P_[i][1] = np.Jets_4P[i].Py();
        np.eventList[0].evJets_4P_[i][2] = np.Jets_4P[i].Pz(); np.eventList[0].evJets_4P_[i][3] = np.Jets_4P[i].E();
    }
    if(np.Jets_4P.size()<10){	
        for(unsigned int i=np.Jets_4P.size(); i<10; i++){
/*            np.integration.evJets_4P_[i][0] = 0; np.integration.evJets_4P_[i][1] = 0;
            np.integration.evJets_4P_[i][2] = 0; np.integration.evJets_4P_[i][3] = 0;*/
            np.eventList[0].evJets_4P_[i][0] = 0; np.eventList[0].evJets_4P_[i][1] = 0;
            np.eventList[0].evJets_4P_[i][2] = 0; np.eventList[0].evJets_4P_[i][3] = 0;
        }
    }/**/
      
    if(np.integration_type == 0){
/*        np.integration.evJet1_4P_[0] = np.Jet1_4P.Px(); np.integration.evJet1_4P_[1] = np.Jet1_4P.Py();
        np.integration.evJet1_4P_[2] = np.Jet1_4P.Pz(); np.integration.evJet1_4P_[3] = np.Jet1_4P.E();
        np.integration.evJet2_4P_[0] = np.Jet2_4P.Px(); np.integration.evJet2_4P_[1] = np.Jet2_4P.Py();
        np.integration.evJet2_4P_[2] = np.Jet2_4P.Pz(); np.integration.evJet2_4P_[3] = np.Jet2_4P.E();*/
        np.eventList[0].evJet1_4P_[0] = np.Jet1_4P.Px(); np.eventList[0].evJet1_4P_[1] = np.Jet1_4P.Py();
        np.eventList[0].evJet1_4P_[2] = np.Jet1_4P.Pz(); np.eventList[0].evJet1_4P_[3] = np.Jet1_4P.E();
        np.eventList[0].evJet2_4P_[0] = np.Jet2_4P.Px(); np.eventList[0].evJet2_4P_[1] = np.Jet2_4P.Py();
        np.eventList[0].evJet2_4P_[2] = np.Jet2_4P.Pz(); np.eventList[0].evJet2_4P_[3] = np.Jet2_4P.E();
    }
    else if(np.integration_type == 1){
/*        np.integration.evJet1_4P_[0] = 0; np.integration.evJet1_4P_[1] = 0;
        np.integration.evJet1_4P_[2] = 0; np.integration.evJet1_4P_[3] = 0;
        np.integration.evJet2_4P_[0] = 0; np.integration.evJet2_4P_[1] = 0;
        np.integration.evJet2_4P_[2] = 0; np.integration.evJet2_4P_[3] = 0;	*/
        np.eventList[0].evJet1_4P_[0] = 0; np.eventList[0].evJet1_4P_[1] = 0;
        np.eventList[0].evJet1_4P_[2] = 0; np.eventList[0].evJet1_4P_[3] = 0;
        np.eventList[0].evJet2_4P_[0] = 0; np.eventList[0].evJet2_4P_[1] = 0;
        np.eventList[0].evJet2_4P_[2] = 0; np.eventList[0].evJet2_4P_[3] = 0;	
    }/**/
      //
/*    np.integration.evRecoMET4P_[0] = np.recoMET_4P.Px(); np.integration.evRecoMET4P_[1] = np.recoMET_4P.Py();
    np.integration.evRecoMET4P_[2] = np.recoMET_4P.Pz(); np.integration.evRecoMET4P_[3] = np.recoMET_4P.E();
      //
    np.integration.evV_[0] = np.recoMETCov[0]; np.integration.evV_[1] = np.recoMETCov[1]; 
    np.integration.evV_[2] = np.recoMETCov[2]; np.integration.evV_[3] = np.recoMETCov[3];     */
    
    // TEMPORATY OUTPUT
    std::cout << "np recoMET_4P.Pt()   : " << np.recoMET_4P.Pt() << std::endl;
    std::cout << "np recoMET_4P.Eta()  : " << np.recoMET_4P.Eta() << std::endl;
    std::cout << "np recoMET_4P.Phi()  : " << np.recoMET_4P.Phi() << std::endl;
    std::cout << "np recoMET_4P.M()    : " << np.recoMET_4P.M() << std::endl;
    // TEMPORATY OUTPUT
    np.eventList[0].evRecoMET4P_[0] = np.recoMET_4P.Px(); np.eventList[0].evRecoMET4P_[1] = np.recoMET_4P.Py(); // WARNING here Px, Py, Pz, E and fill_recoMET_4P use SetPtEtaPhiM !!
    np.eventList[0].evRecoMET4P_[2] = np.recoMET_4P.Pz(); np.eventList[0].evRecoMET4P_[3] = np.recoMET_4P.E();
    // TEMPORATY OUTPUT
    std::cout << "np recoMET_4P.Pt()   : " << np.eventList[0].evRecoMET4P_[0] << std::endl;
    std::cout << "np recoMET_4P.Eta()  : " << np.eventList[0].evRecoMET4P_[1] << std::endl;
    std::cout << "np recoMET_4P.Phi()  : " << np.eventList[0].evRecoMET4P_[2] << std::endl;
    std::cout << "np recoMET_4P.M()    : " << np.eventList[0].evRecoMET4P_[3] << std::endl;
    // TEMPORATY OUTPUT
    //
    np.eventList[0].evV_[0] = np.recoMETCov[0]; np.eventList[0].evV_[1] = np.recoMETCov[1]; 
    np.eventList[0].evV_[2] = np.recoMETCov[2]; np.eventList[0].evV_[3] = np.recoMETCov[3];     /**/
    
    std::cout << "VVVVVVVVVVVVVVV - classical nplet - VVVVVVVVVVVVVVV" << std::endl;
    std::cout << "Lep1_4P (" << np.Lep1_4P.Px()   << ", " << np.Lep1_4P.Py()   << ", " << np.Lep1_4P.Pz()   << ", " << np.Lep1_4P.E()   << ") " << std::endl;
    std::cout << "Lep1_4P (" << np.Lep2_4P.Px()   << ", " << np.Lep2_4P.Py()   << ", " << np.Lep2_4P.Pz()   << ", " << np.Lep2_4P.E()   << ") " << std::endl;
    std::cout << "Had1_4P (" << np.HadSys_4P.Px() << ", " << np.HadSys_4P.Py() << ", " << np.HadSys_4P.Pz() << ", " << np.HadSys_4P.E() << ") " << std::endl;
/*    std::cout << "VVVVVVVVVVVVVVV - integration - VVVVVVVVVVVVVVV" << std::endl;
    std::cout << "Lep1_4P (" << np.integration.evLep1_4P_[0]   << ", " << np.integration.evLep1_4P_[1]   << ", " << np.integration.evLep1_4P_[2]   << ", " << np.integration.evLep1_4P_[3]   << ") " << std::endl;
    std::cout << "Lep1_4P (" << np.integration.evLep2_4P_[0]   << ", " << np.integration.evLep2_4P_[1]   << ", " << np.integration.evLep2_4P_[2]   << ", " << np.integration.evLep2_4P_[3]   << ") " << std::endl;
    std::cout << "Had1_4P (" << np.integration.evHadSys_Tau_4P_[0] << ", " << np.integration.evHadSys_Tau_4P_[1] << ", " << np.integration.evHadSys_Tau_4P_[2] << ", " << np.integration.evHadSys_Tau_4P_[3] << ") " << std::endl;*/
    std::cout << "VVVVVVVVVVVVVVV - eventList - VVVVVVVVVVVVVVV" << std::endl;/**/
    std::cout << "Lep1_4P (" << np.eventList[0].evLep1_4P_[0]   << ", " << np.eventList[0].evLep1_4P_[1]   << ", " << np.eventList[0].evLep1_4P_[2]   << ", " << np.eventList[0].evLep1_4P_[3]   << ") " << std::endl;
    std::cout << "Lep1_4P (" << np.eventList[0].evLep2_4P_[0]   << ", " << np.eventList[0].evLep2_4P_[1]   << ", " << np.eventList[0].evLep2_4P_[2]   << ", " << np.eventList[0].evLep2_4P_[3]   << ") " << std::endl;
    std::cout << "Had1_4P (" << np.eventList[0].evHadSys_Tau_4P_[0] << ", " << np.eventList[0].evHadSys_Tau_4P_[1] << ", " << np.eventList[0].evHadSys_Tau_4P_[2] << ", " << np.eventList[0].evHadSys_Tau_4P_[3] << ") " << std::endl;/**/

      // MPI process source
      //integration.MPIInfo_ = 0;
    scheduler->runNodeScheduler ( np.eventList, 1 ); // (eventList, 1)
    std::cout << "======" << std::endl;
    //std::cout << "nbrOfPermut_ : "<< np.integration.nbrOfPermut_ << std::endl;
    std::cout << "nbrOfPermut_ : "<< np.eventList[0].nbrOfPermut_<< std::endl;
    for( int perm = 0; perm < np.eventList[0].nbrOfPermut_; perm++ ){
        std::cout << "integralttH_[" << perm << "] : " << np.eventList[0].integralttH_[perm] << std::endl;
    }
    std::cout << "======" << std::endl;

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------

void
MEMProducer::beginStream(edm::StreamID)
{
  std::cout 
//    << "cout:STREAM nbrOfProcess = " << nbrOfProcess_ 
//    << "  processID = " << processID_ 
//    << "  processor_name = " << processor_name 
    << std::endl;
}
/**/

// ------------ method called once each stream after processing all runs, lumis and events  ------------

void
MEMProducer::endStream() {
  std::cout << "ttl nb of Electrons : " << runConfig.ttlnbOfElectrons_ << std::endl;
  std::cout << "ttl nb of Muons     : " << runConfig.ttlnbOfMuons_ << std::endl;
  std::cout << "ttl nb of Taus      : " << runConfig.ttlnbOfTaus_ << std::endl;
  std::cout << "nbOfTriplets        : " << runConfig.nbOfTriplets_ << std::endl;
  std::cout << "nbOfNplets          : " << runConfig.nbOfNplets_ << std::endl;
  std::cout << "-----" << std::endl;

  std::cout 
//    << "cout:END STREAM .... processID=" << processID_
    << "\nJe sors" 
    << std::endl;
}
/**/
 
// ------------ method called when starting to processes a luminosity block  ------------

void
MEMProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
/**/
 
// ------------ method called when ending the processing of a luminosity block  ------------

void
MEMProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
/**/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MEMProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

