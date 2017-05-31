// -*- C++ -*-
//
// Package:    MEM/MEMProducer
// Class:      MEM_nplet
// 
/**\class MEM_nplet MEM_nplet.cc MEM/MEMProducer/plugins/MEM_nplet.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Arnaud Chiron
//         Created:  Mon, 1a Apr 2016 14:58:00 GMT
//
//

// system include files
#include "MEM/MEMProducer/plugins/MEM_nplet.h"
#include "MEM/MEMAlgo/interface/Constants.h"

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include <vector>
#include <tuple>
#include <utility>
#include "TMath.h"
#include <typeinfo>
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

using namespace reco;
using namespace std;

MEM_nplet::~MEM_nplet()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

  //std::cout << "Destructor nplet" << std::endl;

}

//
// member functions
//

void MEM_nplet::fill_BJet1_4P(const pat::Jet J1)
{
    BJet1_4P.SetPx( J1.px() );
    BJet1_4P.SetPy( J1.py() );
    BJet1_4P.SetPz( J1.pz() );
    BJet1_4P.SetE( J1.energy() );  
}

void MEM_nplet::fill_BJet2_4P(const pat::Jet J2)
{
    BJet2_4P.SetPx( J2.px() );
    BJet2_4P.SetPy( J2.py() );
    BJet2_4P.SetPz( J2.pz() );
    BJet2_4P.SetE( J2.energy() );  
}

void MEM_nplet::fill_Jet1_4P(const pat::Jet J1)
{
    Jet1_4P.SetPx( J1.px() );
    Jet1_4P.SetPy( J1.py() );
    Jet1_4P.SetPz( J1.pz() );
    Jet1_4P.SetE( J1.energy() );  
}

void MEM_nplet::fill_Jet2_4P(const pat::Jet J2)
{
    Jet2_4P.SetPx( J2.px() );
    Jet2_4P.SetPy( J2.py() );
    Jet2_4P.SetPz( J2.pz() );
    Jet2_4P.SetE( J2.energy() );  
}

void MEM_nplet::fill_Jet2_4P_0()
{
    Jet2_4P.SetPx( 0. );
    Jet2_4P.SetPy( 0. );
    Jet2_4P.SetPz( 0. );
    Jet2_4P.SetE( 0. );  
}

void MEM_nplet::fill_recoMET_4P(float met, float phi)
{
    //std::cout << "inside fill_recoMET_4P : " << met << " - " << phi << std::endl;//
    recoMET_4P.SetPtEtaPhiM( met, 0.,  phi , 0. );
    //std::cout << "inside fill_recoMET_4P : " << recoMET_4P.Pt() << ", " << recoMET_4P.Eta() << ", " << recoMET_4P.Phi() << ", " << recoMET_4P.M() << std::endl;
}

void MEM_nplet::fillEvent()
{
    ;
    eventList[0].MPIInfo_ = 0;
    eventList[0].integration_type_ = integration_type;
    eventList[0].evLep1_4P_[0] = Lep1_4P.Px(); eventList[0].evLep1_4P_[1] = Lep1_4P.Py(); 
    eventList[0].evLep1_4P_[2] = Lep1_4P.Pz(); eventList[0].evLep1_4P_[3] = Lep1_4P.E();
    eventList[0].lepton1_Type_ = lep1_type;
    eventList[0].evLep2_4P_[0] = Lep2_4P.Px(); eventList[0].evLep2_4P_[1] = Lep2_4P.Py(); 
    eventList[0].evLep2_4P_[2] = Lep2_4P.Pz(); eventList[0].evLep2_4P_[3] = Lep2_4P.E();
    eventList[0].lepton2_Type_ = lep1_type;
    eventList[0].evHadSys_Tau_4P_[0] = HadSys_4P.Px(); eventList[0].evHadSys_Tau_4P_[1] = HadSys_4P.Py();
    eventList[0].evHadSys_Tau_4P_[2] = HadSys_4P.Pz(); eventList[0].evHadSys_Tau_4P_[3] = HadSys_4P.E();
    eventList[0].HadtauDecayMode_ =  decayMode;
    
    eventList[0].evBJet1_4P_[0] = BJet1_4P.Px(); eventList[0].evBJet1_4P_[1] = BJet1_4P.Py();
    eventList[0].evBJet1_4P_[2] = BJet1_4P.Pz(); eventList[0].evBJet1_4P_[3] = BJet1_4P.E();
    eventList[0].evBJet2_4P_[0] = BJet2_4P.Px(); eventList[0].evBJet2_4P_[1] = BJet2_4P.Py();
    eventList[0].evBJet2_4P_[2] = BJet2_4P.Pz(); eventList[0].evBJet2_4P_[3] = BJet2_4P.E();
    
    eventList[0].n_lightJets_ = min(10,int(Jets_4P.size()));
    std::cout << "Jets_4P.size() : " << Jets_4P.size()            << std::endl; // TEMPORAIRE
    std::cout << "n_lightJets_      : " << eventList[0].n_lightJets_ << std::endl; // TEMPORAIRE

    for( int i=0; i<eventList[0].n_lightJets_; i++){
        eventList[0].evJets_4P_[i][0] = Jets_4P[i].Px(); eventList[0].evJets_4P_[i][1] = Jets_4P[i].Py();
        eventList[0].evJets_4P_[i][2] = Jets_4P[i].Pz(); eventList[0].evJets_4P_[i][3] = Jets_4P[i].E();
    }
    if(Jets_4P.size()<10){	
        for(unsigned int i=Jets_4P.size(); i<10; i++){
            eventList[0].evJets_4P_[i][0] = 0; eventList[0].evJets_4P_[i][1] = 0;
            eventList[0].evJets_4P_[i][2] = 0; eventList[0].evJets_4P_[i][3] = 0;
        }
    }

    if(integration_type == 0){
        eventList[0].evJet1_4P_[0] = Jet1_4P.Px(); eventList[0].evJet1_4P_[1] = Jet1_4P.Py();
        eventList[0].evJet1_4P_[2] = Jet1_4P.Pz(); eventList[0].evJet1_4P_[3] = Jet1_4P.E();
        eventList[0].evJet2_4P_[0] = Jet2_4P.Px(); eventList[0].evJet2_4P_[1] = Jet2_4P.Py();
        eventList[0].evJet2_4P_[2] = Jet2_4P.Pz(); eventList[0].evJet2_4P_[3] = Jet2_4P.E();
    }
    else if(integration_type == 1){
        eventList[0].evJet1_4P_[0] = 0; eventList[0].evJet1_4P_[1] = 0;
        eventList[0].evJet1_4P_[2] = 0; eventList[0].evJet1_4P_[3] = 0;
        eventList[0].evJet2_4P_[0] = 0; eventList[0].evJet2_4P_[1] = 0;
        eventList[0].evJet2_4P_[2] = 0; eventList[0].evJet2_4P_[3] = 0;	
    }/**/

    eventList[0].evRecoMET4P_[0] = recoMET_4P.Px(); eventList[0].evRecoMET4P_[1] = recoMET_4P.Py(); 
    eventList[0].evRecoMET4P_[2] = recoMET_4P.Pz(); eventList[0].evRecoMET4P_[3] = recoMET_4P.E();

    eventList[0].evV_[0] = recoMETCov[0]; eventList[0].evV_[1] = recoMETCov[1]; 
    eventList[0].evV_[2] = recoMETCov[2]; eventList[0].evV_[3] = recoMETCov[3];     /**/
    
}

void MEM_nplet::fillEvent_num()
{
    ;
    eventList[0].MPIInfo_ = 0;
    eventList[0].integration_type_ = integration_type;
    eventList[0].evLep1_4P_[0] = -4.61837; eventList[0].evLep1_4P_[1] = -68.0647; 
    eventList[0].evLep1_4P_[2] = -30.2518; eventList[0].evLep1_4P_[3] = 74.6279;
    eventList[0].lepton1_Type_ = -13;
    eventList[0].evLep2_4P_[0] = 33.7741; eventList[0].evLep2_4P_[1] = 34.1932; 
    eventList[0].evLep2_4P_[2] = 33.7678; eventList[0].evLep2_4P_[3] = 58.7379;
    eventList[0].lepton2_Type_ = -13;
    eventList[0].evHadSys_Tau_4P_[0] = -28.1687; eventList[0].evHadSys_Tau_4P_[1] = -17.799;
    eventList[0].evHadSys_Tau_4P_[2] = -81.1725; eventList[0].evHadSys_Tau_4P_[3] = 87.7606;
    eventList[0].HadtauDecayMode_ =  1;
    
    eventList[0].evBJet1_4P_[0] = -10.7871; eventList[0].evBJet1_4P_[1] = -29.1924;
    eventList[0].evBJet1_4P_[2] = 31.5088; eventList[0].evBJet1_4P_[3] = 44.6711;
    eventList[0].evBJet2_4P_[0] = -58.8208; eventList[0].evBJet2_4P_[1] = 26.2084;
    eventList[0].evBJet2_4P_[2] = -167.28; eventList[0].evBJet2_4P_[3] = 179.449;
    
    eventList[0].n_lightJets_ = min(10,int(Jets_4P.size()));
    //std::cout << "Jets_4P.size() : " << Jets_4P.size()            << std::endl;
    //std::cout << "n_lightJets_      : " << eventList[0].n_lightJets_ << std::endl;

    for( int i=0; i<eventList[0].n_lightJets_; i++){
        eventList[0].evJets_4P_[i][0] = Jets_4P[i].Px(); eventList[0].evJets_4P_[i][1] = Jets_4P[i].Py();
        eventList[0].evJets_4P_[i][2] = Jets_4P[i].Pz(); eventList[0].evJets_4P_[i][3] = Jets_4P[i].E();
    }
    
    if(Jets_4P.size()<10){	
        eventList[0].evJets_4P_[0][0] = 25.6494; eventList[0].evJets_4P_[0][1] = 37.7786; // TEMPORAIRE
        eventList[0].evJets_4P_[0][2] = -54.7096; eventList[0].evJets_4P_[0][3] = 71.7591;
        eventList[0].evJets_4P_[1][0] = 41.7213; eventList[0].evJets_4P_[1][1] = 0.0626074;
        eventList[0].evJets_4P_[1][2] = -98.0699; eventList[0].evJets_4P_[1][3] = 106.95;
        eventList[0].evJets_4P_[2][0] = -37.962; eventList[0].evJets_4P_[2][1] = 11.7409;
        eventList[0].evJets_4P_[2][2] = -22.0981; eventList[0].evJets_4P_[2][3] = 46.1146;
    }

    if(integration_type == 0){
        eventList[0].evJet1_4P_[0] = 25.6494; eventList[0].evJet1_4P_[1] = 37.7786;
        eventList[0].evJet1_4P_[2] = -54.7096; eventList[0].evJet1_4P_[3] = 71.7591;
        eventList[0].evJet2_4P_[0] = -37.962; eventList[0].evJet2_4P_[1] = 11.7409;
        eventList[0].evJet2_4P_[2] = -22.0981; eventList[0].evJet2_4P_[3] = 46.1146;
        
    }
    else if(integration_type == 1){
        eventList[0].evJet1_4P_[0] = 0; eventList[0].evJet1_4P_[1] = 0;
        eventList[0].evJet1_4P_[2] = 0; eventList[0].evJet1_4P_[3] = 0;
        eventList[0].evJet2_4P_[0] = 0; eventList[0].evJet2_4P_[1] = 0;
        eventList[0].evJet2_4P_[2] = 0; eventList[0].evJet2_4P_[3] = 0;	
    }/**/
    
    eventList[0].evRecoMET4P_[0] = recoMET_4P.Px(); eventList[0].evRecoMET4P_[1] = recoMET_4P.Py(); 
    eventList[0].evRecoMET4P_[2] = recoMET_4P.Pz(); eventList[0].evRecoMET4P_[3] = recoMET_4P.E();

    eventList[0].evV_[0] = recoMETCov[0]; eventList[0].evV_[1] = recoMETCov[1]; 
    eventList[0].evV_[2] = recoMETCov[2]; eventList[0].evV_[3] = recoMETCov[3];     /**/
    
}

TLorentzVector MEM_nplet::fill_temp(const pat::Jet J1)
{
    TLorentzVector temp1;
    temp1.SetPx( J1.px() );
    temp1.SetPy( J1.py() );
    temp1.SetPz( J1.pz() );
    temp1.SetE( J1.energy() );  
//    std::cout << "j1.pt = " << J1.pt() << " - J1.Px() = " << J1.px() << " - temp1.Px() = " << temp1.Px() << std::endl;
    return temp1;
}

double
MEM_nplet::mee(const pat::Jet J1, const pat::Jet J2)
{
    math::XYZTLorentzVector p12 = J1.p4()+J2.p4();
    //std::cout << "\t\t mmm (" << J1.pt() << "," << J2.pt() << ")" << std::endl; // OK, well transfered
    double mass = ( p12.Dot(p12) > 0. ? sqrt( p12.Dot(p12) ) : 0.);
    return mass;
}

bool
MEM_nplet::diff_mass(const pat::Jet J1, const pat::Jet J2)
{
    bool diff_m = false;
    double mass = mee(J1, J2);
    diff_m = ( fabs(mass - Physics::mW) < 20. ? true : false);
    //std::cout << "inside diff : mee= " << mass << std::endl;
    return diff_m;
}

bool
MEM_nplet::diff_mass_2(double m1, double m2)
{
    bool diff_m = false;
    diff_m = ( fabs(m1 - Physics::mW) > fabs(m2 - Physics::mW) ? true : false);
    return diff_m;
}

void
MEM_nplet::covarMET_display()
{
    std::cout << "covariant MET Matrix characteristics : " << std::endl;
    std::cout << "-----" << std::endl;
    std::cout << "covarMET[0] = " << recoMETCov[0] << std::endl ;
    std::cout << "covarMET[1] = " << recoMETCov[1] << std::endl ;
    std::cout << "covarMET[2] = " << recoMETCov[2] << std::endl ;
    std::cout << "covarMET[3] = " << recoMETCov[3] << std::endl ;
    std::cout << "-----" << std::endl;/**/

}

void
MEM_nplet::nplet_display()
{
    std::cout << "          nplet characteristics : " << std::endl;
    std::cout << "          Event ID : " << EventID << std::endl;
    std::cout << "          Lumi  ID : " << LumiID << std::endl;
    std::cout << "          Run   ID : " << RunID << std::endl;
    std::cout << "others characteristics : " << std::endl;
    covarMET_display();
    
}

void
MEM_nplet::eventList_display()
{
    if ( eventList[0].nbrOfPermut_ > 0 ) {
        std::cout << "======" << std::endl;
        std::cout << "nbrOfPermut_ : "<< eventList[0].nbrOfPermut_ << std::endl;
        /*for( int perm = 0; perm < eventList[0].nbrOfPermut_; perm++ ){
            std::cout << "eventList[0].integralttH_[" << perm << "] : " << eventList[0].integralttH_[perm] ; // << std::endl
            std::cout << " - include_perm_ttH_[" << perm << "] : " << eventList[0].include_perm_ttH_[perm] << std::endl;
        }*/
        std::cout << "======" << std::endl;
    }
    else {
        std::cout << "======" << std::endl;
        std::cout << "Pbm with nbrOfPermut_ : "<< eventList[0].nbrOfPermut_ << std::endl;
        std::cout << "======" << std::endl;
    }
}