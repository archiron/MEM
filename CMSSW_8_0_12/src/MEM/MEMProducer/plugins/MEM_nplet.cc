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

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//

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

