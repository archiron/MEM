// -*- C++ -*-
//
// Package:    MEMDataFormats/MEMDataFormats
// Class:      MEMResult
// 
/**\class reco::MEMResult MEMResult.cc MEMDataFormats/MEMDataFormats/plugins/MEMResult.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Arnaud Chiron
//         Created:  Wed, 18 May 2016 17:01:00 GMT
//
//

// system include files
#include "MEMDataFormats/MEMDataFormats/interface/MEMResult.h" 

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include <vector>
#include <iostream>
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
using namespace edm;
using namespace std;
using namespace ROOT::Math;
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MEMResult::MEMResult()
{
   //now do what ever other initialization is needed
   
  //std::cout << "cout::MEMResult Constructor " << std::endl;
  
}

MEMResult::~MEMResult()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

  //std::cout << "MEMResult Destructor. " << std::endl;

}

//
// member functions
//

void MEMResult::setFunction_MR()
{
  std::cout << "cout::MEMResult setFunction call " << std::endl;
  std::cout << "cout::MEMResult subDet2 = " << subDet2_ << std::endl;
  std::cout << "cout::MEMResult valueOne = " << valueOne_ << std::endl;
  valueOne_ = 1.2563;
  subDet2_ = 87;
  std::cout << "cout::MEMResult subDet2 = " << subDet2_ << std::endl;
  std::cout << "cout::MEMResult valueOne = " << valueOne_ << std::endl;

}

void MEMResult::setMemResult1(double v1, double v2)
{
  IntegralResult ir1;

  ir1.setFunction_IR(v1, v2);
  memResult1.push_back(ir1);
  
}
