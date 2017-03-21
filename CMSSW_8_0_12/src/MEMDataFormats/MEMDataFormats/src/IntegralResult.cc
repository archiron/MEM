// -*- C++ -*-
//
// Package:    MEMDataFormats/MEMDataFormats
// Class:      IntegralResult
// 
/**\class IntegralResult IntegralResult.cc MEMDataFormats/MEMDataFormats/interface/IntegralResult.cc

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
#include "MEMDataFormats/MEMDataFormats/interface/IntegralResult.h"

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

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
IntegralResult::IntegralResult()
{
   //now do what ever other initialization is needed
   
  std::cout << "cout::IntegralResult Constructor " << std::endl;
  
}

IntegralResult::~IntegralResult()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

  std::cout << "IntegralResult Destructor " << std::endl;

}

//
// member functions
//

void IntegralResult::setFunction_IR(double v1, double v2)
{
  std::cout << "cout::IntegralResult setFunction call " << std::endl;
  chi2_ = v1;
  integral_ = v2;
}

// ------------ method called to produce the data  ------------


