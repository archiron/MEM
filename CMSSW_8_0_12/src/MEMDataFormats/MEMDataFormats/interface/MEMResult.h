// -*- C++ -*-
//
// Package:    MEMDataFormats/MEMDataFormats
// Class:      MEMResult
// 
/**\class reco::MEMResult MEMResult.h MEMDataFormats/MEMDataFormats/plugins/MEMResult.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Arnaud Chiron
//         Created:  Wed, 18 May 2016 17:01:00 GMT
//
//
#ifndef MEMDataFormats_MEMDataFormats_MEMResult_h 
#define MEMDataFormats_MEMDataFormats_MEMResult_h


// system include files
#include <memory>

// user include files
#include "MEMDataFormats/MEMDataFormats/interface/IntegralResult.h"
#include "MEMDataFormats/MEMDataFormats/interface/MEMResultFwd.h"

#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"

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
#include <typeinfo>
#include "TMath.h"
#include <TMatrixD.h>
#include <Math/SMatrix.h>
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//
// structures needed to retrieve particles characteristics
//

//
// class declaration
//
namespace reco {

class MEMResult : public CompositePtrCandidate  {
//class MEMResult  {
    public:
        MEMResult();
        virtual ~MEMResult();

        // Set additional info
        void setFunction_MR();
        void setMemResult1(double v1, double v2);
        
        int subDet2 = 1;
        float valueOne = 0.;
        
        // Accessors
        int get_subDet2() const { return subDet2_ ; }
        float get_valueOne() const { return valueOne_ ; }

    private:

// ----------member data ---------------------------
        std::vector<IntegralResult> memResult1; //(0: signal; 1: bkg1, 2: bkg2, etc)
        int subDet2_ = 0;
        float valueOne_ = 0.;
};

}

#endif
