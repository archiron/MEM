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
#ifndef MEMDataFormats_MEMDataFormats_IntegralResult_h 
#define MEMDataFormats_MEMDataFormats_IntegralResult_h


// system include files
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
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

# include "MEM/MEMAlgo/interface/Definitions.h"

//
// structures needed to retrieve particles characteristics
//

//
// class declaration
//

namespace reco {

class IntegralResult {
    public:
        IntegralResult();
        virtual ~IntegralResult();

        // Set additional info
        void setFunction_IR(double v1, double v2); // get MEM integral return
        
        double chi2 = 0.;
        
        // Accessors
        double get_chi2() const { return chi2_ ; }

    private:

        // ----------member data ---------------------------
        int processId_ = 0; // matrix element id: signal, bkg1, bkg2, bkg3, signalTESPlus, signalTESMinus, bkg1TESPlus 
        double integral_ = 0.; // result of the integral 
        double chi2_ = 0.; // chi2 of the integral 
        double error_ = 0.; // error on the integral 
        unsigned int niter_ = 0; // number of points used 
        int statusFlag_ = 0; // any other information - integration method 
  
        unsigned long long eventID_;
        bool flagSignal_; // char signal = (flagSignal) ? 'V' : 'D'; flagSignal is boolean
        float POutQuark1_[2];
        float POutQuark2_[2];
        double mTauTau2_;
        double PTauLep_;
        double PTauHad_;
        double cosThetaLep_;
        double phiTauLep_;
        double boost_; 
        double TFJet1_;
        double TFJet2_;
        double TFMET_;
        double TFLepTau_;
        double TFHadTau_;
        double Jac_;
        double wME_;
        double eval_;
        const char *errStr_;
        
};

}

#endif
