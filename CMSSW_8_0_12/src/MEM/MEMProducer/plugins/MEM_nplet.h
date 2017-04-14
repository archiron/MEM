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
//         Created:  Thu, 14 Apr 2016 14:58:00 GMT
//
//          ttH version
//
#ifndef MEM_MEMProducer_MEM_nplet_h 
#define MEM_MEMProducer_MEM_nplet_h


// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>
#include <iostream>
#include <tuple>
#include <utility>
#include <typeinfo>
#include "TMath.h"
#include "TLorentzVector.h"
#include <TMatrixD.h>
#include <Math/SMatrix.h>
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "MEM/MEMAlgo/interface/ThreadScheduler.h"

//
// structures needed to retrieve particles characteristics
// ID : ID in the objects loop ( for (const pat::Muon &muon : *muons ){} )
// E, pt, eta, phi
//

//
// class declaration
//

class MEM_nplet {
    public:
        template <typename T, typename U>
        MEM_nplet(const T& lep1, const U& lep2, const pat::Tau& tau) {
            //std::cout << "fonction template MEM_nplet" << std::endl;

            Lep1_4P.SetPx( lep1.px() );
            Lep1_4P.SetPy( lep1.py() );
            Lep1_4P.SetPz( lep1.pz() );
            Lep1_4P.SetE( lep1.energy() );  
            Lep2_4P.SetPx( lep2.px() );
            Lep2_4P.SetPy( lep2.py() );
            Lep2_4P.SetPz( lep2.pz() );
            Lep2_4P.SetE( lep2.energy() );  
            HadSys_4P.SetPx( tau.px() );
            HadSys_4P.SetPy( tau.py() );
            HadSys_4P.SetPz( tau.pz() );
            HadSys_4P.SetE( tau.energy() );
            
        }
        ~MEM_nplet();
        
        double mee(const pat::Jet J1, const pat::Jet J2);
        bool diff_mass(const pat::Jet J1, const pat::Jet J2);
        bool diff_mass_2(double m1, double m2);
        void fill_BJet1_4P(const pat::Jet J1); 
        void fill_BJet2_4P(const pat::Jet J2); 
        void fill_Jet1_4P(const pat::Jet J1); 
        void fill_Jet2_4P(const pat::Jet J2); 
        void fill_Jet2_4P_0(); 
        void fill_recoMET_4P(float met, float phi);
        void fillEvent();
        void fillEvent_num();
        TLorentzVector fill_temp(const pat::Jet J1); 

        void covarMET_display();
        void nplet_display();
        void eventList_display();
      
        /* source : https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_5_0/doc/html/d1/df0/MultiEventFilter_8cc_source.html */
        unsigned long long EventID; // from iEvent.id().event()
        unsigned int RunID; // from iEvent.id().run()
        unsigned int LumiID; // from iEvent.id().lumi()
        int integration_type = 1; // default Delta m > 20.
        int lep1_type = 0, lep2_type = 0;
        int decayMode = -1;
        float d_eta;
        
        TLorentzVector Lep1_4P, Lep2_4P;
        TLorentzVector HadSys_4P;
        TLorentzVector BJet1_4P, BJet2_4P;
        TLorentzVector Jet1_4P, Jet2_4P;
        std::vector<TLorentzVector> Jets_4P;        
        TLorentzVector recoMET_4P;
        double recoMETCov[4];
        eventList_t eventList;

    private:

      // ----------member data ---------------------------
    
};


#endif
