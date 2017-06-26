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
//         Created:  Thu, 28 Jan 2016 16:43:54 GMT
//
//
#ifndef MEM_MEMProducer_MEMProducer_h 
#define MEM_MEMProducer_MEMProducer_h


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
#include <string>
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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "CondFormats/DataRecord/interface/JetResolutionRcd.h"
#include "FWCore/Framework/interface/eventsetuprecord_registration_macro.h"

#include "MEM/MEMProducer/plugins/MEM_nplet.h"
#include "MEMDataFormats/MEMDataFormats/interface/RunConfig.h"
//#include "MEMDataFormats/MEMDataFormats/interface/PyRun2EventData_t.h"
#include "MEM/MEMAlgo/interface/ThreadScheduler.h"

//
// structures needed to retrieve particles characteristics
// ID : ID in the objects loop ( for (const pat::Muon &muon : *muons ){} )
// E, pt, eta, phi
//

//
// class declaration
//

class MEMProducer : public edm::stream::EDProducer<> {
    public:
        explicit MEMProducer(const edm::ParameterSet&);
        ~MEMProducer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
            
    private:
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        virtual void beginStream(edm::StreamID) override;
        virtual void endStream() override;
      
        virtual void oneMPIProcess(MEM_nplet &np) ; // EventReader<Run1EventData_t> &eventReader, IntegralsOutputs<T> *integralsOutput

        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

        bool tauFilter(const pat::Tau& tauFiltre);
        bool muonFilter(const pat::Muon& muonFiltre);
        bool electronFilter(const pat::Electron& electronFiltre);
        bool jetFilter(const pat::Jet& jetFiltre);

        template <typename T, typename U>
        vector<pat::Jet> npletFilter(const T& lep1, const U& lep2, const pat::Tau& tau, edm::Handle<pat::JetCollection> jets) {
            vector<pat::Jet> cleaned_jets;
  
            /*std::cout << "\t inside template npletFilter" << std::endl;
            std::cout << "type id T = " << typeid(T).name() << std::endl;
            std::cout << "type id U = " << typeid(U).name() << std::endl;*/
            /*for (const pat::Jet &jet : *jets ) {
                std::cout << "npletFilter : " << jet.pt() << std::endl;;
            }*/

            math::XYZTLorentzVector lep1_P4 = lep1.p4();
            math::XYZTLorentzVector lep2_P4 = lep2.p4();
            math::XYZTLorentzVector tau_P4 = tau.p4();

            if ( !(deltaR(lep1_P4,lep2_P4)>0.3 && deltaR(lep1_P4,tau_P4)>0.3 && deltaR(lep2_P4,tau_P4)>0.3) )
                return cleaned_jets ; // null cleaned_jets_ordered

            for (const pat::Jet &jet : *jets ) {
                //std::cout << "npletFilter : " << jet.pt() << std::endl;;
                math::XYZTLorentzVector jet_P4 = jet.p4();
                if(jetFilter(jet) && deltaR(lep1_P4,jet_P4)>0.4 && deltaR(lep2_P4,jet_P4)>0.4 && deltaR(tau_P4,jet_P4)>0.4)
                    cleaned_jets.push_back(jet);
            }
          
            return cleaned_jets;

      }
      
        template <typename T, typename U>
        void JetFilling(const T& lep1, const U& lep2, const pat::Tau& tau, std::vector<pat::Jet> cleanedJets, MEM_nplet &nplet_template, std::string t_pairs) {
            //std::cout << "\t\t fonction template JetFilling for " << t_pairs << std::endl;
          
            //nplet_template.lep1_type = lep1.type;
            /*std::cout << "lep1 type     : " << nplet_template.lep1_type << std::endl;
            std::cout << "lep2 type     : " << nplet_template.lep2_type << std::endl;
            std::cout << "tau decayMode : " << nplet_template.decayMode << std::endl;
            std::cout << "tau decayMode : " << tau.decayMode() << std::endl;
            std::cout << "lep1 type     : " << lep1.pdgId() << std::endl;
            std::cout << "lep2 type     : " << lep2.pdgId() << std::endl;*/
            nplet_template.decayMode = tau.decayMode();
            nplet_template.lep1_type = lep1.pdgId();
            nplet_template.lep2_type = lep2.pdgId();
            /*std::cout << "lep1 type     : " << nplet_template.lep1_type << std::endl;
            std::cout << "lep2 type     : " << nplet_template.lep2_type << std::endl;
            std::cout << "tau decayMode : " << nplet_template.decayMode << std::endl;*/
            
            int nb_jets = cleanedJets.size();
            std::cout << "\t\t nb de jets cleanedJets in template = " << nb_jets << std::endl;
            if ( nb_jets >= 3 ) {
                int i_pos = 0;
                int i = 0;
                std::cout << "\t\t Jet : " << cleanedJets.begin()->pt() << " - " << cleanedJets.begin()->eta() << " - " << cleanedJets.begin()->phi() << std::endl;
                for ( pat::JetCollection::const_iterator i_jet_em = cleanedJets.begin() + 1; i_jet_em < cleanedJets.end(); ++i_jet_em ) {
                    i += 1;
                    if ( i_jet_em->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > cleanedJets[i_pos].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") ) { 
                        i_pos = i;
                    }
                    //std::cout << "\t\t Jet : " << i_jet_em->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << " - " << i_jet_em->pt() << " - " << i_jet_em->p4().px() << std::endl;
                    std::cout << "\t\t Jet : " << i_jet_em->pt() << " - " << i_jet_em->eta() << " - " << i_jet_em->phi() << std::endl;
                }
                std::cout << "\t\t Jet max : " << cleanedJets[i_pos].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;
                nplet_template.fill_BJet1_4P( cleanedJets[i_pos] ); // filling BJet1_4P
                cleanedJets.erase(cleanedJets.begin() + i_pos); // removing BJet1_4P
                //std::cout << "\t\t BJet1_4P (" << nplet_template.BJet1_4P.Px() << ", " << nplet_template.BJet1_4P.Py() << ", " << nplet_template.BJet1_4P.Pz() << ", " << nplet_template.BJet1_4P.E() << ")" << std::endl;
                std::cout << "\t\t BJet1_4P (" << nplet_template.BJet1_4P.Pt() << ", " << nplet_template.BJet1_4P.Eta() << ", " << nplet_template.BJet1_4P.Phi() << ")" << std::endl;
                i_pos = 0;
                i = 0;
                for ( pat::JetCollection::const_iterator i_jet_em = cleanedJets.begin() + 1; i_jet_em < cleanedJets.end(); ++i_jet_em ) {
                    i += 1;
                    if ( i_jet_em->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > cleanedJets[i_pos].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") ) { 
                        i_pos = i;
                    }
                    //std::cout << "\t\t Jet : " << i_jet_em->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << " - " << i_jet_em->pt() << " - " << i_jet_em->p4().px() << std::endl;
                }
                //std::cout << "\t\t Jet max : " << cleanedJets[i_pos].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;
                nplet_template.fill_BJet2_4P( cleanedJets[i_pos] ); // filling BJet2_4P
                cleanedJets.erase(cleanedJets.begin() + i_pos); // removing BJet2_4P
                //std::cout << "\t\t BJet2_4P (" << nplet_template.BJet2_4P.Px() << ", " << nplet_template.BJet2_4P.Py() << ", " << nplet_template.BJet2_4P.Pz() << ", " << nplet_template.BJet2_4P.E() << ")" << std::endl;
                std::cout << "\t\t BJet2_4P (" << nplet_template.BJet2_4P.Pt() << ", " << nplet_template.BJet2_4P.Eta() << ", " << nplet_template.BJet2_4P.Phi() << ")" << std::endl;
                
                if ( cleanedJets.size() > 2 ) { // at least 3 Jets remaining
                    std::cout << "\t\t\t cleanedJets size after BJets filling = " << cleanedJets.size() << std::endl;
                    // looking for the pair
                    double massJets_ref = nplet_template.mee(*cleanedJets.begin(), *(cleanedJets.begin() + 1));
                    auto &jet1 = *(cleanedJets.begin());
                    auto &jet2 = *(cleanedJets.begin() + 1);
                    for ( pat::JetCollection::const_iterator j1 = cleanedJets.begin(); j1 != cleanedJets.end() - 1; ++j1) {
                        for ( pat::JetCollection::const_iterator j2 = j1 + 1; j2 < cleanedJets.end() ; ++j2) {
                            if ( nplet_template.diff_mass(*j1, *j2) && nplet_template.diff_mass_2(massJets_ref, nplet_template.mee(*j1, *j2))) {
                                massJets_ref = nplet_template.mee(*j1, *j2);
                                jet1 = *j1;
                                jet2 = *j2;
                            }
                        }
                    }
                    std::cout << "min of Jets mass = " << massJets_ref << std::endl;
                    if ( nplet_template.diff_mass(jet1, jet2)) {
                        nplet_template.integration_type = 0;
                    }
                    if ( jet1.pt() > jet2.pt() ) {
                        nplet_template.fill_Jet1_4P( jet1 ); // filling Jet1_4P
                        nplet_template.fill_Jet2_4P( jet2 ); // filling Jet2_4P
                    }
                    else {
                        nplet_template.fill_Jet1_4P( jet2 ); // filling Jet4_4P
                        nplet_template.fill_Jet2_4P( jet1 ); // filling Jet2_4P
                    }
                    //cleanedJets.erase(cleanedJets.begin() + 1); // removing Jet2_4P // remove the removing, A.C. 2017/04/14
                    //cleanedJets.erase(cleanedJets.begin() + 0); // removing Jet1_4P // remove the removing, A.C. 2017/04/14
                                        
                    //std::cout << "inside integration_type : " << nplet_template.integration_type << std::endl;
                    
                }
                else if ( cleanedJets.size() == 2 ) { // only 2 Jets remaining
                    std::cout << "\t\t\t cleanedJets size after BJets filling = " << cleanedJets.size() << " Jets array remaining" << std::endl;
                    std::cout << "\t\t\t pt for first Jet  : " << cleanedJets[0].pt() << std::endl;
                    std::cout << "\t\t\t pt for second Jet : " << cleanedJets[1].pt() << std::endl;/**/
                    if ( cleanedJets[0].pt() > cleanedJets[1].pt() ) {
                        nplet_template.fill_Jet1_4P( cleanedJets[0] ); // filling Jet1_4P
                        nplet_template.fill_Jet2_4P( cleanedJets[1] ); // filling Jet2_4P
                    }
                    else {
                        nplet_template.fill_Jet1_4P( cleanedJets[1] ); // filling Jet4_4P
                        nplet_template.fill_Jet2_4P( cleanedJets[0] ); // filling Jet2_4P
                    }
                    //cleanedJets.erase(cleanedJets.begin() + 1); // removing Jet2_4P
                    //cleanedJets.erase(cleanedJets.begin() + 0); // removing Jet1_4P
                }
                else if ( cleanedJets.size() == 1 ) { // only 1 Jet remaining
                    std::cout << "\t\t cleanedJets size after BJets filling = " << cleanedJets.size() << " - keeping Jet2_4P == 0" << std::endl;
                    nplet_template.fill_Jet1_4P( cleanedJets[0] ); // filling Jet1_4P
                    nplet_template.fill_Jet2_4P_0(); // filling Jet2_4P == 0
                    //cleanedJets.erase(cleanedJets.begin() + 0); // removing Jet1_4P
                }
                else { // size = 0 ! if this, there is a BIG pbm !!
                    std::cout << "\t\t Houston, we have a pbm !!! = " << std::endl;
                }
                
                std::cout << "filling Jets_4P" << std::endl;
                for ( const pat::Jet &j1 : cleanedJets) {
                    std::cout << "j1 px= " << j1.px() << " - j1.pt() = " << j1.pt() << std::endl;
                    nplet_template.Jets_4P.push_back(nplet_template.fill_temp(j1));
                    std::cout << "fill_temp.Px() = " << nplet_template.fill_temp(j1).Px() << " - fill_temp.Pt() = " << nplet_template.fill_temp(j1).Pt() << std::endl;
                }
                // print the size of the cleanedJets array resulting
                //std::cout << "\t\t cleanedJets size after Jets filling = " << cleanedJets.size() << std::endl;
                std::cout << "\t\t Jets_4P size after Jets filling = " << nplet_template.Jets_4P.size() << std::endl;/**/
                std::cout << "end filling Jets_4P" << std::endl;

                
            }
            else if ( nb_jets > 0 ) {
                std::cout << "\t\t cleanedJets size = 0, no enough Jet for computation, must return nb_jets = " << nb_jets << std::endl;
            }
            else { // nb_jets = 0
                std::cout << "\t\t cleanedJets size = 0, no Jet for computation, must return nb_jets = 0" << std::endl;
            }
        }
      
        // ----------member data ---------------------------
        
        edm::EDGetTokenT<pat::ElectronCollection> electronToken_; 
        edm::EDGetTokenT<pat::MuonCollection> muonToken_;
        edm::EDGetTokenT<pat::TauCollection> tauToken_;
        edm::EDGetTokenT<pat::METCollection> metToken_;
        edm::EDGetTokenT<pat::JetCollection> jetToken_;/**/
        edm::EDGetTokenT<pat::METCollection> metTokenOne_;
        
        edm::EDGetTokenT<double> theSigTag;
        edm::EDGetTokenT<math::Error<2>::type> theCovTag;
      
};


#endif
