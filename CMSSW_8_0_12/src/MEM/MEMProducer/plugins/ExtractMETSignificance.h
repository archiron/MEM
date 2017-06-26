/*
** class  : ExtractMETSignificance
** author : L. Cadamuro (LLR)
** date   : 4 November 2016
** brief  : takes pat::MET in input and produces a standalone significance and covariance collection from it
**          used as a replacement of the previous standalone significance producer
*/

#ifndef MEM_MEMProducer_ExtractMETSignificance_h 
#define MEM_MEMProducer_ExtractMETSignificance_h

#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/METReco/interface/CommonMETData.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

class ExtractMETSignificance : public edm::stream::EDProducer<> {
    public: 
        /// Constructor
        explicit ExtractMETSignificance(const edm::ParameterSet&);
        /// Destructor
        ~ExtractMETSignificance();  

    private:
//        virtual void beginJob(){};  
        virtual void produce(edm::Event&, const edm::EventSetup&);
//        virtual void endJob(){};
        virtual void beginStream(edm::StreamID) override;
        virtual void endStream() override;

        edm::EDGetTokenT<pat::METCollection> theMETTag;
        
        edm::EDGetTokenT<double> theSigTag;
        edm::EDGetTokenT<math::Error<2>::type> theCovTag;
        
};



#endif