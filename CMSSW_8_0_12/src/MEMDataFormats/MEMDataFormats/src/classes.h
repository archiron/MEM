#include "MEMDataFormats/MEMDataFormats/interface/IntegralResult.h"
#include "MEMDataFormats/MEMDataFormats/interface/MEMResult.h"
#include "MEMDataFormats/MEMDataFormats/interface/MEMResultFwd.h"
#include "MEMDataFormats/MEMDataFormats/interface/RunConfig.h"
#include "MEMDataFormats/MEMDataFormats/interface/PyRun2EventData_t.h"

#include "MEM/MEMProducer/plugins/MEM_nplet.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace { struct dictionary {
    reco::MEMResult      mr;
    reco::IntegralResult ir;
    RunConfig rc;
    
    reco::MEMResultCollection m1;
    edm::Wrapper<reco::MEMResultCollection> w1;
    edm::Ref<reco::MEMResultCollection> r1;
    edm::RefProd<reco::MEMResultCollection> rp1;
    edm::RefVector<reco::MEMResultCollection> rv1;
    
};
}