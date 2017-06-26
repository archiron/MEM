#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/typelookup.h"

#include "MEM/MEMProducer/plugins/MEMProducer.h"
#include "MEM/MEMProducer/plugins/ExtractMETSignificance.h"
#include "MEM/MEMProducer/plugins/MEM_nplet.h"

#include "MEMDataFormats/MEMDataFormats/interface/IntegralResult.h"
#include "MEMDataFormats/MEMDataFormats/interface/MEMResult.h"
#include "MEMDataFormats/MEMDataFormats/interface/RunConfig.h"
#include "MEMDataFormats/MEMDataFormats/interface/PyRun2EventData_t.h"

DEFINE_FWK_MODULE(MEMProducer);
DEFINE_FWK_MODULE(ExtractMETSignificance);