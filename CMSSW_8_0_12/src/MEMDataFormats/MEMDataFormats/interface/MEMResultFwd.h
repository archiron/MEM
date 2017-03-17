#ifndef MEMDataFormats_MEMDataFormats__MEMResultFwd_h
#define MEMDataFormats_MEMDataFormats__MEMResultFwd_h

#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"

namespace reco {
  class MEMResult;
  /// collection of MEMResult objects
  typedef std::vector<MEMResult> MEMResultCollection;
  /// reference to an object in a collection of MEMResult objects
  typedef edm::Ref<MEMResultCollection> MEMResultRef;
  /// reference to a collection of MEMResult objects
  typedef edm::RefProd<MEMResultCollection> MEMResultRefProd;
  /// vector of objects in the same collection of MEMResult objects
  typedef edm::RefVector<MEMResultCollection> MEMResultRefVector;
  /// iterator over a vector of reference to MEMResult objects
  typedef MEMResultRefVector::iterator MEMResult_iterator;
}

#endif