// Cimrman headers
#include "datamodel/Association.h"
#include "datamodel/TrackerHit.h"

namespace cimrman::datamodel {
  
  Association::Association(const ConstTrackerHitHdl & hit)
    : tracker_hit(hit)
  {
    return;
  }
  
} //  end of namespace cimrman::datamodel

