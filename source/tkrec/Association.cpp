// Cimrman headers
#include "tkrec/Association.h"

// Cimrman headers
#include "tkrec/TrackerHit.h"

// ClassImp(tkrec::Association);

namespace tkrec {
  
  using namespace std;

  Association::Association(const ConstTrackerHitHdl & hit)
    : tracker_hit(hit)
  {
    return;
  }
  
} //  end of namespace tkrec

