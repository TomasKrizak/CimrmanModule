#ifndef FALAISE_CIMRMAN_VISU_H
#define FALAISE_CIMRMAN_VISU_H

// Standard headers
#include <functional>

// Bayeux:
#include <bayeux/datatools/bit_mask.h>

// Cimrman headers
#include "tkrec/Geometry.h"
#include "tkrec/Event.h"

// ROOT headers
#include "TCanvas.h"
#include "TColor.h"
#include "TEllipse.h"
#include "TH2D.h"
#include "TAttLine.h"
#include "TGLViewer.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TFile.h"
#include "TROOT.h"
#include "TPolyLine3D.h"
#include "TPolyLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPoint.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraph.h"
#include "TPolyMarker3D.h"

namespace tkrec {

  class Visu
  {
  public:
    
    Visu(const Geometry & geom_);

    bool has_event() const;
    
    void set_event(const Event & _event_);

    // vizualization section
    enum hit_flags {
      show_hit_invalid = datatools::bit_mask::bit00, ///< Show unused hits for reconstruction (yellow) invalid
      show_hit_used = datatools::bit_mask::bit01, ///< Show used hits for reconstruction (red)
      show_hit_associated = datatools::bit_mask::bit02, ///< Show associated hits to track (green)
      show_hit_novertical = datatools::bit_mask::bit03, ///< Differentiate hits with failed vertical position (teal - associatied if turned on, magenta - unassociated)
      show_hit_delayed = datatools::bit_mask::bit04, ///< Differentiate delayed hits (blue)
      show_hit_all = show_hit_invalid | show_hit_used | show_hit_associated | show_hit_novertical | show_hit_delayed
    };
    				 
    // tracking options:
    enum tracking_flags  {
      show_tracking_tracks = datatools::bit_mask::bit00,  ///< Show tracks // TODO? do I need this?
      show_tracking_hit_avalanches = datatools::bit_mask::bit01,  ///< Show hit avalanches
      show_tracking_trajectories = datatools::bit_mask::bit02, ///< Show trajectories
      show_tracking_all = show_tracking_tracks | show_tracking_hit_avalanches | show_tracking_trajectories
    };		 
				 
    void make_top_projection(const uint32_t hits_option_ = show_hit_all,
			     const uint32_t tracking_option_ = show_tracking_all) const;
    
    void build_event(const uint32_t tracking_option_ = show_tracking_all) const;	

    datatools::logger::priority verbosity = datatools::logger::PRIO_FATAL;

  private:


    void make_top_projection(const uint32_t solution_id,
                             const uint32_t hits_option_,
                             const uint32_t tracking_option_) const;
    
    void build_event(const uint32_t solution_id,
                     const uint32_t tracking_option_) const;
        
    const Geometry & _geom_;
    const Event * _event_ = nullptr;
    
  };

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_VISU_H
