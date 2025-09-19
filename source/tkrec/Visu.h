#ifndef FALAISE_CIMRMAN_VISU_H
#define FALAISE_CIMRMAN_VISU_H

// Standard headers
#include <functional>

// Bayeux:
#include <bayeux/datatools/bit_mask.h>

#include "tkrec/Geometry.h"
#include "tkrec/Event.h"

namespace tkrec {

  class Visu
  {
  public:
    
    Visu(const Geometry & geom_);

    bool has_event() const;
    
    void set_event(const Event & _event_);

    // vizualization section
    enum hit_flags {
      show_hit_used = datatools::bit_mask::bit00, ///< Show used hits for reconstruction (red)
      show_hit_unused = datatools::bit_mask::bit01, ///< Show unused hits for reconstruction (yellow)
      show_hit_associated = datatools::bit_mask::bit02, ///< Show associated hits to track (green)
      show_hit_unassociated_vert = datatools::bit_mask::bit03, ///< Show failed vertical position reconstruction but good drift radius (unassociated) (magenta)
      show_hit_associated_vert = datatools::bit_mask::bit04, ///< Show associated hits to track, good vertical position (green)
      show_hit_associated_novert = datatools::bit_mask::bit04, ///< Show associated hits to track, failed vertical position (teal)
      show_hit_all = show_hit_used | show_hit_unused | show_hit_associated | show_hit_unassociated_vert | show_hit_associated_vert | show_hit_associated_novert
    };
    
    // tracker hits options:
    // 0 - no unused hits	
    //	    red	= used hits for reconstruction    
    //
    //	1 - red	= used hits for reconstruction
    //	    yellow 	= unused hits for recontstruction
    //
    //	2 - red	= used hits for reconstruction (unassociated)
    //	    yellow 	= unused hits for recontstruction 
    //	    green 	= associated hits to track
    //
    // 	3 - red	= used hits for reconstruction (unassociated + good vertical position)
    //	    yellow 	= unused hits for recontstruction
    //	    magenta 	= failed vertical position reconstruction but good drift radius (unassociated)
    //	    green 	= associated hits to track, good vertical position
    //	    teal	= associated hits to track, failed vertical position
				 
    // tracking options:
    enum tracking_flags  {
      show_tracking_tracks = datatools::bit_mask::bit00,  ///< Show tracks // TODO? do I need this?
      show_tracking_hit_avalanches = datatools::bit_mask::bit01,  ///< Show hit avalanches
      show_tracking_trajectories = datatools::bit_mask::bit02, ///< Show trajectories
      show_tracking_all = show_tracking_tracks | show_tracking_hit_avalanches | show_tracking_trajectories
    };

    //	0 - only tracks
    //	1 - tracks
    //	    reconstructed tracker hit avalanche origin points
    //	2 - trajectories
    //	    avalanche origin points
    //	3 - tracks
    //	    trajectories
    //	    avalanche origin points				 
				 
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
