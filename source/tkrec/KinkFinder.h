#ifndef FALAISE_CIMRMAN_KINKFINDER_H
#define FALAISE_CIMRMAN_KINKFINDER_H

// Standard headers
#include <vector>
#include <memory>

// Cimrman headers
#include "tkrec/Point.h"

namespace tkrec {
  
  // Forward declarations
  class Track;
  using TrackHdl = std::shared_ptr<Track>;
  using ConstTrackHdl = std::shared_ptr<const Track>;
  
  struct PolylinesConfig;
  struct Geometry;
  
  class KinkFinder
  {		
  public:
  
    enum Status{ UNPROCESSED, CANDIDATE_LOCATED, KINK_FOUND, KINK_NOT_FOUND };
    // KinkFinder is initialized as UNPROCESSED (kink_point is empty),
    // running some connection strategy constructs the potential kink_point (CANDIDATE_LOCATED),
    // after connection strategy finishes, kink is either confirmed (KINK_FOUND)
    // or rejected (KINK_NOT_FOUND)
  
    enum ConnectionStrategy{ VERTICAL_ALIGNMENT, ENDPOINTS_MIDDLE };
    // ConnectionStrategy defines which strategy and set of criteria (which procedure) is chosen by the process function
    // VERTICAL_ALIGNMENT: connection of the track is enforced by changing only vertical parts of the fit
    //                     (2D intersection of the tracks is used)
    // ENDPOINTS_MIDDLE: kink is formed as a middle-point between the two endpoints to be connected

  private:
   
    TrackHdl track1 = nullptr;
    PointHdl endpoint1 = nullptr;
    
    TrackHdl track2 = nullptr;
    PointHdl endpoint2 = nullptr;

    Point kink_point;
    
    Status status = Status::UNPROCESSED;
    const PolylinesConfig & config;
    const Geometry & geom;
    
    
  private:
    
    // angular criteria
    bool check_2D_angle( const double min_angle, const double max_angle ) const;
    bool check_3D_angle( const double min_angle, const double max_angle ) const;
        
    // position criteria
    bool check_is_inside_tracker() const; 
    bool check_distance_to_MW( const double min_distance ) const;
    bool check_distance_to_XW( const double min_distance ) const;
    bool check_distance_to_SF( const double min_distance ) const;

    // associated tracker hit criteria
    bool check_close_hits_2D( const double max_angle ) const;
    bool check_are_ends_close_2D( const double max_distance ) const;
    bool check_are_ends_close_3D( const double max_distance ) const;

    // distance criteria specific to the connection strategies
    bool check_vertical_distance( const double max_distance) const;
    bool check_middle_point_distance( const double max_distance) const;
     
  public:
    
    KinkFinder() = default;

    KinkFinder( TrackHdl _track1, PointHdl _endpoint1,
                TrackHdl _track2, PointHdl _endpoint2,
                const PolylinesConfig & _config,
                const Geometry & _geom);

    // master function 
    void process( ConnectionStrategy strategy );
    
    // process runs one of these procedures based on the chosen connection strategy
    void large_kink_procedure(); 
    void small_kink_procedure();
    
    Status get_status() const;
    Point get_kink_point() const;
    
  };

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_KINKFINDER_H
