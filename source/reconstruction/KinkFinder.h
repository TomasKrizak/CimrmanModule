#ifndef FALAISE_CIMRMAN_KINKFINDER_H
#define FALAISE_CIMRMAN_KINKFINDER_H

// Standard headers
#include <vector>
#include <memory>
#include <functional>

// Cimrman headers
#include "datamodel/Point.h"

namespace cimrman {
  
  // Forward declarations
  namespace datamodel {
    class Track;
    using TrackHdl = std::shared_ptr<Track>;
    using ConstTrackHdl = std::shared_ptr<const Track>;
  }
  
  struct Geometry;
  struct PolylinesConfig;
  
  
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
   
    datamodel::TrackHdl track1 = nullptr;
    datamodel::PointHdl endpoint1 = nullptr;
    
    datamodel::TrackHdl track2 = nullptr;
    datamodel::PointHdl endpoint2 = nullptr;

    datamodel::Point kink_point;
    
    Status status = Status::UNPROCESSED;
    const PolylinesConfig & config;
    const Geometry & geom;
    
    
  private:
    
    // angular criteria
    bool check_angle_2D_after( const double min_angle, const double max_angle ) const;
    bool check_angle_3D_after( const double min_angle, const double max_angle ) const;
    bool check_angle_2D_before( const double min_angle, const double max_angle ) const;
    bool check_angle_3D_before( const double min_angle, const double max_angle ) const;
        
    // position criteria
    bool check_is_inside_tracker() const; 
    bool check_distance_to_MW( const double min_distance ) const;
    bool check_distance_to_XW( const double min_distance ) const;
    bool check_distance_to_SF( const double min_distance ) const;

    // associated tracker hit criteria
    bool check_close_hits_2D( const double max_angle ) const;
    
    // endpoints criteria
    bool check_are_ends_close_2D( const double max_distance ) const;
    bool check_are_ends_close_3D( const double max_distance ) const;

    // distance criteria specific to the connection strategies
    bool check_vertical_distance( const double max_distance) const;
    bool check_middle_point_distance( const double max_distance) const;
    
    
    // checks a list of criteria 
    using Criterion = std::function<bool()>;
    bool check_criteria(std::initializer_list<Criterion> criteria) const; 

    // process runs one of these procedures based on the chosen connection strategy
    void vertical_alignment_procedure(); 
    void endpoints_middle_procedure();
     
  public:
    
    KinkFinder() = default;

    KinkFinder( datamodel::TrackHdl _track1, datamodel::PointHdl _endpoint1,
                datamodel::TrackHdl _track2, datamodel::PointHdl _endpoint2,
                const PolylinesConfig & _config,
                const Geometry & _geom);

    void process( ConnectionStrategy strategy );
    
    Status get_status() const;
    datamodel::Point get_kink_point() const;
    
  };

} //  end of namespace cimrman

#endif // FALAISE_CIMRMAN_KINKFINDER_H
