#ifndef FALAISE_TKRECONSTRUCT_ALGOS_H
#define FALAISE_TKRECONSTRUCT_ALGOS_H

#include <bayeux/datatools/logger.h>
#include <bayeux/datatools/properties.h>

// TK headers
#include "tkrec/Event.h"
#include "tkrec/Geometry.h"
#include "tkrec/Visu.h"

namespace tkrec {

  /// Event tracking reconstruction mode
  enum class EventRecMode
    {
      undefined,
      electron_kinked,  /// electron polyline trajectory reconstruction
      electron_straight,  /// electron line trajectory reconstruction
    };

  EventRecMode event_recmode_from_label(const std::string & label_);

  /// Configuration parameters for event tracking reconstruction
  struct TKEventRecConfig
  {
    datatools::logger::priority verbosity = datatools::logger::PRIO_FATAL;
    EventRecMode mode = EventRecMode::undefined;
    bool reconstruct_alphas = true;
    bool save_sinograms = false;
    bool visualization = false;

    bool force_default_sigma_r = false; ///< Flag to force the default sigma r value for tracker hits


    // For electron reconstruction modes
    double default_sigma_r = 2.0; ///< implicit unit: mm
    double chi_square_threshold = 5.0; ///< dimensionless
    
    //struct clustering //TODO I dont know how to do this
    //{
      double clustering_max_distance = 3.0 * 44.0; ///< mm
      double clustering_hit_association_distance = 6.0; ///< mm
      uint32_t clustering_no_iterations = 2u; ///< dimensionless
      uint32_t clustering_resolution_phi = 100u;
      uint32_t clustering_resolution_r = 250u;
      double clustering_max_initial_precision_r = 6.0; ///< mm
      double clustering_zoom_factor = 10.0; ///< dimensionless
      double clustering_uncertainty = 2.0; ///< mm 
    //};

    //struct polylines
    //{  
      double polylines_max_extention_distance = 120.0;  ///< mm 
      double polylines_max_vertical_distance = 4.0; ///< mm 
      double polylines_min_tracker_hits_distance = 100.0; ///< mm
      double polylines_max_kink_angle = 120; // degrees??
      double polylines_max_trajectories_middlepoint_distance = 10.0;
      double polylines_max_trajectory_endpoints_distance = 75.0;
      double polylines_max_trajectory_connection_angle = 40;
      double polylines_min_distance_from_foil = 75.0;
    //};
    void parse(const datatools::properties & config_);
    
  };

  /// Main cluster/track reconstruction class.
  /// This class implements several algorithms.
  class Algos
  {
  public:

    Algos(const Geometry & geom_);
    ~Algos();

    // Public interface:
    void set_event(Event & event_);
    bool has_event() const;
    
    bool is_initialized() const;
    void initialize(const TKEventRecConfig & evrecconf_);
    void reset();
    void process(Event & event_);

  private:
    
    // step 1: preclustering
    void precluster();
    std::vector<std::vector<ConstTrackerHitHdl>> separate_hits(const std::vector<ConstTrackerHitHdl> & hits);
    bool clusters_close(const std::vector<ConstTrackerHitHdl> & cluster1,
     			              const std::vector<ConstTrackerHitHdl> & cluster2) const;
	
    // step 2: clustering
    void Legendre_transform_cluster_finder();
    void clusterize(std::vector<ConstTrackerHitHdl> & tracker_hits, std::vector<ClusterHdl> & clusters);
    void find_cluster_Legendre(const std::vector<ConstTrackerHitHdl> & hits, 
                               double & phi_estimate,
                               double & r_estimate) const;
    void separate_close_hits_to_line(std::vector<ConstTrackerHitHdl> & hits,
                                     std::vector<ConstTrackerHitHdl> & hits_separated,
                                     const double phi,
                                     const double r,
                                     const double distance_threshold);

    // step 3: MLM line fitting + ambiguity checking and solving
    void make_MLM_fits();
    void detect_ambiguity_type(ClusterHdl cluster);
    void create_mirror_fit(ClusterHdl cluster);
    void make_ML_estimate(LinearFitHdl fit);
    
    // step 4: Linear fits are associated to tracker hits and combined into a precluster solutions
    void combine_into_precluster_solutions();
    
    // step 5: Kink finding and connecting into polyline trajectories
    void create_polyline_trajectories();
    void create_line_trajectories();
    void create_line_trajectory_points(TrajectoryHdl & trajectory);
    std::vector<std::vector<TrackHdl>> find_polyline_candidates(std::vector<TrackHdl> & tracks,
                                                      const int side) const;
    void build_polyline_trajectory(PreclusterSolutionHdl & precluster_solution, TrajectoryHdl & trajectory);    
    void connect_close_trajectories(PreclusterSolutionHdl & precluster_solution);
    //void create_polyline_trajectory_points(TrajectoryHdl & trajectory);
    //void split_fake_polyline_candidates(std::vector<std::vector<TrackHdl>> & trajectory_candidates );
    //void remove_fake_segments(std::vector<std::vector<TrackHdl>> & trajectory_candidates );

    // step 6: Trajectory refinement (MLM, clustering refinement...)
    void refine_trajectories();
    void remove_wrong_hits_associations(PreclusterSolutionHdl & precluster_solution, TrajectoryHdl & trajectory);
    void refine_clustering(PreclusterSolutionHdl & precluster_solution);
    void evaluate_trajectory(TrajectoryHdl & trajectory);
    
    // step 7: Combining precluster solutions into all solutions
    void create_solutions();
    void sort_solutions(std::vector<SolutionHdl> & solutions);
	    
  private:

    void _process_electron_kinked_();
    void _process_electron_straight_();
    void _process_alpha_();


    const Geometry & _geom_; ///< Geometry informations
    TKEventRecConfig _config_; ///< Configuration
    Event * _event_ = nullptr; ///< Working event to be reconstructed
    std::unique_ptr<Visu> _visu_; ///< Visualisation engine
    
  };

} //  end of namespace tkrec

#endif // FALAISE_TKRECONSTRUCT_ALGOS_H
