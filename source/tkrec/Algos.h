#ifndef FALAISE_CIMRMAN_ALGOS_H
#define FALAISE_CIMRMAN_ALGOS_H

#include <bayeux/datatools/logger.h>
#include <bayeux/datatools/properties.h>

// TK headers
#include "tkrec/Event.h"
#include "tkrec/Geometry.h"
#include "tkrec/Visu.h"


namespace tkrec {

  /// Event tracking reconstruction mode
  enum class ElectronRecMode
    {
      undefined,
      polyline,  /// electron polyline trajectory reconstruction
      line,  /// electron line trajectory reconstruction
    };

  ElectronRecMode event_recmode_from_label(const std::string & label_);

  // Config for clustering algorithms
  struct ClusteringConfing
  {
    bool save_sinograms = false;
    double max_distance = 3.0 * 44.0; ///< mm
    double hit_association_distance = 6.0; ///< mm
    uint32_t no_iterations = 2u; ///< dimensionless
    uint32_t resolution_phi = 100u; // No bins
    uint32_t resolution_r = 250u; // No bins
    double max_initial_precision_r = 6.0; ///< mm
    double zoom_factor = 10.0; ///< dimensionless
    double uncertainty = 2.0; ///< mm 
  };

  /// Config for kinked trajectory reconstruction algorithms
  struct PolylinesConfig
  {  
    double max_extention_distance = 120.0;  ///< mm 
    double max_vertical_distance = 4.0; ///< mm 
    double min_tracker_hits_distance = 100.0; ///< mm
    double max_kink_angle = 120.0; // degrees
    double max_trajectories_middlepoint_distance = 10.0; ///< mm 
    double max_trajectory_endpoints_distance = 75.0; ///< mm 
    double max_trajectory_connection_angle = 40.0; ///< degrees
    double min_distance_from_foil = 75.0; ///< mm 
    double min_distance_from_mainwalls = 75.0; ///< mm // TODO implement this
  };
  
  /// Config for kinked trajectory reconstruction algorithms
  struct AlphaConfig
  {  
    uint32_t clustering_resolution_phi = 100u; // No bins  
    bool save_sinograms = false;
    uint32_t resolution_r = 100; // No bins
    double phi_step = 0.5; // degrees 
    double max_r = 30.0; // mm
    double time_step = 100; // ns
    double uncertainty = 2.0; // mm
    double min_possible_drift_time = 0.0; // in nanoseconds
    double max_possible_drift_time = 5000.0; // in nanoseconds
  };

  /// Configuration parameters for event tracking reconstruction
  struct CimrmanAlgoConfig
  {
    datatools::logger::priority verbosity = datatools::logger::PRIO_FATAL;
    ElectronRecMode electron_mode = ElectronRecMode::undefined;
    
    bool use_provided_preclustering = false;
    bool reconstruct_alphas = false;
    bool visualization_2D = false;
    bool visualization_3D = false;
    bool force_default_sigma_r = false; ///< Flag to force the default sigma r value for tracker hits
    double default_sigma_r = 2.0; ///< implicit unit: mm
    double chi_square_threshold = 5.0; ///< dimensionless

    // For electron reconstruction modes
    ClusteringConfing clustering;      
    PolylinesConfig polylines;
    AlphaConfig alphas;

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
    void initialize(const CimrmanAlgoConfig & evrecconf_);
    void reset();
    void process(Event & event_);

  private:
    
    // step 1: preclustering
    void precluster();
    std::vector<std::vector<ConstTrackerHitHdl>> separate_hits(const std::vector<ConstTrackerHitHdl> & hits);


    // step 2: clustering
    void Legendre_transform_clustering();
    void clusterize_precluster(std::vector<ConstTrackerHitHdl> & tracker_hits,
                               std::vector<ClusterHdl> & clusters);
    void find_cluster_Legendre(const std::vector<ConstTrackerHitHdl> & hits, 
                               double & phi_estimate,
                               double & r_estimate) const;
    void separate_close_hits_to_line(std::vector<ConstTrackerHitHdl> & hits,
                                     std::vector<ConstTrackerHitHdl> & hits_separated,
                                     const double phi,
                                     const double r,
                                     const double distance_threshold);

    // step 2: alpha clustering and track estimation
    void alpha_clustering();
    std::pair<double, double> find_alpha_cluster(const std::vector<ConstTrackerHitHdl> & hits);
    void estimate_alpha_track(AlphaClusterHdl alpha_cluster);
    void find_reference_time_bounds(AlphaClusterHdl alpha_cluster);


    // step 3: MLM line fitting + ambiguity checking and solving
    void make_MLM_fits();
    void make_ML_estimate(LinearFitHdl fit);
    void detect_ambiguity_type(ClusterHdl cluster);
    void create_mirror_fit(ClusterHdl cluster);
    
    
    // step 4: Linear fits are associated to tracker hits and combined into a precluster solutions
    void combine_into_precluster_solutions();
    
    
    // step 5: Kink finding and connecting into polyline trajectories
    void create_polyline_trajectories();
    enum TrajectoryType{electron, alpha, both};
    void create_line_trajectories(TrajectoryType type);
    void create_line_trajectory_points(TrajectoryHdl & trajectory);
    std::vector<std::vector<TrackHdl>> find_polyline_candidates(std::vector<TrackHdl> & tracks,
                                                                const int side) const;
    void build_polyline_trajectory(PreclusterSolutionHdl & precluster_solution, 
                                   TrajectoryHdl & trajectory);    
                                   

    // step 6: Trajectory refinement (MLM, clustering refinement...)
    void refine_trajectories();
    void connect_close_trajectories(PreclusterSolutionHdl & precluster_solution);
    void remove_wrong_hits_associations(PreclusterSolutionHdl & precluster_solution, TrajectoryHdl & trajectory);
    void refine_clustering(PreclusterSolutionHdl & precluster_solution);
    void evaluate_trajectory(TrajectoryHdl & trajectory);
    
    
    // step 7: Combining precluster solutions into all solutions
    void create_solutions();
    void sort_solutions(std::vector<SolutionHdl> & solutions);
                   

  private:

    const Geometry & _geom_; ///< Geometry informations
    CimrmanAlgoConfig _config_; ///< Configuration
    Event * _event_ = nullptr; ///< Working event to be reconstructed
    std::unique_ptr<Visu> _visu_; ///< Visualisation engine
    
  };

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_ALGOS_H
