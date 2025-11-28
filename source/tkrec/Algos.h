#ifndef FALAISE_CIMRMAN_ALGOS_H
#define FALAISE_CIMRMAN_ALGOS_H

// Standard headers
#include <iomanip>
#include <algorithm>
#include <vector>
#include <stack>

// Root:
#include <TH2F.h>
#include <TCanvas.h>

// Boost:
#include <boost/multi_array.hpp>

// Bayeux:
#include <bayeux/datatools/exception.h>
#include <bayeux/datatools/clhep_units.h>
#include <bayeux/datatools/logger.h>
#include <bayeux/datatools/properties.h>

// Cimrman headers
#include "tkrec/Event.h"
#include "tkrec/Geometry.h"
#include "tkrec/Visu.h"
#include "tkrec/Sinogram.h"


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
    double max_distance = 130.0 * CLHEP::mm; // 3.0 * 44.0 = 132
    double hit_association_distance = 6.0 * CLHEP::mm;
    uint32_t iterations = 2u;
    uint32_t resolution_phi = 100u; // No bins
    uint32_t resolution_r = 200u; // No bins
    double max_initial_precision_r = 6.0 * CLHEP::mm; 
    double zoom_factor = 10.0;
    double uncertainty = 3.0 * CLHEP::mm;
  };

  /// Config for alpha clustering algorithms
  struct AlphaConfig
  {  
    uint32_t clustering_resolution_phi = 100u; // No bins  
    double min_possible_drift_time = 0.0 * CLHEP::ns; 
    double max_possible_drift_time = 5000.0 * CLHEP::ns; 
    double time_step = 50.0 * CLHEP::ns; 

    bool save_sinograms = false;
    uint32_t iterations = 2u;
    double zoom_factor = 10.0;
    double phi_step = 2.0 * CLHEP::deg;
    uint32_t resolution_r = 30u; // No bins
    double delta_r = 60.0 * CLHEP::mm;
    double uncertainty = 2.0 * CLHEP::mm;
  };

  /// Config for kinked trajectory reconstruction algorithms
  struct PolylinesConfig
  {  
    // sharp kink reconstruction
    double max_vertical_distance = 40.0 * CLHEP::mm; 
    double max_tracker_hits_distance = 100.0 * CLHEP::mm; 
    double max_kink_angle = 120.0 * CLHEP::deg;
    double min_distance_from_foil = 75.0 * CLHEP::mm; 
    double min_distance_from_main_walls = 50.0 * CLHEP::mm;
    double min_distance_from_X_walls = 50.0 * CLHEP::mm;
    
    // clustering refinements
    double max_extention_distance = 120.0 * CLHEP::mm;  

    // small kink reconstruction
    double max_trajectories_middlepoint_distance = 15.0 * CLHEP::mm; 
    double max_trajectory_endpoints_distance = 100.0 * CLHEP::mm;
    double max_trajectory_connection_angle = 40.0 * CLHEP::deg;
  };
  
  /// Configuration parameters for event tracking reconstruction
  struct CimrmanAlgoConfig
  {
    datatools::logger::priority verbosity = datatools::logger::PRIO_FATAL;
    ElectronRecMode electron_mode = ElectronRecMode::undefined;
    
    bool visualization_2D = false;
    bool visualization_3D = false;
    bool use_provided_preclustering = false;
    bool reconstruct_alphas = false;
    bool force_default_sigma_r = false; ///< Flag to force the default sigma r value for tracker hits
    double default_sigma_r = 2.0 * CLHEP::mm;
    double chi_square_threshold = 5.0; ///< dimensionless

    /// Configuration of clustering algorihtms
    ClusteringConfing clustering;
      
    /// Configuration of alpha clustering algorihtms
    AlphaConfig alphas;
      
    /// Configuration of polyline reconstruction algorihtms 
    PolylinesConfig polylines;
  
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
    std::vector<std::vector<TrackerHitHdl>> separate_hits(const std::vector<TrackerHitHdl> & hits);


    // step 2: clustering
    void Legendre_transform_clustering();
    void clusterize_precluster(std::vector<TrackerHitHdl> & tracker_hits,
                               std::vector<ClusterHdl> & clusters);
    void find_cluster_Legendre(const std::vector<TrackerHitHdl> & hits, 
                               double & phi_estimate,
                               double & r_estimate);
                               
    void separate_close_hits_to_line(std::vector<TrackerHitHdl> & hits,
                                     std::vector<TrackerHitHdl> & hits_separated,
                                     const double phi,
                                     const double r,
                                     const double distance_threshold);

    // step 2: alpha clustering and track estimation
    void alpha_clustering();
    std::pair<double, double> find_delayed_cluster(std::vector<TrackerHitHdl> & hits,
                                                   std::vector<TrackerHitHdl> & cluster_hits);
    void estimate_delayed_track(DelayedClusterHdl delayed_cluster);
    void find_reference_time_bounds(DelayedClusterHdl delayed_cluster);


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
                                   

    // step 6: Trajectory kinks and clustering refinements
    void refinements();
    void remove_wrong_hits_associations(PreclusterSolutionHdl & precluster_solution, TrajectoryHdl & trajectory);
    void connect_close_trajectories(PreclusterSolutionHdl & precluster_solution);
    void refine_clustering(PreclusterSolutionHdl & precluster_solution);
    
    // step 7: Evaluating trajectories
    void evaluate_trajectories();
    void evaluate_trajectory(TrajectoryHdl & trajectory);
    
    // step 8: Combining precluster solutions into all solutions
    void create_solutions();
    void sort_solutions(std::vector<SolutionHdl> & solutions);
                   

  private:

    const Geometry & _geom_; ///< Geometry informations
    CimrmanAlgoConfig _config_; ///< Configuration
    Event * _event_ = nullptr; ///< Working event to be reconstructed
    std::unique_ptr<Visu> _visu_; ///< Visualisation engine
    
    Sinogram prompt_sinogram_manager; // support worker class for calculation of prompt hits sinograms 
    Sinogram delayed_sinogram_manager; // support worker class for calculation of delayed hits sinograms 
    
  };

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_ALGOS_H
