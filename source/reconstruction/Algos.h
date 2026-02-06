#ifndef FALAISE_CIMRMAN_ALGOS_H
#define FALAISE_CIMRMAN_ALGOS_H

// Standard headers
#include <iomanip>
#include <vector>
#include <memory>

// Boost:
#include <boost/multi_array.hpp>

// Bayeux:
#include <bayeux/datatools/exception.h>
#include <bayeux/datatools/clhep_units.h>
#include <bayeux/datatools/logger.h>
#include <bayeux/datatools/properties.h>

// Cimrman headers
#include "reconstruction/KinkFinder.h"
#include "reconstruction/Sinogram.h"
#include "geometry/Geometry.h"
#include "visualization/Visu.h"

#include "datamodel/Event.h"
#include "datamodel/Association.h"
#include "datamodel/DelayedCluster.h"
#include "datamodel/TrackerHit.h"
#include "datamodel/Precluster.h"
#include "datamodel/Solution.h"
#include "datamodel/PreclusterSolution.h"
#include "datamodel/LinearFit.h"
#include "datamodel/Track.h"
#include "datamodel/Trajectory.h"
#include "datamodel/Point.h"

namespace cimrman {

namespace dm = datamodel;
  
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
    uint32_t resolution_phi = 100u; // Number of bins
    uint32_t resolution_r = 200u; // Number of bins
    double max_initial_precision_r = 6.0 * CLHEP::mm; 
    double zoom_factor = 10.0;
    double uncertainty = 3.0 * CLHEP::mm;
  };

  /// Config for alpha clustering algorithms
  struct AlphaConfig
  {  
    uint32_t clustering_resolution_phi = 100u; // Number of bins  
    double min_possible_drift_time = 0.0 * CLHEP::ns; 
    double max_possible_drift_time = 5000.0 * CLHEP::ns; 
    double time_step = 50.0 * CLHEP::ns; 

    bool save_sinograms = false;
    uint32_t iterations = 2u;
    double zoom_factor = 10.0;
    double phi_step = 2.0 * CLHEP::deg;
    uint32_t resolution_r = 30u; // Number of bins
    double delta_r = 60.0 * CLHEP::mm;
    double uncertainty = 2.0 * CLHEP::mm;
  };

  /// Config for kinked trajectory reconstruction algorithms
  struct PolylinesConfig
  {  
    // sharp kink reconstruction
    double max_vertical_distance = 40.0 * CLHEP::mm; 
    double max_tracker_hits_distance = 100.0 * CLHEP::mm; 
    double min_kink_angle = 0.0 * CLHEP::deg;
    double max_kink_angle = 120.0 * CLHEP::deg;
    double min_distance_from_foil = 75.0 * CLHEP::mm; 
    double min_distance_from_main_walls = 50.0 * CLHEP::mm;
    double min_distance_from_X_walls = 50.0 * CLHEP::mm;
    
    // clustering refinements
    double max_extention_distance = 120.0 * CLHEP::mm;  

    // small kink reconstruction
    double max_trajectories_middlepoint_distance = 15.0 * CLHEP::mm; 
    double max_trajectory_endpoints_distance = 100.0 * CLHEP::mm;
    double min_trajectory_connection_angle = 0.0 * CLHEP::deg;
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
    bool keep_only_best_solutions = false;
    bool force_default_sigma_r = false; ///< Flag to force the default sigma r value for tracker hits
    double default_sigma_r = 2.0 * CLHEP::mm;
    
    // chi2 limiter is not implemented
    //bool discard_bad_fits = false;
    //double chi_square_threshold = 20.0; ///< dimensionless

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
    void set_event(dm::Event & event_);
    bool has_event() const;
    
    bool is_initialized() const;
    void initialize(const CimrmanAlgoConfig & evrecconf_);
    void reset();
    void process(dm::Event & event_);

  private:
    
    // step 1: preclustering
    void precluster();
    std::vector<std::vector<dm::TrackerHitHdl>> separate_hits(const std::vector<dm::TrackerHitHdl> & hits);


    // step 2: clustering
    void electron_clustering(dm::PreclusterHdl precluster);
    void recursive_clusterizer(std::vector<dm::TrackerHitHdl> & tracker_hits,
                               std::vector<dm::ClusterHdl> & clusters);
    void find_cluster_Legendre(const std::vector<dm::TrackerHitHdl> & hits, 
                               double & phi_estimate,
                               double & r_estimate);
                               
    void separate_close_hits_to_line(std::vector<dm::TrackerHitHdl> & hits,
                                     std::vector<dm::TrackerHitHdl> & hits_separated,
                                     const double phi,
                                     const double r,
                                     const double distance_threshold);
    std::vector<std::pair<double, double>> estimate_three_hit_tracks(const std::vector<dm::TrackerHitHdl>& hits); // different approach (experimental)
                                     

    // step 2: alpha clustering and track estimation
    void alpha_clustering(dm::PreclusterHdl precluster);
    std::pair<double, double> find_delayed_cluster(std::vector<dm::TrackerHitHdl> & hits,
                                                   std::vector<dm::TrackerHitHdl> & cluster_hits);
    void estimate_delayed_track(dm::DelayedClusterHdl delayed_cluster);
    void find_reference_time_bounds(dm::DelayedClusterHdl delayed_cluster);


    // step 3: MLM line fitting + ambiguity checking and solving
    void make_MLM_fits(dm::PreclusterHdl precluster);
    void make_ML_estimate(dm::LinearFitHdl fit);
    void detect_ambiguity_type(dm::ClusterHdl cluster);
    void create_mirror_fit(dm::ClusterHdl cluster);
    
    
    // step 4: Linear fits are associated to tracker hits and combined into a precluster solutions
    void combine_into_precluster_solutions(dm::PreclusterHdl precluster);
    
    
    // step 5: Kink finding and connecting into polyline trajectories
    void create_trajectories(dm::PreclusterHdl precluster);
    
    void create_polyline_trajectories(dm::PreclusterHdl precluster);
    void create_line_trajectories(dm::PreclusterHdl precluster,
                                  dm::Trajectory::Type type = dm::Trajectory::Type::ANY);
    void create_line_trajectory_points(dm::TrajectoryHdl trajectory);
    std::vector<std::vector<dm::TrackHdl>> find_polyline_candidates(std::vector<dm::TrackHdl> & tracks,
                                                                const int side) const;
    void build_polyline_trajectory(dm::PreclusterSolutionHdl precluster_solution, 
                                   dm::TrajectoryHdl trajectory);    
                                   

    // step 6: Trajectory kinks and clustering refinements
    void refinements(dm::PreclusterHdl precluster);
    void remove_out_of_bound_associations(dm::PreclusterSolutionHdl precluster_solution, 
                                          dm::TrajectoryHdl trajectory);
    void connect_trajectories(dm::PreclusterSolutionHdl precluster_solution,
                              KinkFinder::ConnectionStrategy strategy);
    void merge_trajectories_into_base(dm::TrajectoryHdl base_trajectory,
                                      dm::Trajectory::EndPoint base_endpoint, 
                                      dm::TrajectoryHdl other_trajectory, 
                                      dm::Trajectory::EndPoint other_endpoint,
                                      dm::PointHdl kink_point);
    
    void refine_clustering(dm::PreclusterSolutionHdl precluster_solution);
    
    
    
    // step 7: Evaluating trajectories
    void evaluate_trajectories(dm::PreclusterHdl precluster);
    void evaluate_trajectory(dm::TrajectoryHdl trajectory);
    
    
    // step 8: Combining precluster solutions into all solutions
    void create_solutions();
    void combine_precluster_solutions_into_solutions();
    void evaluate_solutions();
    static bool compare_solutions(const dm::SolutionHdl solution1, const dm::SolutionHdl solution2);
    void remove_suboptimal_solutions();

  private:

    const Geometry & _geom_; ///< Geometry informations
    CimrmanAlgoConfig _config_; ///< Configuration
    dm::Event * _event_ = nullptr; ///< Working event to be reconstructed
    std::unique_ptr<Visu> _visu_; ///< Visualisation engine
    
    Sinogram prompt_sinogram_manager; // support worker class for calculation of prompt hits sinograms 
    Sinogram delayed_sinogram_manager; // support worker class for calculation of delayed hits sinograms 
    
    friend class KinkFinder;
  };

} //  end of namespace cimrman

#endif // FALAISE_CIMRMAN_ALGOS_H
