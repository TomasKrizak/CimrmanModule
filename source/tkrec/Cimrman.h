#ifndef FALAISE_CIMRMAN_CIMRMAN_H
#define FALAISE_CIMRMAN_CIMRMAN_H

// Standard headers
#include <memory>
#include <iomanip>

// Interface from Falaise
#include "bayeux/dpp/base_module.h"
#include "bayeux/mctools/simulated_data.h"

// Third party:
// - Bayeux:
#include <bayeux/datatools/clhep_units.h>
#include <bayeux/datatools/logger.h>
#include <bayeux/datatools/handle.h>
#include <bayeux/geomtools/line_3d.h>
#include <bayeux/geomtools/polyline_3d.h>
#include <bayeux/geomtools/box.h>

#include <falaise/snemo/geometry/calo_locator.h>
#include <falaise/snemo/geometry/gveto_locator.h>
#include <falaise/snemo/geometry/xcalo_locator.h>
#include <falaise/snemo/geometry/gg_locator.h>
#include <falaise/snemo/geometry/locator_helpers.h>
#include <falaise/snemo/geometry/locator_plugin.h>

#include "falaise/snemo/processing/module.h"
#include "falaise/snemo/datamodels/calibrated_data.h"
#include "falaise/snemo/datamodels/calibrated_tracker_hit.h"
#include "falaise/snemo/datamodels/data_model.h"
#include "falaise/snemo/datamodels/event_header.h"
#include "falaise/snemo/datamodels/tracker_clustering_data.h"
#include "falaise/snemo/datamodels/tracker_clustering_solution.h"
#include "falaise/snemo/datamodels/tracker_trajectory_data.h"
#include "falaise/snemo/datamodels/tracker_trajectory_solution.h"
#include "falaise/snemo/datamodels/line_trajectory_pattern.h"
#include "falaise/snemo/datamodels/polyline_trajectory_pattern.h"
#include "falaise/snemo/datamodels/particle_track_data.h"
#include "falaise/snemo/services/geometry.h"
#include "falaise/snemo/services/service_handle.h"

#include "tkrec/Event.h"
#include "tkrec/Geometry.h"
#include "tkrec/Algos.h"

namespace tkrec {

  /// Main event reconstruction module
  class Cimrman : public dpp::base_module
  {
  public:

    /// \brief Configuration parameters
    struct config_type
    {
      /// Verbosity
      datatools::logger::priority verbosity = datatools::logger::PRIO_FATAL;
      
      /// Label of the input CD bank
      std::string CD_label = "CD";
      
      /// Label of the output TCD bank
      std::string TCD_label = "TCD";
      
      /// Label of the output TTD bank
      std::string TTD_label = "TTD";
      
      /// Configuration of the Cimrman reconstruct algorihtms (from Algos.h)
      CimrmanAlgoConfig recConfig;
      
    };
    
    ////////////////////////////////////////////////
    // The following PUBLIC methods MUST be defined!
    Cimrman();

    virtual ~Cimrman();

    //! Return a const reference to the module's configuration
    const config_type & config() const;

    //! Read configuration from parameters config
    void read_config(const datatools::properties& config_);

    //! Configure the module
    virtual void initialize(const datatools::properties &myConfig,
			    datatools::service_manager &flServices,
			    dpp::module_handle_dict_type &what);

    //! Reset the module
    virtual void reset();

    //! Process event
    virtual dpp::base_module::process_status process(datatools::things &workItem);
    
 
  private:
 
    // Internal methods 
    
    void _populate_working_event_(const datatools::things & workItem);
  
    void _fill_TCD_bank_(const snemo::datamodel::calibrated_data & falaiseCDbank,
			  snemo::datamodel::tracker_clustering_data & the_tracker_clustering_data) const;

    void _fill_TTD_bank_(snemo::datamodel::tracker_clustering_data & the_tracker_clustering_data,
			  snemo::datamodel::tracker_trajectory_data & the_tracker_trajectory_data) const;
			 
    void _remove_duplicate_clustering_solutions_(snemo::datamodel::tracker_clustering_data & the_tracker_clustering_data,
        snemo::datamodel::tracker_trajectory_data & the_tracker_trajectory_data) const;

    /// Set/initialize default internal resources for the reconstruction module
    void _set_defaults_();

    /// Initialize geometry informations for the TK algos (needs the Falaise's geometry manager, see below)
    void _init_geom_(Geometry & geom_);

    /// Configuration parameters
    config_type _config_;

    // Working resources:
    snemo::service_handle<snemo::geometry_svc> _geoManager_; //!< The geometry manager

    // Working private materials (hidden to the public interface, using the PIMPL idiom)
    struct pimpl_type;
    friend struct pimpl_type;
    std::unique_ptr<pimpl_type> _work_; ///< Embedded resources (data and algo)
    
    DPP_MODULE_REGISTRATION_INTERFACE(Cimrman)
    
  };
  
} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_CIMRMAN_H
