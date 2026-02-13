#ifndef FALAISE_CIMRMAN_GEOMETRY_H
#define FALAISE_CIMRMAN_GEOMETRY_H

// Standard headers
#include <memory>
#include <array>

// Falaise headers
#include <falaise/snemo/geometry/calo_locator.h>
#include <falaise/snemo/geometry/gveto_locator.h>
#include <falaise/snemo/geometry/xcalo_locator.h>
#include <falaise/snemo/geometry/gg_locator.h>
#include <falaise/snemo/geometry/locator_helpers.h>
#include <falaise/snemo/geometry/locator_plugin.h>

#include <falaise/snemo/services/service_handle.h>
#include <falaise/snemo/services/geometry.h>

// Bayeux:
#include <bayeux/datatools/logger.h>
#include <bayeux/datatools/clhep_units.h>
#include <bayeux/geomtools/box.h>

namespace cimrman {
  
  // Forward declarations
  namespace datamodel {
    class Point;
    using PointHdl = std::shared_ptr<Point>;
  }


  class Geometry
  {
  private:
  
    const snemo::geometry::locator_plugin * geo_loc;
    const geomtools::mapping * geoMapping;
    const geomtools::id_mgr * geoIdMgr;
    
    const geomtools::placement * frenchTrkVolPlacement;
    const geomtools::placement * italianTrkVolPlacement;
    const geomtools::box * frenchTrkVolBox;
    const geomtools::box * italianTrkVolBox;
    
  public:    
    
    double tc_radius = 22.0 * CLHEP::mm;

    // Bi sources and OMs are only used by visualization class
    bool has_Bi_source = false;    
    double Bi_source_x = 1.4 * CLHEP::mm;
    double Bi_source_y = 14.2 * CLHEP::mm;
    double Bi_source_z = 21.1 * CLHEP::mm;
    double Bi_source_dist_y = 835.0 * CLHEP::mm;
    double Bi_source_dist_z = 425.0 * CLHEP::mm;

    double mw_sizex = 194.0 * CLHEP::mm;
    double mw_sizey = 256.0 * CLHEP::mm;
    double mw_sizez = 256.0 * CLHEP::mm;

    double gv_sizex = 308.0 * CLHEP::mm;
    double gv_sizey = 310.0 * CLHEP::mm;
    double gv_sizez = 150.0 * CLHEP::mm;

    double xw_sizex = 200.0 * CLHEP::mm;
    double xw_sizey = 150.0 * CLHEP::mm;
    double xw_sizez = 208.5 * CLHEP::mm;
    
    enum Side {IT, FR, ANY};
    
  public:
    
    datatools::logger::priority verbosity = datatools::logger::PRIO_FATAL;
    
    Geometry() = default;
    
    Geometry(snemo::service_handle<snemo::geometry_svc> & geoManager, datatools::logger::priority _verbosity);
    
    double distance_to_MW(const datamodel::Point & point /*, Side side*/) const;
    double distance_to_XW(const datamodel::Point & point /*, Side side*/) const;
    double distance_to_SF(const datamodel::Point & point /*, Side side*/) const;
    bool is_inside_tracker(const datamodel::Point & point, Side side = Side::ANY) const;
    
    std::array<double, 2> get_cell_position(const int SRL[3]) const;
    std::array<double, 3> get_MW_OM_position(const int SWCR[4]) const;
    std::array<double, 3> get_XW_OM_position(const int SWCR[4]) const;
    std::array<double, 3> get_GV_OM_position(const int SWCR[4]) const;
  };

} //  end of namespace cimrman

#endif // FALAISE_CIMRMAN_GEOMETRY_H
