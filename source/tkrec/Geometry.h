#ifndef FALAISE_CIMRMAN_GEOMETRY_H
#define FALAISE_CIMRMAN_GEOMETRY_H

#include <falaise/snemo/geometry/calo_locator.h>
#include <falaise/snemo/geometry/gveto_locator.h>
#include <falaise/snemo/geometry/xcalo_locator.h>
#include <falaise/snemo/geometry/gg_locator.h>
#include <falaise/snemo/geometry/locator_helpers.h>
#include <falaise/snemo/geometry/locator_plugin.h>

namespace tkrec {

  struct Geometry
  {
    
    // dimensions in mm
    // origin in the center of detector

    bool has_Bi_source = false;
    double Bi_source_x = 1.4;
    double Bi_source_y = 14.2;
    double Bi_source_z = 21.1;
    double Bi_source_dist_y = 835.0;
    double Bi_source_dist_z = 425.0;

    // dimensions in mm
    // origin in the center of detector
    double tc_radius = 22.0;

    // OM dimensions in mm
    double mw_sizex = 194.0;
    double mw_sizey = 256.0;
    double mw_sizez = 256.0;

    double gv_sizex = 308.0;
    double gv_sizey = 310.0;
    double gv_sizez = 150.0;

    double xw_sizex = 200.0;
    double xw_sizey = 150.0;
    double xw_sizez = 208.5;
    
    const snemo::geometry::locator_plugin * geo_loc;

  };

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_GEOMETRY_H
