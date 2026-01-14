#ifndef FALAISE_CIMRMAN_SINOGRAM_H
#define FALAISE_CIMRMAN_SINOGRAM_H

// Standard headers
#include <vector>
#include <tuple>
#include <memory>
#include <cmath>

namespace tkrec {
  
  // Forward declarations
  class TrackerHit;
  using TrackerHitHdl = std::shared_ptr<TrackerHit>;
  using ConstTrackerHitHdl = std::shared_ptr<const TrackerHit>;
  
  
  class Sinogram
  {		
  private:
  
    // settings of the peak finder algo
    struct Settings
    {
      bool save_sinograms;
      uint32_t resolution_phi;
      uint32_t resolution_r;
      uint32_t iterations;
      double zoom_factor;
      double sigma;
      //double max_precision_r;
    };
    
    friend class Algos;
    
    Settings settings;
    Settings & get_settings(); 

    // array for the sinogram and tracker hits to investigate
    std::vector<float> dual_space;
    std::vector<TrackerHitHdl> tracker_hits;
    
    // currently investigated parametric space
    double delta_phi = M_PI;
    double delta_r = 2500.0; 
    double center_X = 0.0;
    double center_Y = 0.0;
    
    // current state of the sinogram 
    uint32_t iteration = 0u;
    double peak_phi = 0.0;
    double peak_r = 0.0;
    double peak_value = 0.0;

    // sinogram identifiers for saving png images 
    int run_number = 0; 
    int event_number = 0;
    int sinogram_ID = 0;
    bool prompt;
    
    void add_hit_to_sinogram(TrackerHitHdl & tracker_hit);
    inline double index_to_phi(uint32_t index_phi) const;
    inline double index_to_r(uint32_t index_r) const;

  public:
    
    Sinogram() = default;

    // master function 
    void process();
    
    // returns the (peak_phi, peak_r, peak_value)
    std::tuple<double, double, double> get_peak();
    
    // saves a png of the current state of the sinogram
    void save_sinogram() const;
    void set_run(int _run_number);
    void set_event(int _event_number); 

    void reset_dual_space(); // resets the internal buffer
    void reset_state(double _peak_phi = 0.0, double _peak_r = 0.0); // resets internal state
    void reset_sinogram_counter(); // resets sinogram_ID
    
    void set_delta_phi(double _delta_phi);
    void set_delta_r(double _delta_r);
    void set_center(double _center_X, double _center_Y);
    
    double get_delta_phi() const;
    double get_delta_r() const;
    double get_center_X() const;
    double get_center_Y() const;
    
    std::vector<float> & get_dual_space();
    const std::vector<float> & get_dual_space() const;
    void set_tracker_hits(const std::vector<TrackerHitHdl> & _tracker_hits);
    
  };

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_SINOGRAM_H
