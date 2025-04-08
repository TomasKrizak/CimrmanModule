TKReconstruct is a first implementation of the new tracking as a Falaise module. Currently it is based on TKEvent library. This is only a testing version, not the final implementation!

TKReconstruct module takes CD bank as an input a creates a reconstruction using TKEvent library. The found clustering solution and trajectory solution is stored in TCD and TTD data banks. 

To obtain the final PTD bank you then need to apply Charged Particle Tracking module that extrapolates the verteces and creates the PTD bank.


Set of algorithm parameters:

1. general section:
	"visualization": creating and saving a png image and 3D object root file for each solution of each event

	"save_sinograms": saving the images of sinograms produced during Legendre trasnform based clustering

	"force_default_sigma_r": not using the drift radii uncertainties provided by Falaise and using default value instead

	"default_sigma_r": universal default uncertainty for drift radii

	"chi_square_threshold": currently not implemented. It is intended as a threshold for acceptable linear segment fits.
	
3. clustering section:
	"clustering_max_distance": the maximum distance between two groups of tracker hits to be still considered as a part of one cluster 
		distance between tracker hits i and j = distance of their anodes:
		D(i,j) := sqrt( (x_j - x_i)^2 + (y_j - y_i)^2 )
		distance between two groups of tracker hits G1 and G2 = minimum distance between among all possible pairs of tracker cells from the two groups)
		D(G1, G2) := min( D(i,j) | i from G1, j from G2 )

	"clustering_hit_association_distance": clustering first finds a rough estimate for a line fit. Then associates tracker hits to it based on this distance threshold

	"clustering_no_iterations": number of iterations of the itterative zooming grid search algo used to find the maximum of the sinograms

	"clustering_resolution_phi": number of bins used for the grid search algo. starting range for phi is (0, phi)

	"clustering_resolution_r": number of bins used for the grid search algo. starting range for r is optimized to be minimal for each event.

	"clustering_max_initial_precision_r": If the initial range for r is low enough, number of used bins in r is reduced. In that case one bin = clustering_max_initial_precision_r

	"clustering_zoom_factor": each iteration of the search algo zooms in with this factor to achieve better precision

	"clustering_uncertainty": defines the blurring of the sinograms to compensate for errors and uncertainty caused by scattering
	
5. polyline reconstruction section:
	"polylines_max_extention_distance": trajectory can be extended by associating additional tracker hit that satisfy the "clustering_hit_association_distance". This parameter defines the maximum possible elongation of the trajectory
	
	3. a) finding sharp kinks - only vertical part of segments is changed
	
	"polylines_max_vertical_distance": maximum possible vertical distance of the two segments to be connected (computed in the point of their intersection in the horizontal projection)

	"polylines_min_tracker_hits_distance": tracker hits on both segments must be at least this close to the candidate kink to prevent the creation of distant "stranded" fake kinks

	"polylines_max_kink_angle": limits the 3D kink angle for the kink finder algo

	"polylines_min_distance_from_foil": prevents the kink to be created close to the foil
	
	3. b) finding tiny kinks - all parameters of segments can be slightly changed
				
	"polylines_max_trajectories_middlepoint_distance": maximum allowed horizontal shift of the segments to be connected into a polyline

	"polylines_max_trajectory_endpoints_distance": maximum allowed 3D distance for two endpoints of two trajectories to be connected into polyline

	"polylines_max_trajectory_connection_angle": maximum alowed 3D angle for two trajectories to be connected into one polyline trajectory
	
	
