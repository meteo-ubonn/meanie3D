#ifndef _M3D_VisualizationUtils_H_
#define _M3D_VisualizationUtils_H_

#include <cf-algorithms/defines.h>
#include <cf-algorithms/namespaces.h>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/locale.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <locale>
#include <fstream>
#include <string>

#include <radolan/radolan.h>

#include <meanie3D/types.h>
#include <meanie3D/tracking.h>

namespace m3D { namespace utils {

	using namespace std;
	using ::m3D::Cluster;
	using ::m3D::ClusterList;

    template <typename T>
    class VisitUtils
    {
    public:

    	/** Global variable maintaining the dimension mapping
         */
        static vector<size_t> VTK_DIMENSION_INDEXES;
        
        /** Takes the given list of  dimension names and their re-ordering.
         * Performs some sanity checks (same number of arguments, same variables)
         * and updates the field VTK_DIMENSION_INDEXES accordingly. Also updates 
         * this field in the cf-algorithms counterpart!
         *
         * @param dimension names
         * @param re-ordered dimension names
         */
        static
        void
        update_vtk_dimension_mapping(vector<string> dim_names, vector<string> vtk_dim_names);
        
    	/** Writes out the modes of the clusters. The mode is the
         * point, where the mean-shift graph ends. The resulting
         * file can be used with Visit's "Label" plot type to label
         * clusters
    	 *
    	 */
        static
    	void
        write_cluster_modes_vtk(const string &filename,
                                const typename Cluster<T>::list &list,
                                bool spatial_only=false);
        
        /** Writes out the geometrical centers of the clusters.
         * The resulting file can be used with Visit's "Label" 
         * plot type to label clusters
    	 * @param filename
         * @param cluster list
    	 */
        static
    	void
        write_geometrical_cluster_centers_vtk(const string &filename,
                                              const typename Cluster<T>::list &list);


    	/** Writes the cluster's points out as .vtk file for visit.
    	 *
    	 */
        static
    	void
        write_clusters_vtk(const string &base_name,
                           const typename Cluster<T>::list &list,
                           const vector<T> &bandwidths,
                           bool use_ids = false,
                           bool only_boundary = false);
        
        static
        void
        write_cluster_meanshift_vtk(const string &basename,
                                    const typename Cluster<T>::list &list,
                                    bool use_ids = true,
                                    bool spatial_only = true);

    	/** Writes the cluster's points out as .vtk file for visit.
    	 *
    	 */
        static
        void
        write_clusters_vtk(typename ClusterList<T>::ptr list,
                           std::string infix="_cluster_",
                           bool only_boundary = false);
        
        /** Write out the track centers */
        static
        void
        write_center_tracks_vtk(typename Tracking<T>::trackmap_t &trackmap,
                                const std::string &basename,
                                size_t spatial_dimensions,
                                bool exclude_degenerates = true);
        
        /** Writes the clusters out with the value of the weight response
         * for the whole cluster as value. 
         * @param filename (gets extended by _weight_response_[i].vtk
         * @param the cluster list
         * @param if <code>true</code>, uses the weight response of the mode.
         * If <code>false</code> the average weight response is used.
    	 */
        static
    	void
        write_cluster_weight_response_vtk(const string &base_name,
                                          const typename Cluster<T>::list &list,
                                          WeightFunction<T> *w,
                                          bool useMode=true);

    };
    
}}

#endif
