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

#include <radolan/radolan.h>

#include <meanie3D/types.h>

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

    	/**
    	 *
    	 */
        static
    	void write_cluster_modes_vtk(const string &filename,
    								 const typename Cluster<T>::list &list,
                                     bool spatial_only=false);

    	/**
    	 *
    	 */
        static
    	void write_clusters_vtk(const string &base_name,
    							const typename Cluster<T>::list &list,
    							const vector<T> &bandwidths,
    							bool use_ids = false );
        /**
         */
        static
        void
        write_clusters_vtk( typename ClusterList<T>::ptr list, std::string infix="_cluster_" );

    };
    
}};

#endif
