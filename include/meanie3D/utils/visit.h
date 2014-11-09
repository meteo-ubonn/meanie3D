#ifndef M3D_VISUALIZATIONUTILS_H
#define M3D_VISUALIZATIONUTILS_H

#if WITH_VTK

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/tracking.h>
#include <meanie3D/clustering.h>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/locale.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <vtkSmartPointer.h>
#include <vtkLookupTable.h>
#include <vtkRectilinearGrid.h>

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <locale>
#include <fstream>
#include <string>

#include <radolan/radolan.h>

namespace m3D { 
    
    template <typename T>
    class ConradCluster;
    
    namespace utils {

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

        static void
        get_vtk_image_dimensions(const CoordinateSystem<T> *cs,
                int &nx,int &ny,int &nz);

        static void
        get_vtk_coords(const vector<T>& coords, T &x,T &y,T &z);

        static void
        get_vtk_gridpoint(const vector<int>& gp, int &gx,int &gy,int &gz);

        static
        vtkRectilinearGrid *
        allocate_vtk_rectilinear_grid(const CoordinateSystem<T> *cs, 
                vector<vtkDoubleArray *> &coords);

        static
        vtkIdType
        to_single_index(const CoordinateSystem<T> *cs, 
                int nx, int ny, int nz, 
                int gx, int gy, int gz);

        static size_t index_of(size_t dim)
        {
            return VTK_DIMENSION_INDEXES.empty() || dim > (VTK_DIMENSION_INDEXES.size()-1)
                ? dim
                : VTK_DIMENSION_INDEXES[dim];
        }

        /** Helper function for drawing an ellipse (representing bandwidth and somesuch)
         * in a .curve file format.
         * @param filename
         * @param half axis data
         * @param number of segments in the curve (default 1000)
         * @param pointer to origin vector, defaults to NULL = (0,0)
         */
        static
        void
        write_ellipsis_2d(const string &filename,
                          const vector<T> &half_axis,
                          size_t number_of_segements = 1000,
                          const vector<T> *origin = NULL );

        /** Helper function for drawing an ellipse (representing bandwidth and somesuch)
         * in a .curve file format.
         * @param filename
         * @param half axis data
         * @param number of segments in the curve (default 1000)
         * @param pointer to origin vector, defaults to NULL = (0,0,0)
         */
        static
        void
        write_ellipsis_3d(const string &filename,
                          const vector<T> &half_axis,
                          size_t number_of_segements = 250,
                          const vector<T> *origin = NULL );

        static
        void
        write_featurespace_vtk(const string &filename,
                               FeatureSpace<T> *fs,
                               string variable_name="featurespace");

        /** Writes out the given variables of the featurespace individually
         * @param output filename (without extension)
         * @param the featurespace
         * @param list of variables to include
         * @param if <code>true</code>, the legacy ASCII format is used. 
         *        but if <code>false</code>, xml is used.
         */
        static
        void
        write_featurespace_variables_vtk(const string &filename,
                                         FeatureSpace<T> *fs,
                                         const vector<string> &feature_variables,
                                         const vector<string> &vtk_variables,
                                         bool write_legacy = false);

        /** Writes out the given variables of the featurespace individually
         * as VTK image data files.
         * @param output filename
         * @param the featurespace
         * @param list of variables to include
         */
        static
        void
        write_featurespace_variables_vti(const string &filename,
                                         FeatureSpace<T> *fs,
                                         const vector<NcVar> &vtk_variables);

        static
        void
        write_vectors_vtk(const string& filename,
                          const vector< vector<T> > &origins,
                          const vector< vector<T> > &vectors,
                          string var_name="vectors" );

        static
        void
        write_matrix_vtk(const string &filename,
                         const boost::numeric::ublas::matrix<T> &matrix,
                         string var_name = "matrix" );

        static
        void
        write_ublas_matrix_vtk(const string& filename,
                               const vector< vector<T> > &origins,
                               const vector< vector<T> > &vectors );

        static
        void
        write_pointlist_vtk(const string &filename,
                            typename Point<T>::list *list,
                            size_t dim,
                            string var_name="points" );

        static
        void
        write_pointlist_all_vars_vtk(const string &filename,
                                     typename Point<T>::list *list,
                                     const vector<string> &var_names);

        static
        void
        write_pointlist_vtk(const string &filename,
                            vector< vector<T> > *list,
                            string var_name="points" );

        static
        void
        write_modes_vtk(const string &filename,
                        const vector< vector<T> > &list,
                        const vector<size_t> &trajectory_lenghts,
                        string var_name="mode" );

        static
        void
        write_radolan_vtk(const string &filename,
                          const string &outfile,
                          Radolan::RDDataType* threshold = NULL );

        static
        void
        write_weights(const string &filename,
                      const string &var_name,
                      typename Point<T>::list *list,
                      const vector<T> &weights,
                      bool restrict_to_2D = false );

        static
        void
        write_shift_vectors(const string &filename,
                            FeatureSpace<T> *fs,
                            bool spatial_only=true );

        static
        void
        write_weight_function_response(const string& filename,
                                       FeatureSpace<T> *fs,
                                       WeightFunction<T> *weight_function,
                                       bool write_legacy = false);

        static
        void
        write_multiarray_vtk(const std::string &filename,
                             const std::string &variable_name,
                             const CoordinateSystem<T> *cs,
                             MultiArray<T> *array,
                             bool write_legacy_format=false);

        static
        void
        write_multiarray_vtk(const std::string &filename,
                             const std::string &variable_name,
                             const CoordinateSystem<T> *cs,
                             MultiArray<bool> *array);

        /**
         */
        static vtkSmartPointer<vtkLookupTable> cluster_lookup_table();

        /** Takes the given list of  dimension names and their re-ordering.
         * Performs some sanity checks (same number of arguments, same variables)
         * and updates the field VTK_DIMENSION_INDEXES accordingly. Also updates 
         * this field in the meanie3D counterpart!
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
                                              const typename Cluster<T>::list &list,
                                              bool at_max_height=true);


//    	/** Writes the cluster's points out as .vtk file for visit.
//         *
//    	 * @deprecated
//    	 */
//        static
//    	void
//        write_clusters_vtk(const string &base_name,
//                           const typename Cluster<T>::list &list,
//                           const vector<T> &bandwidths,
//                           bool use_ids = false,
//                           bool only_boundary = false);

        static
        void
        write_cluster_meanshift_vtk(const string &basename,
                                    const typename Cluster<T>::list &list,
                                    bool use_ids = true,
                                    bool spatial_only = true);

        /** Writes the cluster's points out as unstructured grid
         * of 2D quadrilaterals or 3D orthogonal parallelepiped
         * @param cluster list
         * @param coordinate system (required for resolution)
         * @param infix for cluster filenames
         * @param only_bounday: If <code>true</code> only points marked
         * as boundary points are written out. If <code>false</code> all
         * points are written out (default).
         * @param write_ascii If <code>true</code>, the files is written as
         * legacy ascii (expensive). If <code>false</code> (default) the more
         * efficient xml format is used.
         */
        static
        void
        write_clusters_vtu(const ClusterList<T> *list,
                           const CoordinateSystem<T> *cs,
                           const string &base_name,
                           unsigned int max_colors=5,
                           bool use_ids = true,
                           bool only_boundary = false,
                           bool write_xml = false);

        /** Writes the cluster's points out as rectilinear grid
         * @param cluster list
         * @param coordinate system (required for resolution)
         * @param infix for cluster filenames
         * @param only_bounday: If <code>true</code> only points marked
         * as boundary points are written out. If <code>false</code> all
         * points are written out (default).
         * @param write_ascii If <code>true</code>, the files is written as
         * legacy ascii (expensive). If <code>false</code> (default) the more
         * efficient xml format is used.
         */
        static
        void
        write_clusters_vtr(const ClusterList<T> *list,
                           const CoordinateSystem<T> *cs,
                           const string &base_name,
                           bool use_ids=true,
                           bool only_boundary=false,
                           bool write_ascii=false);

        /** Write out the track centers 
         */
        static
        void
        write_center_tracks_vtk(typename Track<T>::trackmap &trackmap,
                                const std::string &basename,
                                size_t spatial_dimensions,
                                bool exclude_degenerates = true);

        /** Write out track centers from CONRAD clusters
         */
        static
        void
        write_center_tracks_vtk(typename ConradCluster<T>::trackmap_t &track_map,
                                const std::string &basename,
                                bool exclude_degenerates);

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

#endif
