#ifndef _M3D_ClusterList_H_
#define _M3D_ClusterList_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>
#include <map>
#include <assert.h>
#include <boost/shared_ptr.hpp>
#include <stdlib.h>
#include <netcdf>

#include <cf-algorithms/cf-algorithms.h>

#include <meanie3D/types/cluster.h>
#include <meanie3D/types/point.h>

namespace m3D {
    
    using namespace netCDF;
	using namespace std;
	using ::cfa::meanshift::Point;
    using ::cfa::meanshift::WeightFunction;

    /** Cluster of points in feature space. A cluster is a point in feature space,
     * where multiple trajectories of original feature space points end. This end
     * point is called the cluster's 'mode'.
     *
     * At the same time, this class represents a serialized form of a list of clusters
     * to be used in further processing, such as tracking. Serialization/Deserialization
     * to and from files are done using the read/write methods.
     */
    template <class T>
    class ClusterList
    {
        
    private:
        
        typedef map< vector<T>, typename Cluster<T>::ptr > ClusterMap;
        
        ClusterMap      m_cluster_map;
        
        
        /** Standard constructor is private
         */
        ClusterList()
        {};
        
        /** Recursive helper method fort find_neighbours
         *
         */
        void
        find_neighbours(typename CoordinateSystem<T>::GridPoint &gridpoint,
                        size_t dim_index,
                        ArrayIndex<T> &index,
                        typename Point<T>::list &list);
        
        typename Point<T>::list
        find_neighbours(const typename CoordinateSystem<T>::GridPoint &gridpoint,
                        ArrayIndex<T> &index);
        
        /** Sanity check / consistency check
         */
        void
        check_clusters(ArrayIndex<T> &index);
        

        void
        aggregate_zeroshifts(FeatureSpace<T> *fs,
                             ArrayIndex<T> &index,
                             bool show_progress);
        
    public:
        
#pragma mark -
#pragma mark Public typedefs
        
        typedef ClusterList<T> *    ptr;
        
#pragma mark -
#pragma mark Public Properties

        // meta-info
        
        NcFile          *ncFile;
        
        vector<NcVar>   feature_variables;            // all variables, including dimension variables
        
        vector<NcDim>   dimensions;                   // all dimensions that were in the featurespace
        
        string          source_file;                  // Name of the file that the clusters were created from
        
        // payload

        typename Cluster<T>::list   clusters;
        
        // tracking help
        
        bool            tracking_performed;
        
        vector<size_t>  tracked_ids;
        
        vector<size_t>  dropped_ids;
        
        vector<size_t>  new_ids;
        
        // control
        
        bool            m_use_original_points_only;
        

#pragma mark -
#pragma mark Constructor/Destructor

        /** @constructor
         * @param all variables used to contruct the points from, in the same order
         * @param all dimensions that were in the featurespace
         * @param name of the file the clusters were created from
         * @param the command line parameters used to cluster
         * @param when searching for predecessors in the meanshift vector graph,
         *          there is a maximum distance in gridpoints, after which the 
         *          neighbourhood search will not accept points. This parameter
         *          sets that distance in #grid-points.
         * @param In the process of linking up meanshift vectors to a graph, the
         *          vicinity of vector endpoints is searched for other points to
         *          link up to. The maximum distance of that search (in grid points)
         *          is determined here.
         * @param When adding points to clusters, this switch can be used to control
         *          if only points from the original data set are added or also points
         *          that resulted from filtering (for example scale-space)
         */
        ClusterList(const vector<NcVar> &variables,
                    const vector<NcDim> &dims,
                    const string& sourcefile,
                    bool use_original_points_only=true)
        : ncFile(NULL)
        , feature_variables(variables)
        , dimensions(dims)
        , source_file(sourcefile)
        , tracking_performed(false)
        , m_use_original_points_only(use_original_points_only)
        {};

        /** @constructor
         * @param a list of clusters
         * @param all variables used to contruct the points from, in the same order
         * @param the spatial dimension (2D/3D)
         * @param name of the file the clusters were created from
         * @param the command line parameters used to cluster
         */
        ClusterList(const typename Cluster<T>::list &list,
                    const vector<NcDim> &dims,
                    const vector<NcVar> &variables,
                    const string& sourcefile )
        : ncFile(NULL)
        , feature_variables(variables)
        , dimensions(dims)
        , source_file(sourcefile)
        , clusters(list)
        , tracking_performed(false)
        , m_use_original_points_only(true)
        {};

        /** Destructor
         */
        ~ClusterList()
        {
            close_ncFile();
        };
        
        void close_ncFile()
        {
            if (ncFile!=NULL)
            {
                delete ncFile;
                ncFile=NULL;
            }
        };
        
#pragma mark -
#pragma mark Accessing the list
        
        /** @return number of clusters in the list
         */
        size_t size();
        
        /** @param index
         * @return cluster at index
         */
        typename Cluster<T>::ptr operator[] (size_t index);
        
#pragma mark -
#pragma mark Clustering by Graph Theory

        // TODO: make the API more consistent by ordering the parameters equally
        // TODO: assign access classifiers (private/protected/public) where needed
        
        /** This requires that the shift has been calculated at each point
         * in feature-space and is stored in the 'shift' property of each
         * point. Aggregates clusters by taking each point, calculate the
         * shift target and find the closest point in feature-space again.
         * If two or more points have the same distance to the shifted
         * coordinate, the point with the steeper vector (=longer) is
         * chosen. If
         */
        void aggregate_cluster_graph(FeatureSpace<T> *fs,
                                     const WeightFunction<T> *weight_function,
                                     bool show_progress=true);
        
        /** Find the boundary points of two clusters.
         * @param cluster 1
         * @param cluster 2
         * @return vector containing the boundary points
         * @param feature-space index
         * @param search range (use cluster resolution)
         */
        typename Point<T>::list
        get_boundary_points(typename Cluster<T>::ptr c1,
                            typename Cluster<T>::ptr c2,
                            ArrayIndex<T> &index);
        
        /** Merges the two clusters into a new cluster and removes the mergees from
         * the list of clusters, while inserting the new cluster. The mode of the merged
         * cluster is the arithmetic mean of the modes of the merged clusters.
         * @param cluster 1
         * @param cluster 2
         * @return pointer to the merged cluster
         */
        typename Cluster<T>::ptr
        merge_clusters(typename Cluster<T>::ptr c1,
                       typename Cluster<T>::ptr c2);

        /** Find all directly adjacent clusters to the given cluster
         *
         * @param cluster
         * @param feature-space index for searching
         * @param resolution search radius for finding neighbours (use cluster_resolution)
         * @return list of neighbouring clusters (can be empty)
         */
        typename Cluster<T>::list
        neighbours_of(typename Cluster<T>::ptr cluster,
                      ArrayIndex<T> &index);

        /** If the two clusters have neighbouring points in the grid,
         * they are counted as neighbours.
         * @param cluster one
         * @param cluster two
         * @param index used for searching boundary points
         * @return <code>true</code> if they are neighbours <code>false</code>
         */
        bool
        are_neighbours(const Cluster<T> *c1,
                       const Cluster<T> *c2,
                       ArrayIndex<T> &index);
#pragma mark -
#pragma mark ID helpers
        
        /** Iterates over all clusters and sets their ID to NO_ID 
         */
        void erase_identifiers();
        
        /** Re-tags sequentially, starting with 0 
         */
        void retag_identifiers();
        
#pragma mark -
#pragma mark Post-Processing
        
        /** Iterate through the list of clusters and trow all with smaller number of
         * points than the given number out.
         * @param min number of points
         * @param show progress indicator
         */
        void apply_size_threshold(unsigned int min_cluster_size,
                                  const bool& show_progress = true);
        
#pragma mark -
#pragma mark I/O
        
        /** Writes out the cluster list into a NetCDF-file.
         * For the format, check documentation at
         * http://git.meteo.uni-bonn.de/projects/meanie3d/wiki/Cluster_File
         * @param full path to filename, including extension '.nc'
         * @param feature space
         * @param parameters used in the run
         */
        void write( const string& path );
        
        /** Static method for reading cluster lists back in.
         * @param path      : path to the cluster file
         * @param list      : contains the clusters after reading
         * @param source    : contains the source file attribute's value after reading
         * @param parameters: contains the run's parameter list after reading
         * @param var_names : contains the list of variables used after reading
         */
        static
        typename ClusterList<T>::ptr
        read(const string& path );
        
        /** Prints the cluster list out to console
         */
        void print();
        /** Counts the number of points in all the clusters. Must be
         * equal to the number of points in the feature space.
         */
        
#pragma mark -
#pragma mark Miscellaneous
        
        /** Sets the cluster property of all points of the given
         * featurespace to NULL
         * @param featurespace
         * TODO: find a better place for this!
         */
        static
        void
        reset_clustering( FeatureSpace<T> *fs );
        
        /** Counts the number of points in all clusters and checks
         * if the number is equal to featureSpace->size()
         */
        void sanity_check( const FeatureSpace<T> *fs );
        
#pragma mark -
#pragma mark Debugging

        /** Used for analyzing the boundaries and signal correlation
         */
        void write_boundaries(const WeightFunction<T> *weight_function,
                              FeatureSpace<T> *fs,
                              PointIndex<T> *index,
                              const vector<T> &resolution );
#if WRITE_MODES
    protected:
        vector< vector<T> >         m_trajectory_endpoints;
        vector<size_t>              m_trajectory_lengths;
    public:
        const vector< vector<T> > &trajectory_endpoints() { return m_trajectory_endpoints; }
        const vector<size_t> &trajectory_lengths() { return m_trajectory_lengths; }
#endif

    };
};
    
#endif
