/* The MIT License (MIT)
 * 
 * (c) Jürgen Simon 2014 (juergen.simon@uni-bonn.de)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#ifndef M3D_CLUSTERLIST_H
#define M3D_CLUSTERLIST_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/clustering/id.h>
#include <meanie3D/featurespace.h>
#include <meanie3D/utils/netcdf_utils.h>

#include <vector>
#include <map>
#include <assert.h>
#include <stdlib.h>
#include <netcdf>

namespace m3D {

    using namespace netCDF;
    using namespace std;

    template<typename T>
    class ArrayIndex;

    /** Cluster of points in feature space. A cluster is a point in feature space,
     * where multiple trajectories of original feature space points end. This end
     * point is called the cluster's 'mode'.
     *
     * At the same time, this class represents a serialized form of a list of clusters
     * to be used in further processing, such as tracking. Serialization/Deserialization
     * to and from files are done using the read/write methods.
     *
     * TODO: decouple the cluster list from the NetCDF layer
     */
    template<class T>
    class ClusterList
    {
#pragma mark -
#pragma mark Public Properties

    public:

        typedef ClusterList<T> *ptr;

        /** Actual list of clusters */
        typename Cluster<T>::list clusters;

        // General attributes

        int time_index; // contains index in 'time' dimension
        timestamp_t timestamp; // contains the variable 'time(time)' value
        
        vector<string> variables; // Variables (value range)
        vector<string> dimensions; // Dimensions 
        vector<string> dimension_variables; // Dimension variables (spatial range)
        
        string source_file; // Name of the file that the clusters were created from
        m3D::id_t highest_id; // Highest used ID
        m3D::uuid_t highest_uuid; // Highest used UUID

        // Tracking related attributes

        bool tracking_performed;
        int tracking_time_difference; // time difference between the two files in [s]
        id_set_t tracked_ids; // IDs that were continues
        id_set_t dropped_ids; // IDs that were discontinued
        id_set_t new_ids; // IDs that were added fresh
        id_map_t splits; // contains a mapping from previous id's to a list of current id's.
        id_map_t merges; // contains a mapping of current id to previous ids

        NcFile *file;
        string filename; // this is filled on read() or write()

#pragma mark -
#pragma mark Private members

    private:
        typedef map<vector<T>, typename Cluster<T>::ptr> ClusterMap;

        bool m_use_original_points_only;
        ClusterMap m_cluster_map;

#pragma mark -
#pragma mark Constructor/Destructor

    public:
        
        /**
         * Standard constructor.
         */
        ClusterList();

        /**
         * Construct an empty cluster list.
         * 
         * @param source
         * @param variables
         * @param dimensions
         * @param dimension_variables
         * @param ti time index
         * @param use_original_points_only
         */
        ClusterList(const string &source,
                const vector<string> &variables,
                const vector<string> &dimensions,
                const vector<string> &dimension_variables,
                int ti = NO_TIME,
                bool use_original_points_only = true);
        /** 
         * @param a list of clusters
         * @param all variables used to construct the points from, in the same order
         * @param the spatial dimension (2D/3D)
         * @param name of the file the clusters were created from
         * @param the command line parameters used to cluster
         */
        ClusterList(const typename Cluster<T>::list &list,
                const string &source,
                const vector<string> &vars,
                const vector<string> &dims,
                const vector<string> &dim_vars,
                int ti = NO_TIME,
                bool orig_pts = true);

        /**
         * Copy constructor
         * @param o
         */
        ClusterList(const ClusterList &o);

        /** Destructor
         */
        ~ClusterList()
        {
            if (file != NULL) {
                delete file;
                file = NULL;
            }
        };

#pragma mark -
#pragma mark Accessing the list

        /** @return number of clusters in the list
         */
        size_t size() const;

        /** @param index
         * @return cluster at index
         */
        typename Cluster<T>::ptr operator[](size_t index);

        /** Clears the cluster list. 
         * 
         * @param deletion_flag if <code>true</code>, the clusters will be
         * properly deleted. That includes their points. Defaults to
         * <code>false</code>
         */
        void clear(bool deletion_flag = false);

#pragma mark -
#pragma mark Timestamp

        /** @return time in seconds since epoch 
         */
        ::units::values::s get_time_in_seconds()
        {
            ::units::values::s seconds(this->timestamp);
            return seconds;
        }

        void set_time_in_seconds(::units::values::s seconds)
        {
            this->timestamp = seconds.get();
        }
        
    private:
        
        /**
         * Sets the value for seconds since epoch from the time dimension
         * variable at the given index.
         */
        void 
        set_time_from_time_index();

#pragma mark -
#pragma mark Clustering by Graph Theory
        
    public:   

        /** 
         * This requires that the shift has been calculated at each point
         * in feature-space and is stored in the 'shift' property of each
         * point. Aggregates clusters by taking each point, calculate the
         * shift target and find the closest point in feature-space again.
         * If two or more points have the same distance to the shifted
         * coordinate, the point with the steeper vector (=longer) is
         * chosen. If
         */
        void aggregate_cluster_graph(FeatureSpace<T> *fs,
                const WeightFunction<T> *weight_function,
                bool coalesceWithStrongestNeighbour,
                bool show_progress = true);

        /** 
         * Find the boundary points of two clusters.
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

        /** 
         * Merges the two clusters into a new cluster and removes the mergees from
         * the list of clusters, while inserting the new cluster. The mode of the merged
         * cluster is the arithmetic mean of the modes of the merged clusters.
         * @param cluster 1
         * @param cluster 2
         * @return pointer to the merged cluster
         */
        typename Cluster<T>::ptr
        merge_clusters(typename Cluster<T>::ptr c1,
                typename Cluster<T>::ptr c2);

        /** 
         * Find all directly adjacent clusters to the given cluster
         *
         * @param cluster
         * @param feature-space index for searching
         * @param resolution search radius for finding neighbours (use cluster_resolution)
         * @return list of neighbouring clusters (can be empty)
         */
        typename Cluster<T>::list
        neighbours_of(typename Cluster<T>::ptr cluster,
                ArrayIndex<T> &index);

        /** 
         * If the two clusters have neighbouring points in the grid,
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
    private:
        
        /** 
         * Estimate the tendency of the weight function at the given
         * grid point with respect to it's direct (9/16) neighbours.
         * @param point
         * @param weight function
         * @param index
         * @return positive: point is a hole, negative: point is a hill, zero: inconclusive
         */
        T weight_function_tendency(typename Point<T>::ptr p,
                const WeightFunction<T> *weight_function,
                const typename Point<T>::list &neighbours,
                ArrayIndex<T> &index);

        /**
         * 
         * @param fs
         * @param weight_function
         * @param index
         * @param coalesceWithStrongestNeighbour
         * @param show_progress
         */
        void aggregate_zeroshifts(FeatureSpace<T> *fs,
                const WeightFunction<T> *weight_function,
                ArrayIndex<T> &index,
                bool coalesceWithStrongestNeighbour,
                bool show_progress);

#pragma mark -
#pragma mark ID helpers

    public:

        /** Iterates over all clusters and sets their ID to NO_ID 
         */
        void erase_identifiers();

#pragma mark -
#pragma mark Post-Processing

        /** Iterate through the list of clusters and trow all with smaller number of
         * points than the given number out.
         * @param min number of points
         * @param show progress indicator
         */
        void apply_size_threshold(unsigned int min_cluster_size,
                const bool &show_progress = true);

#pragma mark -
#pragma mark I/O

        /** Writes out the cluster list into a NetCDF-file.
         *
         * For the format, check documentation at
         * http://git.meteo.uni-bonn.de/projects/meanie3d/wiki/Cluster_File
         * @param full path to filename, including extension '.nc'
         */
        void write(const string &path);

        /** Persists any changes. Requires that the list was read
         * from a file before.
         * @throws runtime_error if the list was not previously 
         * written or read.
         */
        void save();

        /** Persists only the top level attributes. Used to add 
         * tracking information etc. without re-writing the whole 
         * file.
         * @throws runtime_error 
         */
        void save_top_level_attributes();

        /** Static method for reading cluster lists back in.
         * @param path      : path to the cluster file
         * @param pointer to a pointer of coordinate system. 
         * If not null, this is initialized with an instance
         * of coordinate system after the reading
         */
        static
        typename ClusterList<T>::ptr
        read(const string &path, CoordinateSystem<T> **cs_ptr = NULL);

        /** Prints the cluster list out to console
         * @param include point details?
         */
        void print(bool includePoints = false);
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
        reset_clustering(FeatureSpace<T> *fs);

        /** Counts the number of points in all clusters and checks
         * if the number is equal to featureSpace->size()
         */
        void sanity_check(const FeatureSpace<T> *fs);

#pragma mark -
#pragma mark Debugging

        /** Used for analyzing the boundaries and signal correlation
         */
        void write_boundaries(const WeightFunction<T> *weight_function,
                FeatureSpace<T> *fs,
                PointIndex<T> *index,
                const vector<T> &resolution);
        
        size_t rank() {
            return dimension_variables.size() + variables.size();
        }

#if WRITE_MODES
    protected:
        vector< vector<T> > m_trajectory_endpoints;
        vector<size_t> m_trajectory_lengths;
    public:

        const vector< vector<T> > &trajectory_endpoints()
        {
            return m_trajectory_endpoints;
        }

        const vector<size_t> &trajectory_lengths()
        {
            return m_trajectory_lengths;
        }
#endif

    };
}

#endif
