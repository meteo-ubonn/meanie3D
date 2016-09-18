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


#ifndef M3D_FEATURESPACE_H
#define M3D_FEATURESPACE_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <meanie3D/array/multiarray.h>
#include <meanie3D/featurespace/coordinate_system.h>
#include <meanie3D/featurespace/point.h>
#include <meanie3D/featurespace/data_store.h>

#include <boost/progress.hpp>

#include <map>
#include <vector>

namespace m3D {

    // Forward declarations

    template<typename T>
    class PointIndex;

    using std::map;
    using std::vector;

    /** This represents a complete feature space. The template parameter is the
     * data type used (eg. float/double, etc.)
     *
     * TODO: extend the construction of the feature space so, that the validity of
     * each point is ascertained by evaluating dimension variables and other variables,
     * but the actual feature space can be chosen from a subset of those again.
     */
    template<typename T>
    class FeatureSpace
    {
    public:

#pragma mark -
#pragma mark public constants

        static map<int, double> *NO_THRESHOLDS;

#pragma mark -
#pragma mark General Typedefs

        typedef FeatureSpace<T> *ptr;
        typedef map<typename CoordinateSystem<T>::Coordinate, bool> ValidityMap;
        typedef map<typename CoordinateSystem<T>::Coordinate, typename Point<T>::ptr> FeatureSpaceMap;
        typedef typename FeatureSpaceMap::iterator FSMapIterator;
        typedef vector<vector<T> > Trajectory;
        typedef typename Trajectory::iterator TrajectoryIterator;
        typedef MultiArray <T> multi_array_t;


#pragma mark -
#pragma mark Private

    private:

#pragma mark -
#pragma mark Private Member Variables

#pragma mark -
#pragma mark Private Methods

        /** Reads the data from a NetCDF file variable. in constructing the feature space.
         * Calls build_recursive(...).
         */
        void build();

#pragma mark -
#pragma mark Protected

    protected:

#pragma mark -
#pragma mark Protected Member Variables

        /** The data store contains the data used to populate
         * the value range of the feature space 
         */
        const DataStore <T> *m_data_store;

        /** Progress indicator for construction 
         */
        boost::progress_display *m_progress_bar;

        /** Parameter used for construction. Points in feature-space,
         * where one variable's value is less than the given threshold,
         * will not be added.
         */
        map<int, double> m_lower_thresholds;

        /** Parameter used for construction. Points in feature-space,
         * where one variable's value is more than the given threshold,
         * will not be added.
         */
        map<int, double> m_upper_thresholds;

        /** Parameter used for construction. Points in feature-space,
         * where one variable's value is out of valid range are usually
         * ignored. By providing fill values, it is possible to replace
         * the missing value with the value provided, thus allowing
         * use of that point in spite of missing values.
         */
        map<int, double> m_replacement_values;

        // bool flag array, indicating where the original data sets
        // had '_fillValue' or values outside of valid_range.
        // Note: only applicaple for gridded data!
        MultiArray<bool> *m_off_limits;

        // Factual limits

        map<size_t, T> m_min;
        map<size_t, T> m_max;

#pragma mark-
#pragma mark Public Member Variables

    public:

        typename Point<T>::list points;
        size_t dimension;
        const CoordinateSystem <T> *coordinate_system;

#if WRITE_MEANSHIFT_WEIGHTS
        vector< vector<T> > weight_sample_points;
#endif


#pragma mark -
#pragma mark Protected Methods

        /** Construct a map of Point objects, representing the feature space.
         * The variables used to do so are the dimension variables plus the
         * additional variables handed in (at least one)
         * 
         * @param show_progress
         */
        void construct_featurespace(bool show_progress = true);

#pragma mark -
#pragma mark Public

    public:

#pragma mark -
#pragma mark Consctructor/Destructor

        /** Constructs a feature-space from the variables in the
         * given data store.
         * 
         * @param coordinate_system
         * @param dataStore
         * @param lower_thresholds
         * @param upper_thresholds
         * @param replacement_values
         * @param show_progress
         */
        FeatureSpace(const CoordinateSystem <T> *coordinate_system,
                     const DataStore <T> *dataStore,
                     const map<int, double> &lower_thresholds,
                     const map<int, double> &upper_thresholds,
                     const map<int, double> &replacement_values,
                     const bool &show_progress = true);

        // Making copies

        /** Copy constructor
         * @param other
         * @param without_points if true, the points are not copied
         */
        FeatureSpace<T>(const FeatureSpace<T> &other,
                        bool with_points = true);

        /** Copy constructor
         * @param pointer to other 
         * @param without_points if true, the points are not copied
         */
        FeatureSpace<T>(const FeatureSpace<T> *other,
                        bool with_points = true);

        /** Copy operator
         * @param other
         * @return pointer to a copy
         */
        FeatureSpace<T> operator=(const FeatureSpace<T> &other);

        /** Destructor
         */
        ~FeatureSpace();

#pragma mark -
#pragma mark Accessors

        /** @return the data store used to construct this featurespace
         */
        const DataStore <T> *data_store() const {
            return m_data_store;
        }

        /** @return number of points in featurespace
         */
        const size_t size() const {
            return points.size();
        }

        /** @return a pointer to the points array 
         */
        typename Point<T>::list *get_points() {
            return &(this->points);
        }

        /** @return rank of featurespace, which is rank of the
         * spatial range plus rank of the value range.
         */
        const size_t rank() const {
            return coordinate_system->rank() + data_store()->rank();
        }

        /** @return the number of components in the spatial
         * range of the featurespace
         */
        inline const size_t spatial_rank() const {
            return coordinate_system->rank();
        }

        /** @return the number of components in the value
         * range of the featurespace. Note: when value range 
         * is omitted, this is always 0.
         */
        inline const size_t value_rank() const {
            return m_data_store->rank();
        }

        /** Lower thresholds used in the construction of the featurespace. 
         * Only variables from the value range are used.
         */
        const map<int, double> &lower_thresholds() const {
            return m_lower_thresholds;
        }

        /** Upper thresholds used in the construction of the featurespace.
         * Only variables from the value range are used.
         */
        const map<int, double> &upper_thresholds() const {
            return m_upper_thresholds;
        }

        /** A map of bool values, which indicates where in the construction
         * points were encountered, that are 'off limits' in the sense of 
         * not valid.
         */
        const MultiArray<bool> *off_limits() const {
            return m_off_limits;
        }

        /** Contains the inf of all points in the featurespace.
         */
        const map<size_t, T> &min() const {
            return m_min;
        }

        /** Contains the sup of all points in the featurespace
         */
        const map<size_t, T> &max() const {
            return m_max;
        }


#pragma mark -
#pragma mark Utilities

        /** print details, including points, out to console
         */
        void print() const;

        /** Clears the points array and releases the memory.
         */
        void clear();

        /** @return vector containing the indexes of the
         * spatial range component
         */
        vector<size_t> spatial_range_indexes() {
            vector<size_t> result(this->coordinate_system->rank());
            for (size_t i = 0; i < this->coordinate_system->rank(); ++i)
                result[i] = i;
            return result;
        }

        /** @return vector containing the indexes of the
         * value range component. If omit_value_range is set, 
         * the resulting vector is empty
         */
        vector<size_t> range_indexes() {
            vector<size_t> result;

            if (this->m_omit_value_range) {
                result = vector<size_t>(this->data_store()->rank());
                for (size_t i = 0; i < this->data_store()->rank(); ++i)
                    result[i] = this->coordinate_system->rank() + i;
            }

            return result;
        }

        /** Returns a subvector, containing the spatial range component 
         * of the given feature space value only.
         *
         * @param point in feature space
         * @return it's spatial component
         */
        typename CoordinateSystem<T>::Coordinate
        spatial_component(const vector<T> &value) const;

        /**  
         * @param p
         * @param index
         * @return the value at the n-th component of the spatial
         */
        T
        get_spatial_component_at(typename Point<T>::ptr p, size_t index) const;

        /** 
         * @param point
         * @param resolution
         * @return the value at the n-th component of the value
         */
        T
        get_value_component_at(typename Point<T>::ptr p, size_t index) const;

        /** Changes the given vector so, that it's components are multiples
         * of the grid resolution.
         *
         * @param point to round
         * @param resolution to round to (for example (5km,5km,10dbZ))
         */
        void
        round_to_resolution(vector<T> &point,
                            const vector<T> &resolution) const;

        /** Count the number of points in the featurespace
         * that were part of the original make up (and not
         * filtered out yet)
         * @return number of points with 'isOriginalPoint' = true
         */
        size_t
        count_original_points() const;

        /** Debug tool. Checks if any of the points are out of the 
         * ordinary. Throws an exception if so.
         */
        void
        sanity_check();

    };
}

#endif
