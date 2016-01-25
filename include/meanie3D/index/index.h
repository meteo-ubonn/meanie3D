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

#ifndef M3D_FEATURESPACE_INDEX_H
#define M3D_FEATURESPACE_INDEX_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <meanie3D/featurespace/point.h>
#include <meanie3D/index/search_parameters.h>

#include <vector>

namespace m3D {

    template <typename T>
    class FeatureSpace;

    /** Abstract base class. This interface abstracts the various implementations
     * used for indexing the points, such as kdtree.c, flann etc.
     * 
     * TODO: remove the feature-space entanglement
     *
     * @abstract
     */
    template <typename T>
    class PointIndex
    {
    protected:

        // Member Variables

        typename Point<T>::list *m_points;

        /** A reference to the feature-space.
         */
        FeatureSpace<T> *m_fs;

        /** Stores the indicees of the variables used for
         *  building the index.
         */
        vector<size_t> m_index_variable_indexes;

        /** Debugging method. Writes out the search window and found points to files.
         */
        void
        write_search(const vector<T>& x,
                const vector<T> &ranges,
                typename Point<T>::list *result);

    public:

#pragma mark -
#pragma mark Statics

        /** Controls the 'logging' of index searches
         */
        static bool write_index_searches;

#pragma mark -
#pragma mark Definitions

        /** Enumerates the various concrete implementations available
         */
        typedef enum
        {
            /** Standard linear search. Exact but very slow
             */
            IndexTypeLinear,

            /** Indexing via kdtree.c
             */
            IndexTypeKDTree,

            /** Approximate KNN using FLANN. Many options available
             */
            IndexTypeFLANN,

            /** Special quick lookup index for rectilinear coordinate systems 
             */
            IndexTypeRectilinearGrid

        } IndexType;

#pragma mark -
#pragma mark Defaults

        /** Default
         */
        static const IndexType DefaultIndexType = IndexTypeFLANN;

#pragma mark -
#pragma mark Public Abstract Methods

        /** Searches the index according to the given search parameters.
         * @abstract
         * @param x
         * @param search parameters
         * @param Some search algorithms also return the distances. Highly implementation specific!
         * @return vector of the n closest feature space points.
         */
        virtual
        typename Point<T>::list *
        search(const vector<T> &x, const SearchParameters *params, vector<T> *distances = NULL) = 0;

        /** Add a new point to the index. If the point already exists, it is
         * replaced with the new point
         * @param feature-space point
         */
        virtual
        void
        add_point(typename Point<T>::ptr p) = 0;

        /** Remove a point from the index. 
         * @param feature-space point
         */
        virtual
        void
        remove_point(typename Point<T>::ptr p) = 0;


#pragma mark -
#pragma mark Public Methods

        /** Performs a 1-nn search and checks, if the result matches the
         * given point. If necessary, the indexed components are extracted
         * from point->values first.
         *
         * @param feature-space point
         * @return true or false
         */
        bool
        has_point(typename Point<T>::ptr p) const;

        /** Performs a 1-nn search and checks, if the result matches the
         * given point. The value is taken as is, if the points were indexed
         * using a subset of indices, the caller is responsible for the
         * correct mapping.
         *
         * @param vector
         * @return true or false
         */
        bool
        has_value(vector<T> &value);

        /** Picks the component defined by the index variables 
         * from the given point's values.
         * @param feature-space point
         * @return vector for search
         */
        vector<T>
        indexed_components(typename Point<T>::ptr p);

        /** Accesor
         * @return feature space
         */
        inline
        const FeatureSpace<T> *feature_space()
        {
            return m_fs;
        };

        inline
        const vector<size_t> index_variable_indexes()
        {
            return m_index_variable_indexes;
        }

        /** @return dimensionality of index
         */
        inline
        const size_t dimension()
        {
            return m_index_variable_indexes.size();
        }

        /** @return size of index (number of indexed points)
         */
        inline
        size_t size()
        {
            return this->m_points->size();
        }

#pragma mark -
#pragma mark Destructor

        /** Destructor
         */
        virtual ~PointIndex()
        {
        };

        /** Factory method
         */
#pragma mark -
#pragma mark Factories

        /** Creates an index for the given points, 
         * restricting the
         *
         * @param point list
         * @param point dimension
         * @param index type
         */
        static PointIndex<T> *
        create(typename Point<T>::list *points,
                size_t dimension,
                IndexType index_type = DefaultIndexType);

        /** Creates an index for the given points by using a subset of variables
         * as indicated by their indices in the point->values vector.
         * @param point list
         * @param indices
         * @param index type
         */
        static PointIndex<T> *
        create(typename Point<T>::list *points,
                const vector<size_t> &indexes,
                IndexType index_type = DefaultIndexType);


#pragma mark -
#pragma mark Protected Constructors

    protected:

        // Constructors

        /** Constructor
         * @param pointer to feature space
         * @param dimension of the points 
         */
        PointIndex(typename Point<T>::list *points, size_t dimension)
        : m_points(points)
        , m_fs(NULL)
        {
            this->m_index_variable_indexes = vector<size_t>(dimension);

            for (size_t var_index = 0; var_index < dimension; var_index++) {
                this->m_index_variable_indexes[var_index] = var_index;
            }
        };

        /** Constructor
         * @param pointer to feature space
         * @param dimension of the points
         * @param indices to use for building index
         */
        PointIndex(typename Point<T>::list *points,
                const vector<size_t> &indexes)
        : m_points(points)
        , m_fs(NULL)
        , m_index_variable_indexes(indexes)
        {
        };

        /** Constructor
         * @param pointer to feature space
         */
        PointIndex(FeatureSpace<T> *fs)
        : m_points(&fs->points)
        , m_fs(fs)
        {
            this->retrieve_variables_indexes();
        };

        /** Copy constructor
         */
        PointIndex(const PointIndex<T> &o)
        : m_points(o.m_points)
        , m_fs(o.m_fs)
        , m_index_variable_indexes(o.index_variable_indexes())
        {
        };

        /** This method turns the given list of variables into a list of indexes to be used
         * TODO: remove this when disentangling index and featurespace
         */
        void
        retrieve_variables_indexes();

        /** Triggers the construction of an index of the underlying
         * featurespace.
         * @abstract
         * @param ranges many KD-Tree implementations work with one radius parameter only, which will require pre-whitening
         *               of the data, which requires the 'ranges' parameter.  Implementations who don't have this problem,
         *               can ignore this parameter.
         */
        virtual
        void
        build_index(const vector<T> &ranges) = 0;
    };
}

#endif
