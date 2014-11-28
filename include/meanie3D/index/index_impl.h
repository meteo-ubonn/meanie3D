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

#ifndef M3D_FEATURESPACE_INDEXIMPL_H
#define M3D_FEATURESPACE_INDEXIMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/index.h>
#include <meanie3D/utils/visit.h>

#include <exception>
#include <algorithm>

#include <boost/lexical_cast.hpp>

namespace m3D { 

#if WITH_VTK
    using utils::VisitUtils;
#endif

    template <typename T>
    bool PointIndex<T>::write_index_searches = false;

    template <typename T>
    PointIndex<T> *
    PointIndex<T>::create( typename Point<T>::list *points, size_t dimension, IndexType index_type )
    {
        PointIndex<T> *instance = NULL;

        switch ( index_type )
        {
            case IndexTypeLinear:
                instance = new LinearIndex<T>( points, dimension );
                break;

            case IndexTypeFLANN:
                instance = new FLANNIndex<T>( points, dimension );
                break;

            case IndexTypeKDTree:
                instance = new KDTreeIndex<T>( points, dimension );
                break;

            case IndexTypeRectilinearGrid:
                instance = new RectilinearGridIndex<T>( points, dimension );
                break;
        }

        return instance;
    }

    template <typename T>
    PointIndex<T> *
    PointIndex<T>::create(typename Point<T>::list *points,
                          const vector<size_t> &indexes,
                          IndexType index_type)
    {
        PointIndex<T> *instance = NULL;

        switch ( index_type )
        {
            case IndexTypeLinear:
                instance = new LinearIndex<T>( points, indexes );
                break;

            case IndexTypeFLANN:
                instance = new FLANNIndex<T>( points, indexes );
                break;

            case IndexTypeKDTree:
                instance = new KDTreeIndex<T>( points, indexes );
                break;

            case IndexTypeRectilinearGrid:
                instance = new RectilinearGridIndex<T>( points, indexes );
                break;
        }

        return instance;
    }

    template <typename T>
    void
    PointIndex<T>::write_search( const vector<T>& x, const vector<T> &ranges, typename Point<T>::list *result )
    {
        static size_t search_count = 0;

        size_t dim = x.size();

        // Write search window

        // TODO: find some generic way of writing out search windows

        string extension = ( dim == 2 ? ".curve" : ".3D" );

        string fn = "search_window_" + boost::lexical_cast<string>(search_count) + extension;

        if ( dim == 3 )
        {
#if WITH_VTK
            VisitUtils<T>::write_ellipsis_3d( fn, ranges, 250, &x );
#endif
        }
        else if ( dim == 2 )
        {
#if WITH_VTK
            VisitUtils<T>::write_ellipsis_2d( fn, ranges, 250, &x );
#endif
        }

#if WITH_VTK
        fn = "search_result_" + boost::lexical_cast<string>(search_count) + ".vtk";
        VisitUtils<T>::write_pointlist_vtk( fn, result, dim, "search_result" );
#endif
        search_count++;
    }

    template <typename T>
    void
    PointIndex<T>::retrieve_variables_indexes()
    {
        // Default behaviour: if no indexes are set expressively,
        // use the whole range in order

        if (m_index_variable_indexes.empty())
        {
            for ( size_t i=0; i < m_fs->rank(); i++)
            {
                m_index_variable_indexes.push_back(i);
            }
        }
    }

    template <typename T>
    vector<T>
    PointIndex<T>::indexed_components( typename Point<T>::ptr p )
    {
        vector<T> result( m_index_variable_indexes.size() );

        for ( size_t i=0; i < m_index_variable_indexes.size(); i++ )
        {
            result[i] = p->values[ m_index_variable_indexes[i] ];
        }

        return result;
    }

    template <typename T>
    bool
    PointIndex<T>::has_point( typename Point<T>::ptr p ) const
    {
        bool have_point = false;

        KNNSearchParams<T> knn(1);

        vector<T> x = this->indexed_components( p );

        typename Point<T>::list *result = this->search( x, &knn );

        if ( result->size() > 0 )
        {
            have_point = ( result->front() == p );
        }

        return have_point;
    }

    template <typename T>
    bool
    PointIndex<T>::has_value( vector<T> &x )
    {
        bool have_point = false;

        KNNSearchParams<T> knn(1);

        typename Point<T>::list *result = this->search( x, &knn );

        if ( result->size() > 0 )
        {
            vector<T> xr = this->indexed_components( result->front() );

            have_point = ( xr == x );
        }

        return have_point;
    }
}
    
#endif
