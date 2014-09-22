#ifndef _M3D_Operation_Cluster_Task_H_
#define _M3D_Operation_Cluster_Task_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <stdlib.h>
#include <cf-algorithms/cf-algorithms.h>

#if WITH_TBB
#include <tbb/tbb.h>
#endif

#if WITH_BOOST_THREADS
#include <boost/thread.hpp>
#endif

#if WITH_OPENMP
#include <omp.h>
#endif

namespace m3D {

	using namespace ::std;
	using ::cfa::meanshift::Point;
	using ::cfa::meanshift::FeatureSpace;
	using ::cfa::meanshift::PointIndex;
	using ::cfa::meanshift::SearchParameters;
	using ::cfa::meanshift::Kernel;
	using ::cfa::meanshift::MeanshiftOperation;

    template <typename T> 
    struct ClusterTask
    {
        FeatureSpace<T>         *m_fs;
        PointIndex<T>           *m_index;
        ClusterOperation<T>     *m_op;
        ClusterList<T>          *m_cluster_list;
        
        const Kernel<T>         *m_kernel;
        const WeightFunction<T> *m_weight_function;
        const double            m_drf_threshold;
        
        const SearchParameters  *m_search_params;
        
        const size_t            m_start_index;
        const size_t            m_end_index;
        
        const bool              m_show_progress;
        
    public:
        
        /** Constructor for single-threaded or boost threads 
         */
        ClusterTask( FeatureSpace<T> *fs,  
                     PointIndex<T> *index,
                     ClusterOperation<T> *op,
                     ClusterList<T> *cs,
                     const size_t start_index,
                     const size_t end_index,
                     const SearchParameters *params,
                     const Kernel<T> *kernel,
                     const WeightFunction<T> *weight_function,
                     const double &drf_threshold,
                     const bool show_progress = true ) :
        m_fs(fs), 
        m_index(index),
        m_op(op),
        m_cluster_list(cs),
        m_kernel( kernel),
        m_weight_function(weight_function),
        m_drf_threshold(drf_threshold),
        m_search_params( params ), 
        m_start_index(start_index), 
        m_end_index(end_index), 
        m_show_progress(show_progress)
        {
        };
        
        /** Copy contructor */
        ClusterTask( const ClusterTask &other )
        : m_fs(other.m_fs),
        m_index(other.m_index),
        m_op(other.m_op),
        m_cluster_list(other.m_cluster_list),
        m_kernel(other.m_kernel),
        m_weight_function(other.m_weight_function),
        m_drf_threshold(other.m_drf_threshold),
        m_search_params(other.m_search_params),
        m_start_index(other.m_start_index),
        m_end_index(other.m_end_index),
        m_show_progress(other.m_show_progress)
        {
        }
        
#if WITH_TBB
        
        /** Constructor for TBB
         */
        ClusterTask( FeatureSpace<T> *fs,  
                     PointIndex<T> *index,
                     ClusterOperation<T> *op,
                     ClusterList<T> *cs,
                     const SearchParameters *params,
                     const Kernel<T> *kernel,
                     const WeightFunction<T> *weight_function,
                     const double &drf_threshold,
                     const bool show_progress = true ) :
        m_fs(fs), 
        m_index(index),
        m_op(op),
        m_cluster_list(cs),
        m_start_index(0), 
        m_end_index(0), 
        m_search_params( params ),
        m_kernel( kernel),
        m_weight_function(weight_function),
        m_drf_threshold( drf_threshold ),
        m_show_progress(show_progress)
        {};
        
        void 
        operator()( const tbb::blocked_range<size_t>& r ) const
        {
            MeanshiftOperation<T> ms_op( m_fs, m_index );

            for ( size_t index = r.begin(); index != r.end(); index++ )
            {
                if ( m_show_progress )
                {
                    m_op->increment_cluster_progress();
                }
                
                typename Point<T>::ptr x = m_fs->points[ index ];
                
                bool normalize = !(GRID_ROUNDING_METHOD_NONE == 1);
                
                x->shift = ms_op.meanshift( x->values, m_search_params, m_kernel, m_weight_function, normalize );
                
                vector<T> spatial_shift = m_fs->spatial_component(x->shift);
                
                if (normalize)
                {
                    x->gridded_shift = m_fs->coordinate_system->rounded_gridpoint(spatial_shift);
                }
            }
            
            m_op->report_done();
        }
#endif
        
        void 
        operator()()
        {
            MeanshiftOperation<T> ms_op( m_fs, m_index );
            
               
            #if WITH_OPENMP
            // Prime the index
            ms_op.prime_index(m_search_params);

            #pragma omp parallel for schedule(dynamic) 
            #endif
            for ( size_t index = m_start_index; index < m_end_index; index++ )
            {
                if ( m_show_progress )
                {
                    #if WITH_OPENMP
                    #pragma omp critical
                    #endif
                    m_op->increment_cluster_progress();
                }
                
                typename Point<T>::ptr x = m_fs->points[ index ];
                
                x->shift = ms_op.meanshift( x->values, m_search_params, m_kernel, m_weight_function );
                
                vector<T> spatial_shift = m_fs->spatial_component(x->shift);
                
                x->gridded_shift = m_fs->coordinate_system->rounded_gridpoint(spatial_shift);
            }
            
            m_op->report_done();
        }
    };
}

#endif
