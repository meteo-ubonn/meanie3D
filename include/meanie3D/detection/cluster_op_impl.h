#ifndef _M3D_Operation_Cluster_Impl_H_
#define _M3D_Operation_Cluster_Impl_H_

#include "cluster_task.h"

#include <cf-algorithms/cf-algorithms.h>
#include <meanie3D/utils/visit.h>

namespace m3D {

	using namespace cfa::meanshift;

#if WRITE_MODES

    static size_t s_pass_counter = 0;

    template <typename T>
    size_t
    ClusterOperation<T>::pass_counter()
    {
        return s_pass_counter;
    };

    
    template <typename T>
    void
    ClusterOperation<T>::reset_pass_counter()
    {
        s_pass_counter=0;
    };

    template <typename T>
    void
    ClusterOperation<T>::increment_pass_counter()
    {
        s_pass_counter++;
    }
    
#endif

    
    template <typename T>
    void 
    ClusterOperation<T>::increment_cluster_progress() const
    {
#if PROVIDE_THREADSAFETY
        boost::mutex::scoped_lock( m_cluster_progress_mutex );
#endif
        if ( m_progress_bar )
        {
            m_progress_bar->operator++();
        }
    }
    
    template <typename T> 
    void 
    ClusterOperation<T>::report_done() 
    { 
#if PROVIDE_THREADSAFETY
        boost::mutex::scoped_lock( m_cluster_threadcount_mutex );
#endif            
        m_cluster_threadcount++;
    };

    
#pragma mark -
#pragma mark Clustering Code
    
    template <typename T> 
    ClusterList<T> 
    ClusterOperation<T>::cluster( const SearchParameters *params,
                                  const Kernel<T> *kernel,
                                  const int weight_index,
                                  const double &drf_threshold,
                                  const bool show_progress_bar )
    {
        using namespace cfa::utils::timer;
        using namespace cfa::utils::visit;
        
        vector<T> resolution;
        
        if ( params->search_type() == SearchTypeRange )
        {
            RangeSearchParams<T> *p = (RangeSearchParams<T> *) params;
            resolution = p->bandwidth;
        }
        else
        {
            KNNSearchParams<T> *p = (KNNSearchParams<T> *) params;
            resolution = p->resolution;
        }

#if WRITE_TRAJECTORIES
        static size_t trajectory_number = 0;
#endif  

        if ( show_progress_bar )
        {
            cout << endl << "Creating meanshift vector graph ...";
            
            start_timer();
            
            m_progress_bar = new boost::progress_display( this->feature_space->size() );
        }

        ClusterList<T> cluster_list;
        
        // Guard against empty feature-space
        
        if ( this->feature_space->points.size() == 0 )
        {
            cerr << "FeatureSpace is empty" << endl;
            
            return cluster_list;
        }
        
#if WITH_BOOST_THREADS
        
        // TODO: find some dynamic way to this number
        static const size_t NUMBER_OF_CLUSTERING_THREADS = 20;
        
        m_cluster_threadcount = 0;
        size_t num_threads = 0;
        size_t points_per_thread = points.size() / NUMBER_OF_CLUSTERING_THREADS;
        
        // create the threads
        
        for ( size_t i = 0; i < NUMBER_OF_CLUSTERING_THREADS; i++ )
        {
            size_t start_index = i * points_per_thread;
            
            size_t end_index = (( i < NUMBER_OF_CLUSTERING_THREADS - 1 ) ? ( i + 1 ) * points_per_thread - 1 : points.size());
            
            boost::thread thread( ClusterTask<T>(this->feature_space,
                                                 this->feature_space_index,
                                                 this, &cluster_list,
                                                 weight_index,
                                                 this->feature_space->size(),
                                                 params,
                                                 kernel,
                                                 weight_index,
                                                 drf_threshold,
                                                 show_progress_bar) );
            thread.detach();
            
            num_threads++;
        }
        
        // Wait until cluster progress is 100%
        while ( m_cluster_threadcount < num_threads ) {};
        
#elif WITH_TBB
        
        parallel_for(tbb::blocked_range<size_t>(0,points.size()) ,
                     ClusterTask<T>(this->feature_space,
                                    this->feature_space_index,
                                    this, &cluster_list,
                                    weight_index,
                                    this->feature_space->size(),
                                    params,
                                    kernel,
                                    weight_index,
                                    drf_threshold,
                                    show_progress_bar),
                     tbb::auto_partitioner() );
        
#else
        // Single threaded execution
        
        ClusterTask<T> ct(this->feature_space,
                          this->feature_space_index,
                          this, &cluster_list,
                          0,
                          this->feature_space->size(),
                          params,
                          kernel,
                          weight_index,
                          drf_threshold,
                          show_progress_bar);
        ct();
        
#endif
        
#if WRITE_MEANSHIFT_VECTORS
        VisitUtils<T>::write_shift_vectors( "meanshift_vectors.vtk", this->feature_space, false );
        VisitUtils<T>::write_shift_vectors( "meanshift_vectors-spatial.vtk", this->feature_space, true );
#endif
        
        if ( show_progress_bar )
        {
            cout << "done. (" << stop_timer() << "s)" << endl;
            
            delete m_progress_bar;

            m_progress_bar = NULL;
        }
        
        // Analyse the graph and create clusters
        
        cluster_list.aggregate_cluster_graph( (size_t)weight_index, this->feature_space, resolution, show_progress_bar );
        
#if WRITE_BOUNDARIES
        cluster_list.write_boundaries( (size_t) weight_index, this->feature_space, this->feature_space_index, resolution );
#endif

        // Analyze the clusters to create objects
        cluster_list.aggregate_clusters_by_boundary_analysis( (size_t)weight_index, this->feature_space_index, resolution, drf_threshold, show_progress_bar );
        
#if WRITE_MODES
        
        size_t min_size = std::numeric_limits<size_t>::max();
        size_t max_size =std::numeric_limits<size_t>::min();
        for (size_t i=0; i<cluster_list.size(); i++)
        {
            if ( cluster_list[i]->size() < min_size ) { min_size = cluster_list[i]->size(); }
            if ( cluster_list[i]->size() > max_size ) { max_size = cluster_list[i]->size(); }
        }
        cout << "Cluster sizes range: [" << min_size << "," << max_size << "]" << endl;
        
        std::string fn = this->feature_space->filename()
        + "-modes-" + boost::lexical_cast<string>( pass_counter() )
        + ".vtk";
        m3D::utils::VisitUtils<T>::write_cluster_modes_vtk( fn, cluster_list.clusters );
        
        fn = this->feature_space->filename()
        + "-raw-modes-" + boost::lexical_cast<string>( pass_counter() )
        + ".vtk";
        cfa::utils::VisitUtils<T>::write_modes_vtk( fn, cluster_list.trajectory_endpoints(), cluster_list.trajectory_lengths() );
        m3D::utils::VisitUtils<T>::write_cluster_modes_vtk( fn, cluster_list.clusters );
        
        increment_pass_counter();
        
#endif

        return cluster_list;
    }
    
#if 0
    template <typename T> 
    typename Cluster<T>::list 
    ClusterOperation<T>::create_cluster_graph( vector<T>& bandwidth, 
                                          Kernel<T> &kernel, 
                                          WeightFunction<T> *weight_function,
                                          bool show_progress_bar )
    {
        typename Cluster<T>::list result;
        
        if ( show_progress_bar )
        {
            m_progress_bar = new boost::progress_display( 2 * points.size() );
        }
        
        typename Point<T>::listIterator pi;
        
        vector<T> resolution = coordinate_system->resolution();
        
        // calculate all gradient estimates
        
        for ( pi = points.begin(); pi != points.end(); pi++ )
        {
            // Grab the sample
            
            typename Point<T>::list *sample = m_index->range_search( (*pi)->values, bandwidth );
            
            //            cout << "Sample around "; 
            //            print_vector( (*pi)->values );
            //            cout << " contains " << sample->size() << " points." << endl;
            
            (*pi)->shift = sample_meanshift( (*pi)->values, sample, bandwidth, kernel, weight_function );
            
            //            cout << "-> meanshift = ";
            //            print_vector( (*pi)->shift );
            //            cout << endl;
            
            if ( show_progress_bar )
            {
                increment_cluster_progress();
            }
            
            delete sample;
        }
        
        for ( pi = points.begin(); pi != points.end(); pi++ )
        {
            // Grab the sample
            typename Point<T>::list *sample = m_index->range_search( (*pi)->values, bandwidth );
            
            // Find the vector xi in the sample, where the angle between | x - xi | 
            // is closest to the sample meanshift vector
            // (see Fukunaga, Introduction to statistical pattern recognition
            //  second edition, pp539 ff)
            
            Point<T> *steepest_ascending_point = NULL;
            
            T min_angle = std::numeric_limits<T>::max();
            
            //            T max_alignment = std::numeric_limits<T>::min();
            
            typename Point<T>::listIterator si;
            
            vector<T> dx( (*pi)->values.size() );
            
            for ( si = sample->begin(); si != sample->end(); si++ )
            {
                vector_subtract_r( (*si)->values, (*pi)->values , dx );
                
                T angle = vector_angle( (*si)->shift, dx );
                
                if ( angle < min_angle )
                {
                    min_angle = angle;
                    
                    steepest_ascending_point = *si;
                }
                
                //                T alignment = scalar_product( shift, dx );
                //                
                //                if ( alignment > max_alignment )
                //                {
                //                    steepest_ascending_point = *si;
                //                    
                //                    max_alignment = alignment;
                //                }
                
            }
            
            // No self-reference allowed here
            
            if ( steepest_ascending_point != (*pi) )
            {
                (*pi)->next = steepest_ascending_point;
            }
            
            if ( show_progress_bar )
            {
                increment_cluster_progress();
            }
            
            delete sample;
        }
        
        if ( show_progress_bar )
        {
            delete m_progress_bar;
            
            m_progress_bar = NULL;
        }
        
        return result;
    };    

    template <typename T> 
    void 
    ClusterOperation<T>::aggregate_clusters( vector<T>& fuzziness )
    {
        vector< Cluster<T> > aggregates;
        
        while ( m_clusters.size() > 0 )
        {
            Cluster<T> cluster = m_clusters.front();
            
            // find all candidates in the vicinity
            
            vector< Cluster<T> > candidates;
            
            candidates.push_back( cluster );
            
            // some vars
            
            vector<T> mode = cluster.mode();
            
            vector<T> mean_mode( mode.size(), 0 );
            
            if ( m_clusters.size() > 1 )
            {
                for ( size_t j = 1; j < m_clusters.size(); j++ )
                {
                    Cluster<T> candidate = m_clusters[j];
                    
                    vector<T> candidate_mode = candidate.mode();
                    
                    if ( within_range( mode, candidate_mode, fuzziness ) )
                    {
                        candidates.push_back( candidate );
                    }
                }
            }
            
            // calculate arithmetic mean of candidate modes
            
            for ( size_t j = 0; j < candidates.size(); j++ )
            {
                Cluster<T> candidate = candidates[j];
                
                for ( size_t ri = 0; ri < mode.size(); ri++ )
                {
                    mean_mode[ri] += candidate.mode()[ri];
                }
            }
            
            vector_mult_scalar( mean_mode, 1.0 / ((T) candidates.size()) );
            
            // Create aggregate cluster
            
            Cluster<T> new_cluster( mean_mode );
            
            for ( size_t j = 0; j < candidates.size(); j++ )
            {
                new_cluster.add_points( candidates[j].points() );
            }
            
            aggregates.push_back( new_cluster );
            
            // remove the aggregated clusters from the original set
            
            for ( size_t j = 0; j < candidates.size(); j++ )
            {
                Cluster<T> c = candidates[ j ];
                
                typename vector< Cluster<T> >::iterator it = find( m_clusters.begin(), m_clusters.end(), c );
                
                m_clusters.erase( it );
            }
            
            cout << "Aggregating " << candidates.size() << " clusters." << endl;
            cout << "mean mode = ";
            print_vector( mean_mode );
            cout << endl;
            cout << "#clusters remaning = " << m_clusters.size() << endl;
        }
        
        m_clusters = aggregates;
    }
    
#endif

}

#endif
