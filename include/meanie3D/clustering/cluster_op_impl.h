#ifndef M3D_OPERATION_CLUSTER_IMPL_H
#define M3D_OPERATION_CLUSTER_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/parallel.h>
#include <meanie3D/utils.h>

#include <vector>

#include "cluster_task.h"

namespace m3D {

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
    }


#pragma mark -
#pragma mark Clustering Code

    template <typename T> 
    ClusterList<T> 
    ClusterOperation<T>::cluster(const SearchParameters *params,
                                 const Kernel<T> *kernel,
                                 const WeightFunction<T> *weight_function,
                                 const bool coalesceWithStrongestNeighbour,
                                 const bool write_meanshift_vectors,
                                 const bool show_progress_bar)
    {
        using namespace m3D::utils::vectors;
        
        vector<T> resolution;

        if ( params->search_type() == SearchTypeRange )
        {
            RangeSearchParams<T> *p = (RangeSearchParams<T> *) params;

            // Physical grid resolution in the
            // spatial range
            resolution = this->feature_space->coordinate_system->resolution();

            resolution = ((T)4.0) * resolution;

            // Supplement with bandwidth values for
            // the value range
            for (size_t i=resolution.size(); i < p->bandwidth.size(); i++)
            {
                resolution.push_back(p->bandwidth[i]);
            }
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

        // Compile list of dimension names

        const CoordinateSystem<T> *cs = this->feature_space->coordinate_system;

        NcFile file(m_data_store->filename(), NcFile::read);

        vector<NcVar> vars(cs->dimension_variables());

        for (size_t i=0; i < m_data_store->rank(); i++)
        {
            NcVar variable = file.getVar(m_data_store->variable_names()[i]);
            vars.push_back(variable);
        }

        ClusterList<T> cluster_list(vars,cs->dimensions(),m_data_store->filename());

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
                                                 this->point_index,
                                                 this, &cluster_list,
                                                 weight_index,
                                                 this->feature_space->size(),
                                                 params,
                                                 kernel,
                                                 weight_function,
                                                 drf_threshold,
                                                 show_progress_bar) );
            thread.detach();

            num_threads++;
        }

        // Wait until cluster progress is 100%
        while ( m_cluster_threadcount < num_threads ) {};

#elif WITH_TBB

        // perform one search operation on the index
        // in order to pre-empt the index building

        typename Point<T>::ptr x = this->feature_space->points[0];
        typename Point<T>::list *list = this->point_index->search(x->values, params);
        delete list;

        // Run parallelized version

        // tbb::task_scheduler_init init(1);

        parallel_for(tbb::blocked_range<size_t>(0,this->point_index->size()),

                     ClusterTask<T>(this->feature_space,
                                    this->point_index,
                                    this, &cluster_list,
                                    0,
                                    this->feature_space->size(),
                                    params,
                                    kernel,
                                    weight_function,
                                    show_progress_bar),
                     tbb::auto_partitioner() );

#else
        // Single threaded execution

        ClusterTask<T> ct(this->feature_space,
                          this->point_index,
                          this, &cluster_list,
                          0,
                          this->feature_space->size(),
                          params,
                          kernel,
                          weight_function,
                          show_progress_bar);
        ct();

#endif

        if (write_meanshift_vectors)
        {
            VisitUtils<T>::write_shift_vectors( "meanshift_vectors.vtk", this->feature_space, false );
            VisitUtils<T>::write_shift_vectors( "meanshift_vectors-spatial.vtk", this->feature_space, true );
        }

        if ( show_progress_bar )
        {
            cout << "done. (" << stop_timer() << "s)" << endl;

            delete m_progress_bar;

            m_progress_bar = NULL;
        }

        // Analyse the graph and create clusters

        cluster_list.aggregate_cluster_graph(this->feature_space,weight_function,coalesceWithStrongestNeighbour,show_progress_bar);

        #if WRITE_BOUNDARIES
            cluster_list.write_boundaries( weight_function, this->feature_space, this->point_index, resolution );
        #endif

        #if WRITE_MODES
            size_t min_size = std::numeric_limits<size_t>::max();
            size_t max_size =std::numeric_limits<size_t>::min();
            for (size_t i=0; i<cluster_list.size(); i++)
            {
                if ( cluster_list[i]->size() < min_size ) { min_size = cluster_list[i]->size(); }
                if ( cluster_list[i]->size() > max_size ) { max_size = cluster_list[i]->size(); }
            }
            cout << "Cluster sizes range: [" << min_size << "," << max_size << "]" << endl;

            NetCDFDataStore<T> *ds = (NetCDFDataStore<T> *) this->feature_space->data_store();
            std::string fn = ds->filename()+"-modes-"+boost::lexical_cast<string>(pass_counter())+".vtk";
            VisitUtils<T>::write_cluster_modes_vtk( fn, cluster_list.clusters );

            fn = ds->filename()+"-raw-modes-"+boost::lexical_cast<string>(pass_counter())+".vtk";
            VisitUtils<T>::write_modes_vtk( fn, cluster_list.trajectory_endpoints(), cluster_list.trajectory_lengths() );
            VisitUtils<T>::write_cluster_modes_vtk( fn, cluster_list.clusters );
        #endif

        return cluster_list;
    }
}

#endif