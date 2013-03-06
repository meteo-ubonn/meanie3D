#ifndef _M3D_ScaleSpaceFilter_Impl_H_
#define _M3D_ScaleSpaceFilter_Impl_H_

#include <exception>
#include <stdexcept>
#include <cmath>
#include <map>

#include <boost/progress.hpp>

#include <cf-algorithms/cf-algorithms.h>

namespace m3D {

	using namespace std;
	using namespace ::cfa::meanshift;
	using ::cfa::utils::VisitUtils;
	using ::cfa::utils::coords::CoordinateSystem;

	using namespace std;
 
    template <typename T>
    ScaleSpaceFilter<T>::ScaleSpaceFilter(T scale, bool show_progress) : FeatureSpaceFilter<T>(show_progress)
    , m_index(NULL)
    , m_progress_bar(NULL)
    , m_range(NULL)
    {
        if ( scale < 0 )
        {
            throw logic_error("scale can not be less than zero");
        }
        
        m_scale = scale;
    };

    template <typename T>
    ScaleSpaceFilter<T>::~ScaleSpaceFilter()
    {
        if ( m_index != NULL )
        {
            delete m_index;
            m_index = NULL;
        }
        
        if ( m_progress_bar != NULL )
        {
            delete m_progress_bar;
            m_progress_bar = NULL;
        }
        
        if ( m_range != NULL )
        {
            delete m_range;
            m_range = NULL;
        }
    };
    
#pragma mark -
#pragma mark Accessors
    
    template <typename T>
    void
    ScaleSpaceFilter<T>::set_scale( const double scale )
    {
        if ( scale < 0 )
        {
            throw std::logic_error("scale can not be less than zero");
        }

        m_scale = scale;
    };
    
    template <typename T>
    double
    ScaleSpaceFilter<T>::scale() { return m_scale; };
    
#pragma mark -
#pragma mark Abstract filter method
    
    template <typename T>
    void
    ScaleSpaceFilter<T>::applyWithoutNewPoints( FeatureSpace<T> *fs )
    {
      using namespace cfa::utils::timer;
      using namespace cfa::utils::vectors;
 
#if WRITE_FEATURESPACE
        std::string fn = fs->filename() + "_scale_0.vtk";
        VisitUtils<T>::write_featurespace_vtk( fn, fs );
#endif
        if ( m_scale > 0.0 )
        {
            // Replace the value of each point in feature-space with the average
            // of the neighbours within a range determined by the scale parameter
            
            double filter_width = sqrt( ceil( -2.0 * m_scale * log(0.01) ) ) / 2;
            
            vector<T> spatial_range( fs->coordinate_system->size(), filter_width );
            
            boost::progress_display *progress_bar = NULL;
            
            if ( this->show_progress() )
            {
                cout << endl << "Applying scale filter t=" << m_scale << " (ranges " << spatial_range << ") ...";

                progress_bar = new boost::progress_display( fs->size() );

                start_timer();
            }
            
            RangeSearchParams<T> params( spatial_range );
            
            // Make a copy of the feature space
            
            FeatureSpace<T> fs_copy = FeatureSpace<T>( *fs );
            
            // create a 'spatial only' index of the feature-space points
            
            PointIndex<T> *index = PointIndex<T>::create( &fs_copy, fs->coordinate_system->dimension_variables() );

            // Initialize limit trackers
            
            map<int,T> min;
            map<int,T> max;

            // Pre-calculate the gauging factor for the gaussian kernel
            
            double gauging_factor = 1.0 / std::pow( (double)sqrt(2.0 * M_PI * m_scale), (double)fs->coordinate_system->size() );

            // Iterate over feature-space

            typename Point<T>::list::iterator pit;
            
            for ( pit = fs->points.begin(); pit != fs->points.end(); pit++ )
            {
                if ( this->show_progress() )
                {
                    progress_bar->operator++();
                }
                
                Point<T> *p = *pit;
                
                // search in range filter_width around p in spatial index
                
                typename Point<T>::list *sample = index->search( p->coordinate, &params );
                
                // Iterate over the range variables in the sample
                
                for ( size_t var_index = fs->coordinate_system->size(); var_index < fs->feature_variables().size(); var_index++ )
                {
                    // replace the point's value at var_index with the gauss-weighted (by distance)
                    // weighed sum of it's neighbours
                    
                    T sum = 0.0;
                    
                    typename Point<T>::list::iterator sample_iter;
                    
                    for ( sample_iter = sample->begin(); sample_iter != sample->end(); sample_iter++ )
                    {
                        typename Point<T>::ptr sample_point = *sample_iter;
                        
                        // distance from origin?
                        
                        T r = vector_norm( p->coordinate - sample_point->coordinate );
                        
                        double gauss = gauging_factor * exp(-(r*r)/(2*m_scale));
                        
                        // weighted sum
                        
                        sum += gauss * sample_point->values[var_index];
                    }
                    
                    // replace the value by the kernel-weighted sum
                    
                    p->values[var_index] = sum;
                    
                    // make sure limit trackers are initialized
                    
                    typename map<int,T>::iterator m;
                    
                    if ( (m = min.find(var_index)) == min.end() )
                    {
                        min[var_index] = std::numeric_limits<T>::max();
                    }
                    
                    if ( (m = max.find(var_index)) == max.end() )
                    {
                        max[var_index] = std::numeric_limits<T>::min();
                    }
                    
                    // track limits
                    
                    if (p->values[var_index] < min[var_index] )
                    {
                        min[var_index] = p->values[var_index];
                    }
                    
                    if (p->values[var_index] > max[var_index] )
                    {
                        max[var_index] = p->values[var_index];
                    }
                }
                
                delete sample;
            }
            
            delete index;
            
            if ( this->show_progress() )
            {
                delete progress_bar;
                cout << "done. (" << stop_timer() << "s)" << endl;
                cout << "Value ranges of smoothed features:" << endl;
                for ( size_t var_index = fs->coordinate_system->size(); var_index < fs->feature_variables().size(); var_index++ )
                {
                    cout << "\t#" << var_index << " : [" << min[var_index] << "," << max[var_index] << "]" << endl;
                }
            }
        }
        
#if WRITE_FEATURESPACE
        fn = fs->filename() + "_scale_" + boost::lexical_cast<string>(m_scale) + ".vtk";
        VisitUtils<T>::write_featurespace_vtk( fn, fs );
#endif
    }
    
    
    template <typename T>
    void
    ScaleSpaceFilter<T>::applyWithNewPointsRecursive( FeatureSpace<T> *fs,
                                                      size_t dimensionIndex,
                                                      typename cfa::utils::CoordinateSystem<T>::GridPoint &gridpoint )
    {
        using namespace std;
        
        NcDim dim = fs->coordinate_system->dimensions()[dimensionIndex];
        
        // iterate over dimensions
        
        for ( int index = 0; index < dim.getSize(); index++ )
        {
            // at each point for each variable

            gridpoint[dimensionIndex] = index;
            
            if ( dimensionIndex < (gridpoint.size()-1) )
            {
                // recurse
                
                applyWithNewPointsRecursive(fs, dimensionIndex+1, gridpoint);
            }
            else
            {
                // Get the points in radius
                if ( this->show_progress() )
                {
                    m_progress_bar->operator++();
                }
                
                vector<T> x(gridpoint.size(),0);
                
                fs->coordinate_system->lookup(gridpoint,x);
                
                // search in range filter_width around p in spatial index
                // Iterate over the range variables in the sample
                // calculate the result from the points found in the index

                typename Point<T>::list *sample = m_index->search( x, m_range );
                
                if ( sample->size() == 0 )
                {
                    delete sample;
                    
                    continue;
                }

                vector<T> values( fs->feature_variables().size(), 0.0 );
                
                typename Point<T>::list::iterator sample_iter;
                
                for ( sample_iter = sample->begin(); sample_iter != sample->end(); sample_iter++ )
                {
                    typename Point<T>::ptr sample_point = *sample_iter;
                    
                    // distance from origin?
                    
                    // T r = vector_norm( x - sample_point->coordinate );
                    
                    vector<T> dx = x - sample_point->coordinate;
                    
                    double gauss = m_gauging_factor * exp( - (dx * dx) / 2 * m_scale );
                    
                    for ( size_t var_index = fs->coordinate_system->size(); var_index < fs->feature_variables().size(); var_index++ )
                    {
                        values[var_index] += gauss * sample_point->values[var_index];
                    }
                }
                
                // make sure limit trackers are initialized
                
                typename map<int,T>::iterator m;
                
                for ( size_t var_index = fs->coordinate_system->size(); var_index < fs->feature_variables().size(); var_index++ )
                {
                    if ( (m = m_min.find(var_index)) == m_min.end() )
                    {
                        m_min[var_index] = std::numeric_limits<T>::max();
                    }
                    
                    if ( (m = m_max.find(var_index)) == m_max.end() )
                    {
                        m_max[var_index] = std::numeric_limits<T>::min();
                    }
                    
                    // track limits
                    
                    if (values[var_index] < m_min[var_index] )
                    {
                        m_min[var_index] = values[var_index];
                    }
                    
                    if (values[var_index] > m_max[var_index] )
                    {
                        m_max[var_index] = values[var_index];
                    }
                }
                
                delete sample;

                // Create point
                
                for ( size_t i = 0; i < fs->coordinate_system->size(); i++ )
                {
                    values[i] = x[i];
                }
            
                typename Point<T>::ptr p = PointFactory<T>::get_instance()->create(x,values);
                
                M3DPoint<T> *mp = (M3DPoint<T> *)p;
                
                mp->setIsOriginalPoint(false);
                
                fs->points.push_back(p);
            }
        }
    }
    
    
    template <typename T>
    void
    ScaleSpaceFilter<T>::applyWithNewPoints( FeatureSpace<T> *fs )
    {
        using namespace cfa::utils::timer;
        using namespace cfa::utils::vectors;
        
        using namespace cfa::utils::timer;
        using namespace cfa::utils::vectors;
        
#if WRITE_FEATURESPACE
        std::string fn = fs->filename() + "_scale_0.vtk";
        VisitUtils<T>::write_featurespace_vtk( fn, fs );
#endif
        if ( m_scale > 0.0 )
        {
            // Replace the value of each point in feature-space with the average
            // of the neighbours within a range determined by the scale parameter
            
            m_filter_width = sqrt( ceil( -2.0 * m_scale * log(0.01) ) ) / 2;
            
            // Pre-calculate the gauging factor for the gaussian kernel
            
            m_gauging_factor = 1.0 / std::pow( 2.0 * M_PI * m_scale, (double)fs->coordinate_system->size()/2 );

            vector<T> spatial_range( fs->coordinate_system->size(), m_filter_width );
            
            m_range = new RangeSearchParams<T>( spatial_range );
            
            if ( this->show_progress() )
            {
                cout << endl << "Applying scale filter t=" << m_scale << " (ranges " << spatial_range << ") ...";
                
                long numPoints = 1;
                
                for ( size_t i=0; i < fs->coordinate_system->dimensions().size(); i++)
                {
                    numPoints *= fs->coordinate_system->dimensions()[i].getSize();
                }
                
                m_progress_bar = new boost::progress_display( numPoints );
                
                start_timer();
            }
            
            // create a 'spatial only' index of the feature-space points
            
            m_index = PointIndex<T>::create( fs, fs->coordinate_system->dimension_variables() );
            
            // Iterate over the dimensions
            
            typename cfa::utils::CoordinateSystem<T>::GridPoint gp = fs->coordinate_system->newGridPoint();
            
            FeatureSpace<T> *filtered_featurespace = new FeatureSpace<T>(*fs,false);
            
            this->applyWithNewPointsRecursive(filtered_featurespace, 0, gp);
            
            // replace points in the featurespace with
            // the filtered results
            
            fs->clear_points();
            
            fs->points = filtered_featurespace->points;
            
            // Finish 

            if ( this->show_progress() )
            {
                cout << "done. (" << stop_timer() << "s)" << endl;
                cout << "Value ranges of smoothed features:" << endl;
                for ( size_t var_index = fs->coordinate_system->size(); var_index < fs->feature_variables().size(); var_index++ )
                {
                    cout << "\t#" << var_index << " : [" << m_min[var_index] << "," << m_max[var_index] << "]" << endl;
                }

                delete m_progress_bar;
                m_progress_bar = NULL;
            }
#if WRITE_FEATURESPACE
            fn = fs->filename() + "_scale_" + boost::lexical_cast<string>(m_scale) + ".vtk";
            VisitUtils<T>::write_featurespace_vtk( fn, fs );
#endif
        }
    }
    
    template <typename T>
    void
    ScaleSpaceFilter<T>::apply( FeatureSpace<T> *fs )
    {
        this->applyWithNewPoints( fs );
        //this->applyWithoutNewPoints( fs );
    }
    
};

#endif
