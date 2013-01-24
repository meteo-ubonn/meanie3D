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
    {
        if ( scale < 0 )
        {
            throw logic_error("scale can not be less than zero");
        }
        
        m_scale = scale;
    };

    template <typename T>
    ScaleSpaceFilter<T>::~ScaleSpaceFilter() {};
    
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
    ScaleSpaceFilter<T>::apply( FeatureSpace<T> *fs )
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
    
    
};

#endif
