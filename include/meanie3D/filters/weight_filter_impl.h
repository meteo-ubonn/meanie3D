#ifndef _M3D_WeightThresholdFilter_Impl_H_
#define _M3D_WeightThresholdFilter_Impl_H_

#include <cf-algorithms/utils.h>

namespace m3D {
    
    template <typename T>
    WeightThresholdFilter<T>::WeightThresholdFilter(WeightFunction<T> *w,
                                                    T lower,
                                                    T upper,
                                                    bool show_progress)
    : FeatureSpaceFilter<T>(show_progress)
    , m_weight_function(w)
    , m_lower_threshold(lower)
    , m_upper_threshold(upper)
    {}
    
    template <typename T>
    WeightThresholdFilter<T>::~WeightThresholdFilter() {};
    
#pragma mark -
#pragma mark Abstract filter method
    
    template <typename T>
    void
    WeightThresholdFilter<T>::apply( FeatureSpace<T> *fs )
    {
        using namespace std;
        using namespace cfa::utils::timer;
        
        boost::progress_display *progress_bar = NULL;
        
        if ( this->show_progress() )
        {
            cout << endl << "Applying weight function filter ...";
            
            progress_bar = new boost::progress_display( fs->size() );
            
            start_timer();
        }
        
        typename Point<T>::list::iterator pit;
        
        for ( pit = fs->points.begin(); pit != fs->points.end(); )
        {
            if ( this->show_progress() )
            {
                progress_bar->operator++();
            }
            
            Point<T> *p = *pit;
            
            T w = m_weight_function->operator()(p);
            
            if ( w < m_lower_threshold || w > m_upper_threshold)
            {
                fs->points.erase( pit );
                
                delete p;
            }
            else
            {
                pit++;
            }
        }

        if ( this->show_progress() )
        {
            delete progress_bar;
            cout << "done. (" << stop_timer() << "s)" << std::endl;
        }

        
    }

};

#endif
