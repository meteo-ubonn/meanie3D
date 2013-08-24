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
    WeightThresholdFilter<T>::~WeightThresholdFilter() {}
    
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
        
        vector< Point<T> * > accepted;
        vector< Point<T> * > erased;
        
        for ( size_t k=0; k < fs->points.size(); k++)
        {
            if ( this->show_progress() )
            {
                progress_bar->operator++();
            }
            
            Point<T> *p = fs->points.at(k);
            
            T w = m_weight_function->operator()(p);
            
            if ( w < m_lower_threshold || w > m_upper_threshold)
            {
                erased.push_back(p);
            }
            else
            {
                accepted.push_back(p);
            }
        }
        
        fs->points = accepted;
        
        for (size_t i=0; i<erased.size(); i++)
        {
            Point<T> *p = erased[i];
            delete p;
        }

        if ( this->show_progress() )
        {
            delete progress_bar;
            cout << "done. (" << stop_timer() << "s)" << std::endl;
        }
    }

}

#endif
