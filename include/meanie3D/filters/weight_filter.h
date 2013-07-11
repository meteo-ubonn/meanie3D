#ifndef _M3D_WeightThresholdFilter_H_
#define _M3D_WeightThresholdFilter_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

#include <cf-algorithms/featurespace/weight_function.h>
#include <meanie3D/filters/filter.h>

namespace m3D {
    
	using namespace std;
    using cfa::meanshift::WeightFunction;

    /** Filters the values in the feature-space based on a set of 
     * thresholds, one for each feature-space dimension. Points that
     * are smaller than the theshold are dropped from the feature-space.
     */
    template <class T>
    class WeightThresholdFilter : public FeatureSpaceFilter<T>
    {
    private:
        
        WeightFunction<T>   *m_weight_function;
        T                   m_lower_threshold;
        T                   m_upper_threshold;
        
    public:
        
#pragma mark -
#pragma mark Constructor/Destructor
        
        /** Constructor
         * @param weight function
         * @param lower boundary
         * @param upper boundary
         * @param show progress indicator while filtering (default no)
         * @throws logic_error if |thresholds| = 0
         */
        WeightThresholdFilter(WeightFunction<T> *weight_function,
                              T lower_threshold=numeric_limits<T>::min(),
                              T upper_threshold=numeric_limits<T>::max(),
                              bool show_progress=false);
        
        virtual ~WeightThresholdFilter();

#pragma mark -
#pragma mark Abstract filter method
        
        virtual void apply( FeatureSpace<T> *fs );
    };
    
};

#endif
