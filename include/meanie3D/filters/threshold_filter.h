#ifndef _M3D_ThresholdFilter_H_
#define _M3D_ThresholdFilter_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

#include <meanie3D/filters/filter.h>

namespace m3D {
    
	using namespace std;

    /** Filters the values in the feature-space based on a set of 
     * thresholds, one for each feature-space dimension. Points that
     * are smaller than the theshold are dropped from the feature-space.
     */
    template <class T>
    class ThresholdFilter : public FeatureSpaceFilter<T>
    {
        
    private:
        
        vector<T>  m_thresholds;
        
    public:
        
#pragma mark -
#pragma mark Constructor/Destructor
        
        /** Constructor
         * @param thresholds
         * @param show progress indicator while filtering (default no)
         * @throws logic_error if |thresholds| = 0
         */
        ThresholdFilter(const vector<T> &thresholds);
        
        virtual ~ThresholdFilter();

#pragma mark - 
#pragma mark Accessors
        
        vector<T> thresholds();

        /** Sets new threshold values
         * @param thresholds
         * @throws logic_error if |thresholds| = 0
         */
        void set_thresholds(const vector<T> &);
        
#pragma mark -
#pragma mark Abstract filter method
        
        virtual void apply( FeatureSpace<T> *fs );
    };
    
}

#endif
