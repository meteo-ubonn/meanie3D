#ifndef _M3D_FeatureSpaceFilter_H_
#define _M3D_FeatureSpaceFilter_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <cf-algorithms/cf-algorithms.h>

namespace m3D {

    /** Abstract base class for feature-space filters.
     */
    template <class T>
    class FeatureSpaceFilter
    {
        
    private:
        
        bool    m_show_progress;
        
    public:
        
#pragma mark -
#pragma mark Constructor/Destructor
        
        /** Default constructor
         */
        FeatureSpaceFilter(bool show_progress) : m_show_progress(show_progress) {};
        
        /** Destructor 
         */
        virtual ~FeatureSpaceFilter() {};

#pragma mark -
#pragma mark Accessors
        
        bool show_progress() { return m_show_progress; };

#pragma mark -
#pragma mark Abstract filter method

        /** The 'apply' method modifies the given feature-space (destructive).
         * If you want to hold on to the data prior to filtering, you need to
         * make a copy.
         * @param feature space
         */
        virtual void apply( FeatureSpace<T> *fs ) = 0;
    };
    
};

#endif
