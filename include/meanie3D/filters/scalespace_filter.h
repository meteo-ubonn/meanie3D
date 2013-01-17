#ifndef _M3D_ScaleSpaceFilter_H_
#define _M3D_ScaleSpaceFilter_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <meanie3D/filters/filter.h>

namespace m3D {

    /** Smoothes the data with a scale-space filter. This filter does NOT create
     * new points. Only the existing points are smoothed out.
     */
    template <class T>
    class ScaleSpaceFilter : public FeatureSpaceFilter<T>
    {
        
    private:
        
        double   m_scale;

    public:
        
#pragma mark -
#pragma mark Constructor/Destructor
        
        /** Constructor
         * @param scale parameter
         * @param show progress indicator while filtering (default no)
         * @throws logic_error if scale < 0
         */
        ScaleSpaceFilter(T scale, bool show_progress = false);
        
        /** Destructor
         */
        virtual ~ScaleSpaceFilter();
        
#pragma mark -
#pragma mark Accessors
        
        /** Sets the scale parameter
         * @param new scale
         * @throws logic_error if scale < 0
         */
        void set_scale( const double scale );
        
        double scale();
        
#pragma mark -
#pragma mark Abstract filter method
        
        virtual void apply( FeatureSpace<T> *fs );
    };
    
};

#endif
