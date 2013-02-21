#ifndef _M3D_ScaleSpaceFilter_H_
#define _M3D_ScaleSpaceFilter_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <meanie3D/filters/filter.h>
#include <cf-algorithms/cf-algorithms.h>
#include <boost/progress.hpp>

namespace m3D {

    /** Smoothes the data with a scale-space filter. This filter does NOT create
     * new points. Only the existing points are smoothed out.
     */
    template <class T>
    class ScaleSpaceFilter : public FeatureSpaceFilter<T>
    {
        
    private:

        double          m_scale;            /// Sole input parameter
        
        PointIndex<T>   *m_index;           /// Index of a copy of the FS
        
        map<int,T>      m_min;              /// minimum tracker
        
        map<int,T>      m_max;              /// maximum tracker
        
        double          m_filter_width;     /// filter size, stored for efficiency
        
        double          m_gauging_factor;   /// multiplier, stored for efficiency
        
        cfa::meanshift::RangeSearchParams<T>    *m_range;           /// Search range, stored for efficiency
        
        boost::progress_display                 *m_progress_bar;    /// Progress meter

        /** Shortcut method. Does not create new points. Violates
         * the integrational constraint of the scale-space theory,
         * where the integral over the whole domain should not change
         */
        void applyWithoutNewPoints( FeatureSpace<T> *fs );
        
        /** Proper method. Preserves the integral.
         */
        void applyWithNewPoints( FeatureSpace<T> *fs );

        /** Recursive helper method for the proper method 
         */
        void applyWithNewPointsRecursive( FeatureSpace<T> *fs,
                                          size_t dimensionIndex,
                                          typename cfa::utils::CoordinateSystem<T>::GridPoint &gridpoint );

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
