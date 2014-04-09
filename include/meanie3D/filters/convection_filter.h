#ifndef _M3D_ReflectivityConvectionFilter_H_
#define _M3D_ReflectivityConvectionFilter_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

#include <meanie3D/filters/filter.h>

namespace m3D {
    
	using namespace std;

    /** Implementes the stratiform/convective scheme based on the paper:
     * "Climatological Characterization of Three-Dimensional Storm Structure 
     * from Operational Radar and Rain Gauge Data" Steiner, Houze & Yuter, 1995
     * DOI 10.1175/1520-0450(1995)034<1978:CCOTDS>2.0.CO;2
     *
     */
    template <class T>
    class ConvectionFilter : public FeatureSpaceFilter<T>
    {
    private:
        
        vector<T>   m_bandwidth;
        size_t      m_index_of_z;
        T           m_convective_threshold;
        T           m_critical_delta_z;
        T           m_convective_radius_factor;
        bool        m_erase_non_convective;
        
    public:
        
#pragma mark -
#pragma mark Constructor/Destructor

        /** Constructor
         * @param bandwidth     search radius for the background reflectivity
         * @param index_of_z    in the list of feature-space variables, which is Z ?
         * @param z_convective  reflectivity [dbZ] above which everything is
         *                      considered convective
         * @param critical_delta_z  the deltaZ that indicates convection in the 
         *                      cases where Z is below z_convective
         * @param convective_radius_factor  fraction of the bandwidth that determines
         *                      the radius where points are marked as convective around
         *                      a peak point found
         * @param erase_non_convective  if true, points are deleted instead being set to zero
         * @param show_progress if true, progress bar and status messages are displayed
         */
        ConvectionFilter(const vector<T> &bandwidth,
                         const size_t index_of_z,
                         const bool show_progress=false,
                         const T z_convective = 40.0,
                         const T critical_delta_z = 4.5,
                         const T convective_radius_factor = 0.2,
                         const bool erase_non_convective = false);
        
        virtual ~ConvectionFilter();

#pragma mark -
#pragma mark Abstract filter method
        
        virtual void apply( FeatureSpace<T> *fs );
    };
    
}

#endif
