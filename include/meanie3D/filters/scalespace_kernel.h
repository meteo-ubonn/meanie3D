#ifndef _M3D_ScaleSpaceKernel_H_
#define _M3D_ScaleSpaceKernel_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

namespace m3D {
    
	using namespace std;

    /** Kernel profile g for scale-space filtering. Allows for 
     * pre-rendering at a number of points for speed-up.
     */
    template <class T>
    class ScaleSpaceKernel
    {
        
    private:
        
        T          m_t;
        T          m_gauging_factor;
        vector<T>  m_values;
        
    public:
        
#pragma mark -
#pragma mark Constructor/Destructor
        
        /** Constructor
         * @param scale parameter t
         */
        ScaleSpaceKernel(const T &t);

        /** Constructor. Pre-calculates the values at the given
         * distances. 
         *
         * @param scale parameter t
         * @param distances
         */
        ScaleSpaceKernel(const T &t, const vector<T> &distances);
        
        /** Copy constructor
         * @param other kernel
         */
        ScaleSpaceKernel(const ScaleSpaceKernel<T> &o);
        
        /** Destructor
         */
        virtual ~ScaleSpaceKernel();
        
        /** Copy operator 
         */
        ScaleSpaceKernel<T>
        operator = ( const ScaleSpaceKernel<T>& other );


#pragma mark -
#pragma mark Accessors
        
        /** @return pre-calculated values.
         */
        const vector<T> &values();
        
        /** @return true if the kernel was constructed with pre-sampled values, false else.
         */
        bool isPreSampled();
        
        const T scale() {return m_t;};

#pragma mark -
#pragma mark Calculation
        
        /** Obtains the pre-calculated value at index i
         * @param index
         * @return value
         * @throws
         */
        const T value(size_t index);
        
        /** Calculates the value at distance r
         * @param distance r
         * @return calculated value;
         */
        const T value(T r);
    };    
};

#endif
