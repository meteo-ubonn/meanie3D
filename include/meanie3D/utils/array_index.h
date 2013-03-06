#ifndef _M3D_ArrayIndex_H_
#define _M3D_ArrayIndex_H__

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>
#include <netcdf>

#include <cf-algorithms/cf-algorithms.h>

namespace m3D {
    
    // Forward declarations

    using namespace std;
    using namespace netCDF;
    using cfa::meanshift::Point;
    using cfa::meanshift::FeatureSpace;
    using cfa::meanshift::utils::CoordinateSystem;

    template <class T>
    class ArrayIndex
    {
        
#pragma mark -
#pragma mark Attributes

    private:
        
        FeatureSpace<T>     *m_fs;
        void                *data;
        
        /** Helper method for constructing the arbitrary deep
         * vector structure. Guided by the feature-spaces's
         * coordinate system.
         * @param index of current dimension
         */
        void
        construct_array_recursive(size_t dim_index);
        
        void
        destroy_array_recursive(size_t dim_index);
        
#pragma mark -
#pragma mark Constructors/Destructors

    public:

        /** Creates an array index for the given feature-space
         * @param feature space
         */
        ArrayIndex(const FeatureSpace<T> *fs);
        
        /** Destructor 
         */
        virtual ~ArrayIndex();
    
#pragma mark -
#pragma mark Accessors
    
        typename Point<T>::ptr get(const GridPoint &gp);
    
        void set(const GridPoint &gp, typename Point<T>::ptr p);
    }

};
    
#endif
