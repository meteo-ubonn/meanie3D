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
    using cfa::utils::coords::CoordinateSystem;

    template <class T>
    class ArrayIndex
    {
        
    protected:
        
        typedef vector<void*>   array_t;
        typedef array_t*        array_t_ptr;
        
#pragma mark -
#pragma mark Attributes

    private:
        
        const FeatureSpace<T>       *m_fs;
        array_t                     *m_data;
        
        /** Helper method for constructing the arbitrary deep
         * vector structure. Guided by the feature-spaces's
         * coordinate system.
         * @param index of current dimension
         * @param pointer to data
         * @param current grid point
         */
        void
        construct_array_recursive(size_t dim_index, array_t **array, typename CoordinateSystem<T>::GridPoint &gridpoint );
        
        /** Helper method for destroying the data. Does not free
         * the points, only the index arrays.
         * @param index of current dimension
         * @param pointer to data
         * @param current grid point
         */
        void
        destroy_array_recursive(size_t dim_index, array_t **array, typename CoordinateSystem<T>::GridPoint &gridpoint );
        
        /** Constructs the actual index.
         */
        void
        build_index();
        
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
    
        typename Point<T>::ptr get(const typename CoordinateSystem<T>::GridPoint &gp);
    
        void set(const typename CoordinateSystem<T>::GridPoint &gp, typename Point<T>::ptr p);
    };

};
    
#endif
