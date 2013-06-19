#ifndef _M3D_ScalarIndex_H_
#define _M3D_ScalarIndex_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>
#include <netcdf>

#include <cf-algorithms/cf-algorithms.h>

namespace m3D {
    
    // Forward declarations

    using namespace std;
    using namespace netCDF;
    using cfa::utils::coords::CoordinateSystem;

    template <class T>
    class ScalarIndex
    {
        
    protected:
        
        typedef vector<void*>   array_t;
        typedef array_t*        array_t_ptr;
        
#pragma mark -
#pragma mark Attributes

    private:

        const CoordinateSystem<T>   *m_coordinate_system;
        array_t                     *m_data;
        
        /** Called by build_index, copy constructor
         */
        void
        construct_recursive(size_t dim_index,
                            array_t **array,
                            typename CoordinateSystem<T>::GridPoint &gridpoint,
                            ScalarIndex<T> *other = NULL);
        
        /** Called by destructor.
         */
        void
        destroy_recursive(size_t dim_index,
                          array_t **array,
                          typename CoordinateSystem<T>::GridPoint &gridpoint);
        
        void
        copy_recursive(ScalarIndex<T> *otherIndex,
                       size_t dim_index,
                       typename CoordinateSystem<T>::GridPoint &gridpoint);
        
#pragma mark -
#pragma mark Constructors/Destructors

    public:

        /** Constructs an array index for the given coordinate system.
         * @param coordinate system
         */
        ScalarIndex(CoordinateSystem<T> *cs);
        
        /** Copy constructor on pointer
         * @param pointer to array index
         */
        ScalarIndex(ScalarIndex<T> *other);
        
        /** Destructor
         */
        virtual ~ScalarIndex();
        
#pragma mark -
#pragma mark Accessors

        /** Gets value at given gridpoint
         * @param grid point
         */
        T get(const typename CoordinateSystem<T>::GridPoint &gp) const;
        
        /** Gets value at given gridpoint
         * @param coordinate
         */
        T get(const typename CoordinateSystem<T>::Coordinate &c) const;

    
        /** Sets a value in the index at a given grid point.
         * @param grid point
         * @param value
         */
        void
        set(const typename CoordinateSystem<T>::GridPoint &gp, T value);
        
        /** Sets a value in the index at a given grid point.
         * @param coordinate
         * @param value
         */
        void
        set(const typename CoordinateSystem<T>::Coordinate &gp, T value);

    };

};
    
#endif
