#ifndef _M3D_ArrayIndex_H_
#define _M3D_ArrayIndex_H_

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

        typename Point<T>::list     m_points;
        const CoordinateSystem<T>   *m_coordinate_system;
        array_t                     *m_data;
        
        /** Called by build_index, copy constructor
         */
        void
        construct_array_recursive(size_t dim_index,
                                  array_t **array,
                                  typename CoordinateSystem<T>::GridPoint &gridpoint,
                                  ArrayIndex<T> *other = NULL);
        
        /** Called by destructor.
         */
        void
        destroy_array_recursive(size_t dim_index,
                                array_t **array,
                                typename CoordinateSystem<T>::GridPoint &gridpoint);
        
        /** Called by replace_points.
         */
        void
        replace_points_recursive(typename Point<T>::list &points,
                                 size_t dim_index,
                                 typename CoordinateSystem<T>::GridPoint &gridpoint);
        
         /** Constructs the actual index.
         */
        void
        build_index();
        
        void
        count_recursive(size_t dim_index,
                        array_t *array,
                        typename CoordinateSystem<T>::GridPoint &gridpoint,
                        size_t &count);
        
        
#pragma mark -
#pragma mark Constructors/Destructors

    public:

        /** Creates an array index for the given feature-space
         * @param feature space
         */
        ArrayIndex(const FeatureSpace<T> *fs);
        
        /** Constructs an array index using coordinate system
         * and point list instead.
         * @param coordinate system
         * @param point list
         */
        ArrayIndex(CoordinateSystem<T> *cs, typename Point<T>::list &points);
        
        /** Copy constructor on pointer
         * @param pointer to array index
         */
        ArrayIndex(ArrayIndex<T> *other);
        
        /** Destructor
         */
        virtual ~ArrayIndex();
    
#pragma mark -
#pragma mark Accessors
    
        typename Point<T>::ptr
        get(const typename CoordinateSystem<T>::GridPoint &gp);
    
        void
        set(const typename CoordinateSystem<T>::GridPoint &gp,
            typename Point<T>::ptr p);
        
#pragma mark -
#pragma mark Misc
        
        /** Remove the points from the given feature-space and replace them
         * by the points in this index
         */
        void
        replace_points(typename Point<T>::list &points);
        
        /** Iterates over the structure and counts the number of non-null
         * entries in the index.
         */
        size_t count();

    };

};
    
#endif
