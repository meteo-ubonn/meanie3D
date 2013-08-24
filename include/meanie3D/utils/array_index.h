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
        
        const CoordinateSystem<T>   *m_coordinate_system;
        array_t                     *m_data;
        bool                        m_make_copies;
        
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
        
        void
        copy_points_recursive(ArrayIndex<T> *otherIndex,
                              size_t dim_index,
                              typename CoordinateSystem<T>::GridPoint &gridpoint);
        
        void
        count_recursive(size_t dim_index,
                        array_t *array,
                        typename CoordinateSystem<T>::GridPoint &gridpoint,
                        size_t &count,
                        bool originalPointsOnly=false);
        
        void
        clear_recursive(size_t dim_index,
                        array_t *array,
                        typename CoordinateSystem<T>::GridPoint &gridpoint,
                        bool delete_points);
        
#pragma mark -
#pragma mark Constructors/Destructors

    public:

        /** Constructs an array index for the given coordinate system.
         * @param coordinate system
         * @param if <code>true</code>, the index makes copies of the
         *        original points. If <code>false</code> it simply holds
         *        references
         */
        ArrayIndex(CoordinateSystem<T> *cs,
                   bool make_copies);
        
        /** Constructs an array index vor the given coordinate system
         * and indexes the point list.
         * @param coordinate system
         * @param point list
         * @param if <code>true</code>, the index makes copies of the
         *        original points. If <code>false</code> it simply holds
         *        references
         */
        ArrayIndex(CoordinateSystem<T> *cs,
                   const typename Point<T>::list &points,
                   bool make_copies);
        
        /** Copy constructor on pointer
         * @param pointer to array index
         * @param if <code>true</code>, the index makes copies of the
         *        original points. If <code>false</code> it simply holds
         *        references
         */
        ArrayIndex(ArrayIndex<T> *other);
        
        /** Destructor
         */
        virtual ~ArrayIndex();
        
#pragma mark -
#pragma mark Indexing operation

        /** Indexes the point list. All points are added to the
         * index (and are copied in the process). Existing points
         * are overwritten (and released).
         * @param point list
         */
        void
        index(const typename Point<T>::list &list);

#pragma mark -
#pragma mark Accessors
    
        typename Point<T>::ptr
        get(const typename CoordinateSystem<T>::GridPoint &gp);
    
        /** Sets a point in the index at a given grid point. If a point exists
         * at the given coordinate, it is released and replaced.
         *
         * @param grid point vector
         * @param pointer to the point object
         * @param if <code>true</code>, a copy is made, if <code>false</code> 
         *        the pointer is used as is. Defaults to <code>true</code>
         */
        void
        set(const typename CoordinateSystem<T>::GridPoint &gp,
            typename Point<T>::ptr p, bool copy=true);
        
        const CoordinateSystem<T> *
        coordinate_system() {return m_coordinate_system;}
        
        /** Clear the index.
         * @param If <code>true</code> indexed points are deleted. 
         *        If <code>false</code> they are left alone.
         */
        void
        clear(bool delete_points);
        
#pragma mark -
#pragma mark Misc
        
        /** Remove the points from the given feature-space and replace them
         * by the points in this index. The points in the original list are
         * released. 
         */
        void
        replace_points(typename Point<T>::list &points);
        
        /** Iterates over the structure and counts the number of non-null
         * entries in the index.
         */
        size_t count(bool originalPointsOnly=false);

    };

}
    
#endif
