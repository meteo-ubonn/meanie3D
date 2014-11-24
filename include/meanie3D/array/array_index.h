#ifndef M3D_ARRAYINDEX_H
#define M3D_ARRAYINDEX_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils.h>
#include <meanie3D/featurespace/point.h>

#include <vector>
#include <netcdf>

namespace m3D { 
    
    using std::vector;

    template <class T>
    class ArrayIndex
    {

    protected:

        typedef vector<void*>   array_t;
        typedef array_t*        array_t_ptr;

#pragma mark -
#pragma mark Attributes

    private:

        vector<size_t>              m_dimensions;
        array_t                     *m_data;
        bool                        m_make_copies;

        /** Called by build_index, copy constructor
         */
        void
        construct_array_recursive(size_t dim_index,
                                  array_t **array,
                                  vector<int> &gridpoint,
                                  ArrayIndex<T> *other = NULL);

        /** Called by destructor.
         */
        void
        destroy_array_recursive(size_t dim_index,
                                array_t **array,
                                vector<int> &gridpoint);

        /** Copy all points from the final dimension into the list provided
         */
        void
        add_points_to_list(typename Point<T>::list &points,
                           vector<int> &gridpoint);

        /** Called by replace_points.
         */
        void
        replace_points_recursive(typename Point<T>::list &points,
                                 size_t dim_index,
                                 vector<int> &gridpoint);

        void
        get_points_from_index(ArrayIndex<T> *otherIndex,
                              vector<int> &gridpoint);
        void
        copy_points_recursive(ArrayIndex<T> *otherIndex,
                              size_t dim_index,
                              vector<int> &gridpoint);

        void
        count_recursive(size_t dim_index,
                        array_t *array,
                        vector<int> &gridpoint,
                        size_t &count,
                        bool originalPointsOnly=false);

        void
        clear_recursive(size_t dim_index,
                        array_t *array,
                        vector<int> &gridpoint,
                        bool delete_points);
        
        void
        find_neighbours_recursive(vector<int> &gridpoint,
                                  size_t dim_index,
                                  typename Point<T>::list &list,
                                  size_t reach);

#pragma mark -
#pragma mark Constructors/Destructors

    public:

        /** Constructs an array index for the given dimensions
         * @param dimensions
         * @param if <code>true</code>, the index makes copies of the
         *        original points. If <code>false</code> it simply holds
         *        references
         */
        ArrayIndex(const vector<size_t> &dimensions,
                   bool make_copies);

        /** Constructs an array index vor the given dimensions
         * and indexes the point list.
         * @param coordinate system
         * @param point list
         * @param if <code>true</code>, the index makes copies of the
         *        original points. If <code>false</code> it simply holds
         *        references
         */
        ArrayIndex(const vector<size_t> &dimensions,
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
        get(const vector<int> &gp);

        /** Sets a point in the index at a given grid point. If a point exists
         * at the given coordinate, it is released and replaced.
         *
         * @param grid point vector
         * @param pointer to the point object
         * @param if <code>true</code>, a copy is made, if <code>false</code> 
         *        the pointer is used as is. Defaults to <code>true</code>
         */
        void
        set(const vector<int> &gp,
            typename Point<T>::ptr p,
            bool copy=true);

        /** @return the rank (# of dimensions) of this array index
         */
        const size_t rank();

        /** @return the dimensions of this array index
         */
        const vector<size_t> &dimensions();

        /** Clear the index.
         * @param If <code>true</code> indexed points are deleted. 
         *        If <code>false</code> they are left alone.
         */
        void
        clear(bool delete_points);

#pragma mark -
#pragma mark Misc

        /** Replace the points in the given list by the points in this index. 
         * The points in the original list are released.
         */
        void
        replace_points(typename Point<T>::list &points);

        /** Iterates over the structure and counts the number of non-null
         * entries in the index.
         */
        size_t count(bool originalPointsOnly=false);
        
        /** Find the points in the neighborhood of the given gridpoint.
         * @param grid point
         * @param neighborhood size (in #grid points)
         * @return list of pointers to points (no copies)
         */
        typename Point<T>::list
        find_neighbours(const vector<int> &gridpoint, size_t reach=1);

    };
}

#endif
