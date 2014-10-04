#ifndef M3D_ARRAYINDEX_IMPL_H
#define M3D_ARRAYINDEX_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>
#include <stdexcept>

#include "array_index.h"

namespace m3D {

    template <typename T>
    const size_t
    ArrayIndex<T>::rank() {
        return m_dimensions.size();
    };

    template <typename T>
    const vector<size_t> &
    ArrayIndex<T>::dimensions() {
        return m_dimensions;
    };

    template <typename T>
    ArrayIndex<T>::ArrayIndex(const vector<size_t> &dimensions,
                              bool make_copies)
    : m_dimensions(dimensions)
    , m_data(NULL)
    , m_make_copies(make_copies)
    {
        vector<int> gp = vector<int>(dimensions.size(),0);
        this->construct_array_recursive( 0, &m_data, gp );
    }

    template <typename T>
    ArrayIndex<T>::ArrayIndex(const vector<size_t> &dimensions,
                              const typename Point<T>::list &points,
                              bool make_copies)
    : m_dimensions(dimensions)
    , m_data(NULL)
    , m_make_copies(make_copies)
    {
        vector<int> gp = vector<int>(dimensions.size(),0);
        this->construct_array_recursive( 0, &m_data, gp );
        this->index(points);
    }

    template <typename T>
    ArrayIndex<T>::ArrayIndex(ArrayIndex<T> *o)
    : m_dimensions(o->m_dimensions)
    , m_data(NULL)
    , m_make_copies(o->m_make_copies)
    {
        vector<int> gp = vector<int>(o->m_dimensions.size(),0);
        this->copy_points_recursive(o,0,gp);
    }

    template <typename T>
    ArrayIndex<T>::~ArrayIndex()
    {
        vector<int> gp = vector<int>(m_dimensions.size(),0);
        this->destroy_array_recursive( 0, &m_data, gp );
    }

    template <typename T>
    void
    ArrayIndex<T>::add_points_to_list(typename Point<T>::list &points,
                                      vector<int> &gridpoint)
    {
        // obtain the size of the last dimension

        size_t dimSize = m_dimensions[gridpoint.size()-1];

        // copy all points over into the given point list
        // using the point factory's copy method

        vector<int> gIter = gridpoint;

        for (size_t i=0; i<dimSize; i++)
        {
            gIter[gridpoint.size()-1] = i;

            typename Point<T>::ptr p = this->get(gIter);

            if ( p != NULL )
            {
                typename Point<T>::ptr copy = PointFactory<T>::get_instance()->copy(p);

                points.push_back(copy);
            }
        }
    }

    template <typename T>
    void
    ArrayIndex<T>::get_points_from_index(ArrayIndex<T> *otherIndex,
                                         vector<int> &gridpoint)
    {
        // obtain the size of the last dimension

        size_t dimSize = m_dimensions[gridpoint.size()-1];

        vector<int> gIter = gridpoint;

        for (size_t i=0; i<dimSize; i++)
        {
            gIter[gridpoint.size()-1] = i;

            typename Point<T>::ptr p = otherIndex->get(gIter);

            // TODO: obtain the pointer to the array directly
            // instead of using 'set'

            if ( p != NULL )
            {
                if (this->m_make_copies)
                {
                    this->set(gIter,PointFactory<T>::get_instance()->copy(p));
                }
                else
                {
                    this->set(gIter,p);
                }
            }
        }
    }

    template <typename T>
    void
    ArrayIndex<T>::replace_points_recursive(typename Point<T>::list &points,
                                            size_t dim_index,
                                            vector<int> &gridpoint )
    {
        size_t dimSize = m_dimensions[dim_index];

        if (m_dimensions.size()==1)
        {
            this->add_points_to_list(points,gridpoint);
        }
        else if (dim_index < (m_dimensions.size()-1) )
        {
            for ( size_t index = 0; index < dimSize; index++ )
            {
                gridpoint[dim_index] = index;

                replace_points_recursive(points,dim_index+1,gridpoint);
            }
        }
        else
        {
            this->add_points_to_list(points,gridpoint);
        }
    }


    template <typename T>
    void
    ArrayIndex<T>::copy_points_recursive(ArrayIndex<T> *otherIndex,
                                         size_t dim_index,
                                         vector<int> &gridpoint)
    {
        size_t dimSize = m_dimensions[dim_index];

        if (m_dimensions.size()==1)
        {
            this->get_points_from_index(otherIndex,gridpoint);
        }
        else if (dim_index < (m_dimensions.size()-1) )
        {
            for ( size_t index = 0; index < dimSize; index++ )
            {
                gridpoint[dim_index] = index;

                copy_points_recursive(otherIndex,dim_index+1,gridpoint);
            }
        }
        else
        {
            this->get_points_from_index(otherIndex,gridpoint);
        }
    }

    template <typename T>
    void
    ArrayIndex<T>::replace_points(typename Point<T>::list &points)
    {
        // clean the original list out and
        // release all the points

        for (size_t i=0; i < points.size(); i++)
        {
            typename Point<T>::ptr p = points[i];

            delete p;

            points[i] = NULL;
        }

        points.clear();

        vector<int> gp = vector<int>(m_dimensions.size(),0);

        // Add copies of points from this index

        this->replace_points_recursive(points,0,gp);
    }

    template <typename T>
    void
    ArrayIndex<T>::construct_array_recursive(size_t dim_index,
                                             array_t **array,
                                             vector<int> &gridpoint,
                                             ArrayIndex<T> *other)
    {
        size_t dimSize = m_dimensions[dim_index];

        if (m_dimensions.size()==1)
        {
            *array = (array_t *) new vector<typename Point<T>::ptr>(dimSize,NULL);

            if (other != NULL)
            {
                this->get_points_from_index(other,gridpoint);
            }
        }
        else if (dim_index < (m_dimensions.size()-1))
        {
            // create array

            if ( dim_index == 0 )
            {
                vector<void *> *new_array = new vector<void *>(dimSize,NULL);

                *array = new_array;

                for ( size_t index = 0; index < dimSize; index++ )
                {
                    gridpoint[dim_index] = index;

                    construct_array_recursive( dim_index+1, array, gridpoint, other);
                }
            }
            else
            {
                vector<void *> *super_array = *array;

                size_t super_index = gridpoint[dim_index-1];

                vector<void *> *new_array = new vector<void *>(dimSize,NULL);

                super_array->at(super_index) = new_array;

                for ( size_t index = 0; index < dimSize; index++ )
                {
                    gridpoint[dim_index] = index;

                    construct_array_recursive( dim_index+1, &new_array, gridpoint);
                }
            }
        }
        else
        {
            vector<void *> *super_array = *array;

            size_t super_index = gridpoint[dim_index-1];

            vector<typename Point<T>::ptr> *new_array = new vector<typename Point<T>::ptr>(dimSize,NULL);

            super_array->at(super_index) = new_array;

            if ( other != NULL )
            {
                this->get_points_from_index(other,gridpoint);
            }
        }
    }

    template <typename T>
    void
    ArrayIndex<T>::destroy_array_recursive(size_t dim_index, array_t **array, vector<int> &gridpoint)
    {
        size_t dimSize = m_dimensions[dim_index];

        if (m_dimensions.size()==1)
        {
            delete m_data;
        }
        else if (dim_index < (m_dimensions.size()-1) )
        {
            if ( dim_index == 0 )
            {
                vector<void *> *the_array = *array;

                for ( size_t index = 0; index < dimSize; index++ )
                {
                    gridpoint[dim_index] = index;

                    destroy_array_recursive( dim_index+1, array, gridpoint);
                }

                delete the_array;

                *array = NULL;
            }
            else
            {
                size_t super_index = gridpoint[dim_index-1];

                vector<void *> *super_array = *array;

                vector<void *> *the_array = (vector<void *> *) super_array->at(super_index);

                for ( size_t index = 0; index < dimSize; index++ )
                {
                    gridpoint[dim_index] = index;

                    destroy_array_recursive( dim_index+1, &the_array, gridpoint);
                }

                delete the_array;

                super_array->at(super_index) = NULL;
            }
        }
        else
        {
            vector<void *> *super_array = *array;

            size_t super_index = gridpoint[dim_index-1];

            vector<typename Point<T>::ptr> *the_array = (vector<typename Point<T>::ptr> *) super_array->at(super_index);

            if (this->m_make_copies)
            {
                for (size_t i = 0; i < the_array->size(); i++)
                {
                    typename Point<T>::ptr p = the_array->at(i);

                    if (p!=NULL)
                    {
                        delete p;
                        the_array->at(i) = NULL;
                    }
                }
            }

            super_array->at(super_index) = NULL;
        }
    }

    template <typename T>
    void
    ArrayIndex<T>::index(const typename Point<T>::list &list)
    {
        for (size_t i=0; i < list.size(); i++)
        {
            typename Point<T>::ptr p = list[i];

            this->set( p->gridpoint, p, this->m_make_copies );
        }
    }

#pragma mark -
#pragma mark Accessors

    template <typename T>
    typename Point<T>::ptr
    ArrayIndex<T>::get(const vector<int> &gp)
    {
        vector<void *> *array = m_data;

        typename Point<T>::ptr result = NULL;

        for (size_t dim_index = 0; dim_index < gp.size(); dim_index++)
        {
            if ( dim_index < gp.size()-1 )
            {
                int index = gp[dim_index];

                if (index >= 0 && index < array->size())
                {
                    array = (vector<void *> *) array->at(index);
                }
                else
                {
                    throw std::invalid_argument("index parameter out of range");
                }
            }
            else
            {
                vector<typename Point<T>::ptr> *points = (vector<typename Point<T>::ptr> *) array;

                int index = gp[dim_index];

                if (index >= 0 && index < points->size())
                {
                    result = points->at(index);
                }
            }
        }

        return result;
    }

    template <typename T>
    void
    ArrayIndex<T>::set(const vector<int> &gp, typename Point<T>::ptr p, bool copy)
    {
        vector<void *> *array = m_data;

        for (size_t dim_index = 0; dim_index < gp.size(); dim_index++)
        {
            if ( dim_index < gp.size()-1 )
            {
                int index = gp[dim_index];

                if (index >= 0 && index < array->size())
                {
                    array = (vector<void *> *) array->at(index);
                }
                else
                {
                    throw std::invalid_argument("index parameter out of range");
                }
            }
            else
            {
                vector<typename Point<T>::ptr> *points = (vector<typename Point<T>::ptr> *) array;

                int index = gp[dim_index];

                if (index >= 0 && index < array->size())
                {
                    typename Point<T>::ptr existingPoint = points->at(index);

                    if (existingPoint!=NULL)
                    {
                        delete existingPoint;
                    }

                    if (copy)
                    {
                        Point<T> *c = PointFactory<T>::get_instance()->copy(p);
                        points->at(index) = c;
                        if (p->isOriginalPoint != c->isOriginalPoint)
                        {
                            cerr << "Copy ERROR" << endl;
                        }
                    }
                    else
                    {
                        points->at(index) = p;
                    }
                }
            }
        }
    }

#pragma mark -
#pragma mark Clear Index

    template <typename T>
    void
    ArrayIndex<T>::clear_recursive(size_t dim_index,
                                   array_t *array,
                                   vector<int> &gridpoint,
                                   bool delete_points)
    {
        size_t dimSize = m_dimensions[dim_index];

        if (dim_index < (m_dimensions.size()-1) )
        {
            for ( size_t index = 0; index < dimSize; index++ )
            {
                gridpoint[dim_index] = index;

                vector<void *> *a = (vector<void *> *) array->at(index);

                clear_recursive( dim_index+1, a, gridpoint, delete_points);
            }
        }
        else
        {
            vector<typename Point<T>::ptr> *points = (vector<typename Point<T>::ptr> *) array;

            for ( size_t index = 0; index < dimSize; index++ )
            {
                typename Point<T>::ptr p = points->at(index);

                if ( p != NULL && delete_points )
                {
                    delete p;
                }
            }

            points->clear();
        }
    }


    template <typename T>
    void
    ArrayIndex<T>::clear(bool delete_points)
    {
        vector<int> gp = vector<int>(m_dimensions.size(),0);

        clear_recursive(0, m_data, gp, delete_points);
    }

#pragma mark -
#pragma mark Counting

    template <typename T>
    void
    ArrayIndex<T>::count_recursive(size_t dim_index,
                                   array_t *array,
                                   vector<int> &gridpoint,
                                   size_t &count,
                                   bool originalPointsOnly)
    {
        size_t dimSize = m_dimensions[dim_index];

        if (dim_index < (m_dimensions.size()-1) )
        {
            for ( size_t index = 0; index < dimSize; index++ )
            {
                gridpoint[dim_index] = index;

                vector<void *> *a = (vector<void *> *) array->at(index);

                count_recursive( dim_index+1, a, gridpoint, count, originalPointsOnly);
            }
        }
        else
        {
            vector<typename Point<T>::ptr> *points = (vector<typename Point<T>::ptr> *) array;

            for ( size_t index = 0; index < dimSize; index++ )
            {
                typename Point<T>::ptr p = points->at(index);

                if ( p != NULL)
                {
                    if ((originalPointsOnly && p->isOriginalPoint) || !originalPointsOnly)
                    {
                        count++;
                    }
                }
            }
        }
    }


    template <typename T>
    size_t
    ArrayIndex<T>::count(bool originalPointsOnly)
    {
        size_t count = 0;

        vector<int> gp = vector<int>(m_dimensions.size(),0);

        count_recursive(0, m_data, gp, count, originalPointsOnly);

        return count;
    }
}

#endif
