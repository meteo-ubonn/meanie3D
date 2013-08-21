#ifndef _M3D_ArrayIndex_Impl_H_
#define _M3D_ArrayIndex_Impl_H_

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include <meanie3D/types/point.h>

namespace m3D {

	using namespace std;
	using cfa::meanshift::Point;
    
    template <typename T>
    ArrayIndex<T>::ArrayIndex(CoordinateSystem<T> *cs,bool make_copies)
    : m_coordinate_system(cs)
    , m_data(NULL)
    , m_make_copies(make_copies)
    {
        typename CoordinateSystem<T>::GridPoint gp = m_coordinate_system->newGridPoint();
        
        this->construct_array_recursive( 0, &m_data, gp );
    }
    
    template <typename T>
    ArrayIndex<T>::ArrayIndex(CoordinateSystem<T> *cs,
                              const typename Point<T>::list &points,
                              bool make_copies)
    : m_coordinate_system(cs)
    , m_data(NULL)
    , m_make_copies(make_copies)
    {
        typename CoordinateSystem<T>::GridPoint gp = m_coordinate_system->newGridPoint();
        
        this->construct_array_recursive( 0, &m_data, gp );
        
        this->index(points);
    }
    
    template <typename T>
    ArrayIndex<T>::ArrayIndex(ArrayIndex<T> *o)
    : m_coordinate_system(o->m_coordinate_system)
    , m_data(NULL)
    , m_make_copies(o->m_make_copies)
    {
        typename CoordinateSystem<T>::GridPoint gp = m_coordinate_system->newGridPoint();

        this->copy_points_recursive(o,0,gp);
    }
    
    template <typename T>
    ArrayIndex<T>::~ArrayIndex()
    {
        typename CoordinateSystem<T>::GridPoint gp = m_coordinate_system->newGridPoint();

        this->destroy_array_recursive( 0, &m_data, gp );
    }
    
    template <typename T>
    void
    ArrayIndex<T>::replace_points_recursive(typename Point<T>::list &points,
                                            size_t dim_index,
                                            typename CoordinateSystem<T>::GridPoint &gridpoint )
    {
        NcDim dim = m_coordinate_system->dimensions()[dim_index];
        
        size_t dimSize = dim.getSize();
        
        if (dim_index < (m_coordinate_system->size()-1) )
        {
            for ( size_t index = 0; index < dimSize; index++ )
            {
                gridpoint[dim_index] = index;
                
                replace_points_recursive(points,dim_index+1,gridpoint);
            }
        }
        else
        {
            // Create a copy, since the original will be released
            // with the working copy
            
            typename CoordinateSystem<T>::GridPoint gIter = gridpoint;
            
            for (size_t i=0; i<dimSize; i++)
            {
                gIter[dim_index] = i;
                
                typename Point<T>::ptr p = this->get(gIter);
                
                if ( p != NULL )
                {
                    typename Point<T>::ptr copy = PointFactory<T>::get_instance()->copy(p);
                    
                    points.push_back(copy);
                }
            }
        }
    }
    
    template <typename T>
    void
    ArrayIndex<T>::copy_points_recursive(ArrayIndex<T> *otherIndex,
                                         size_t dim_index,
                                         typename CoordinateSystem<T>::GridPoint &gridpoint)
    {
        NcDim dim = m_coordinate_system->dimensions()[dim_index];
        
        size_t dimSize = dim.getSize();
        
        if (dim_index < (m_coordinate_system->size()-1) )
        {
            for ( size_t index = 0; index < dimSize; index++ )
            {
                gridpoint[dim_index] = index;
                
                copy_points_recursive(otherIndex,dim_index+1,gridpoint);
            }
        }
        else
        {
            typename CoordinateSystem<T>::GridPoint gIter = gridpoint;
            
            for (size_t i=0; i<dimSize; i++)
            {
                gIter[dim_index] = i;
                
                typename Point<T>::ptr p = this->otherIndex(gIter);
                
                if ( p != NULL )
                {
                    this->set(gIter,p,true);
                }
            }
        }
    }

    template <typename T>
    void
    ArrayIndex<T>::replace_points(typename Point<T>::list &points)
    {
        // clean the list out
        
        for (size_t i=0; i < points.size(); i++)
        {
            typename Point<T>::ptr p = points[i];
            
            delete p;
            
            points[i] = NULL;
        }
        
        points.clear();

        typename CoordinateSystem<T>::GridPoint gp = m_coordinate_system->newGridPoint();
        
        this->replace_points_recursive(points,0,gp);
    }

    template <typename T>
    void
    ArrayIndex<T>::construct_array_recursive(size_t dim_index,
                                             array_t **array,
                                             typename CoordinateSystem<T>::GridPoint &gridpoint,
                                             ArrayIndex<T> *other)
    {
        NcDim dim = m_coordinate_system->dimensions()[dim_index];
        
        size_t dimSize = dim.getSize();
        
        if (dim_index < (m_coordinate_system->size()-1) )
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
                typename CoordinateSystem<T>::GridPoint gIter = gridpoint;
                
                for (size_t i=0; i<dimSize; i++)
                {
                    gIter[dim_index] = i;
                    
                    typename Point<T>::ptr p = other->get(gIter);
                    
                    if ( p != NULL )
                    {
                        if (this->m_make_copies)
                        {
                            new_array->at(i) = PointFactory<T>::get_instance()->copy(p);
                        }
                        else
                        {
                            new_array->at(i) = p;
                        }
                    }
                }
            }
        }
    }
    
    template <typename T>
    void
    ArrayIndex<T>::destroy_array_recursive(size_t dim_index, array_t **array, typename CoordinateSystem<T>::GridPoint &gridpoint)
    {
        NcDim dim = m_coordinate_system->dimensions()[dim_index];
        
        size_t dimSize = dim.getSize();
        
        if (dim_index < (m_coordinate_system->size()-1) )
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
    ArrayIndex<T>::get(const typename CoordinateSystem<T>::GridPoint &gp)
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
            }
            else
            {
                vector<typename Point<T>::ptr> *points = (vector<typename Point<T>::ptr> *) array;
                
                int index = gp[dim_index];
                
                if (index >= 0 && index < array->size())
                {
                    result = points->at(index);
                }
            }
        }
        
        return result;
    }
    
    template <typename T>
    void
    ArrayIndex<T>::set(const typename CoordinateSystem<T>::GridPoint &gp, typename Point<T>::ptr p, bool copy)
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
                                   typename CoordinateSystem<T>::GridPoint &gridpoint,
                                   bool delete_points)
    {
        NcDim dim = m_coordinate_system->dimensions()[dim_index];
        
        size_t dimSize = dim.getSize();
        
        if (dim_index < (m_coordinate_system->size()-1) )
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
        typename CoordinateSystem<T>::GridPoint gp = m_coordinate_system->newGridPoint();
        
        clear_recursive(0, m_data, gp, delete_points);
    }
    
#pragma mark -
#pragma mark Counting
    
    template <typename T>
    void
    ArrayIndex<T>::count_recursive(size_t dim_index,
                                   array_t *array,
                                   typename CoordinateSystem<T>::GridPoint &gridpoint,
                                   size_t &count,
                                   bool originalPointsOnly)
    {
        NcDim dim = m_coordinate_system->dimensions()[dim_index];
        
        size_t dimSize = dim.getSize();

        if (dim_index < (m_coordinate_system->size()-1) )
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
        
        typename CoordinateSystem<T>::GridPoint gp = m_coordinate_system->newGridPoint();
        
        count_recursive(0, m_data, gp, count, originalPointsOnly);
        
        return count;
    }

    
}; //namespace

#endif
