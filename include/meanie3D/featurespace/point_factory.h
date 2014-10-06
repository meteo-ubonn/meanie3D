#ifndef M3D_POINT_FACTORY_H
#define M3D_POINT_FACTORY_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

namespace m3D { 

    /** Abstract factory pattern coupled with a singleton pattern
     * to solve the problem of feature-space being able to produce
     * different types of point objects without having to know
     * about them.
     *
     * TODO: consider using boost::shared_ptr for storing the instance pointer
     */
    template <typename T>
    class PointFactory
    {
    private:

        // concrete factory instance
        static PointFactory<T> * m_instance_ptr;

    public:

        /** This method sets the concrete instance of PointFactory to be used.
         * @param pointer to point factory
         */
        static void set_instance(PointFactory<T> *factory_ptr)
        {
            m_instance_ptr = factory_ptr;
        };

        /** This method gets the concrete instance of PointFactory. If no instance
         * is set when this is called, the default is returned.
         * @returns pointer to point factory
         */
        static PointFactory<T> * get_instance()
        {
            return m_instance_ptr;
        };

        /** This is an abstract method used to create new instances of Point.
         */
        virtual
        Point<T> *
        create() = 0;

        /** This is an abstract method used to create new instances of Point.
         * @param gridpoint
         * @param coordinates
         * @param value
         */
        virtual 
        Point<T> * 
        create( vector<int> &gridpoint, vector<T> &coord, vector<T> &value ) = 0;

        /** This is an abstract method used to create new instances of Point.
         * @param coordinates
         * @param value
         * @deprecated
         */
        virtual
        Point<T> *
        create( vector<T> &coord, vector<T> &value ) = 0;


        /** Copy 
         * @param original
         * @return copy
         */
        virtual
        Point<T> *
        copy( const Point<T> *p) = 0;
    };

    // initialize with default factory

    template <typename T>
    PointFactory<T>* PointFactory<T>::m_instance_ptr = new PointDefaultFactory<T>();
}

#endif
