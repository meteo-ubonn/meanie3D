/* The MIT License (MIT)
 * 
 * (c) Jürgen Simon 2014 (juergen.simon@uni-bonn.de)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#ifndef M3D_POINT_FACTORY_H
#define M3D_POINT_FACTORY_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <meanie3D/featurespace/point_default_factory.h>

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
        create(vector<int> &gridpoint, vector<T> &coord, vector<T> &value) = 0;

        /** This is an abstract method used to create new instances of Point.
         * @param coordinates
         * @param value
         * @deprecated
         */
        virtual
        Point<T> *
        create(vector<T> &coord, vector<T> &value) = 0;


        /** Copy 
         * @param original
         * @return copy
         */
        virtual
        Point<T> *
        copy(const Point<T> *p) = 0;
    };

    // initialize with default factory

    template <typename T>
    PointFactory<T>* PointFactory<T>::m_instance_ptr = new PointDefaultFactory<T>();
}

#endif
