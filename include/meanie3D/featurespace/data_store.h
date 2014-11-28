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


#ifndef M3D_DATASTORE_H
#define M3D_DATASTORE_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/parallel.h>
#include <meanie3D/array.h>
#include <meanie3D/utils.h>

#include <boost/progress.hpp>

#include <map>
#include <vector>
#include <netcdf>

namespace m3D { 

    /** The DataStore represents a set of input data for the method.
     */
    template <typename T>
    class DataStore
    {
        
    public:
        
#pragma mark -
#pragma mark Typedefs
        
        typedef DataStore<T> * ptr;

#pragma mark -
#pragma mark for_each functor

        class ForEachFunctor
        {
        public:
            virtual
            void
            operator()(DataStore<T> *store,
                       const size_t variable_index,
                       vector<int> &gridpoint,
                       T& value,
                       bool &is_valid) = 0;
        };

#pragma mark -
#pragma mark Consctructor/Destructor

        DataStore() {};

        /** Destructor
         */

        virtual ~DataStore() {};

#pragma mark -
#pragma mark Public methods

        /** Retrieves a value.
         *
         * @param index of the variable
         * @param grid point
         * @param reference to bool variable, which will indicate if the
         *        value is within valid range or not after the call
         * @return (<b>unpacked</b>) variable value
         */
        virtual
        T get(size_t variable_index,
              const vector<int> &gridpoint,
              bool &is_valid) const = 0;

        /** Gets a point by it's linear index. The index may run
         * from 0 ... (N-1) where N is the total number of points
         * in the grid.
         */
        virtual 
        T get(size_t variable_index,
              size_t index,
              bool &is_valid) const = 0;

        /** @return number of variables in the data store
         */
        virtual
        const size_t rank() const = 0;

        /** @return number of total points in the grid. 
         * (basically the product of dimension sizes)
         */
        virtual 
        const size_t size() const = 0;

        /** @return extent in the dimensions
         */
        virtual
        const vector<size_t> get_dimension_sizes() const = 0;

        /** @return min value of a variable
         */
        virtual
        const T min(size_t index) const = 0;

        /** @return max value of a variable
         */
        virtual
        const T max(size_t index) const = 0;

        /** Iterates over each point in the data store 
         * and calls the given callback handler.
         */
        virtual
        void
        for_each(size_t variable_index, ForEachFunctor *f) = 0;
    };
}

#endif
