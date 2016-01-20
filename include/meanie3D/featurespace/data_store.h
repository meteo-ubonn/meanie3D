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
    private:

        std::vector<std::string> m_variables;
        std::vector<std::string> m_dimensions;
        std::vector<std::string> m_dimension_variables;

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

        /**
         * Main constructor. Constructs the data store given the variables.
         * dimensions and dimension variables. 
         * 
         * @param variables
         * @param dimensions
         * @param dimension_variables
         */
        DataStore(const std::vector<std::string> &variables,
                const std::vector<std::string> &dimensions,
                const std::vector<std::string> &dimension_variables)
        : m_variables(variables)
        , m_dimensions(dimensions)
        , m_dimension_variables(dimension_variables) {};

        /** Destructor
         */
        virtual ~DataStore() {};

#pragma mark -
#pragma mark Public methods

        virtual 
        const std::vector<std::string> &variables() const = 0;

        /**
         * Implementing classes should know how to construct a 
         * coordinate system from the dimensions and dimension
         * variables used to configure the data store.
         * 
         * @param variable_index
         * @param index
         * @param is_valid
         * @return 
         */
        virtual 
        CoordinateSystem<T> *coordinate_system() const = 0;

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

        /** 
         * @return number of variables in the data store
         */
        virtual
        const size_t rank() const = 0;

        /** @return number of total points in the grid. 
         * (basically the product of dimension sizes)
         */
        virtual
        const size_t size() const = 0;

        /** 
         * @param index variable index
         * @return min encountered value of a variable.
         */
        virtual
        const T min(size_t index) const = 0;

        /** 
         * @param index variable index
         * @return max encountered value of a variable.
         */
        virtual
        const T max(size_t index) const = 0;

        /** Iterates over each point in the data store 
         * and calls the given callback handler.
         */
        virtual
        void
        for_each(size_t variable_index, ForEachFunctor *f) = 0;

        // META data

        /**
         * Lower valid bound of variable values.
         * @param index variable index
         * @return 
         */
        virtual
        const T valid_min(size_t index) const = 0;

        /**
         * Upper valid bound of variable values.
         * @param index
         * @return 
         */
        virtual
        const T valid_max(size_t index) const = 0;

        /**
         * Scale factor for unpacking packed data.
         * 
         * @param index
         * @return 
         */
        virtual
        const T scale_factor(size_t index) const = 0;

        /**
         * Offset for unpacking packed data.
         * 
         * @param index
         * @return 
         */
        virtual
        const T add_offset(size_t index) const = 0;

        /**
         * Value representing missing value. 
         * 
         * @param index
         * @return 
         */
        virtual
        const T fill_value(size_t index) const = 0;

    };
}

#endif
