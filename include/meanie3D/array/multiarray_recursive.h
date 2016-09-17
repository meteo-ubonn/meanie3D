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


#ifndef M3D_MULTIARRAY_RECURSIVE_H
#define M3D_MULTIARRAY_RECURSIVE_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <iostream>
#include <vector>
#include <netcdf>
#include <fstream>

namespace m3D {

    // Forward declarations

    using namespace std;

    template<typename T>
    class MultiArrayRecursive : public MultiArray<T>
    {
    protected:

        typedef vector<void *> array_t;
        typedef array_t *array_t_ptr;

#pragma mark -
#pragma mark Attributes

    private:

        array_t *m_data;

    private:

        void
        construct_recursive(size_t dim_index,
                            array_t **array,
                            vector<int> &gridpoint,
                            const T *default_value = NULL,
                            MultiArrayRecursive<T> *other = NULL) {
            size_t dimSize = this->m_dims[dim_index];

            if (dim_index < (this->m_dims.size() - 1)) {
                // create array

                if (dim_index == 0) {
                    vector<void *> *new_array = new vector<void *>(dimSize, NULL);

                    *array = new_array;

                    for (size_t index = 0; index < dimSize; index++) {
                        gridpoint[dim_index] = index;

                        construct_recursive(dim_index + 1, array, gridpoint, default_value, other);
                    }
                } else {
                    vector<void *> *super_array = *array;

                    size_t super_index = gridpoint[dim_index - 1];

                    vector<void *> *new_array = new vector<void *>(dimSize, NULL);

                    super_array->at(super_index) = new_array;

                    for (size_t index = 0; index < dimSize; index++) {
                        gridpoint[dim_index] = index;

                        construct_recursive(dim_index + 1, &new_array, gridpoint, default_value, other);
                    }
                }
            } else {
                vector<void *> *super_array = *array;

                if (super_array == NULL) {
                    // 1D !!
                } else {
                    size_t super_index = boost::numeric_cast<size_t>(gridpoint[dim_index - 1]);

                    vector<T> *new_array = (default_value == NULL)
                                           ? new vector<T>(dimSize)
                                           : new vector<T>(dimSize, *default_value);

                    super_array->at(super_index) = new_array;

                    if (other != NULL) {
                        vector<int> gIter = gridpoint;

                        for (size_t i = 0; i < dimSize; i++) {
                            gIter[dim_index] = i;

                            new_array->at(i) = other->get(gIter);
                        }
                    }
                }
            }
        };

        void
        destroy_recursive(size_t dim_index,
                          array_t **array,
                          vector<int> &gridpoint) {
            size_t dimSize = this->m_dims[dim_index];

            if (dim_index < (this->m_dims.size() - 1)) {
                if (dim_index == 0) {
                    vector<void *> *the_array = *array;

                    for (size_t index = 0; index < dimSize; index++) {
                        gridpoint[dim_index] = index;

                        destroy_recursive(dim_index + 1, array, gridpoint);
                    }

                    delete the_array;

                    *array = NULL;
                } else {
                    size_t super_index = gridpoint[dim_index - 1];

                    vector<void *> *super_array = *array;

                    vector<void *> *the_array = (vector<void *> *) super_array->at(super_index);

                    for (size_t index = 0; index < dimSize; index++) {
                        gridpoint[dim_index] = index;

                        destroy_recursive(dim_index + 1, &the_array, gridpoint);
                    }

                    delete the_array;

                    super_array->at(super_index) = NULL;
                }
            } else {
                vector<void *> *super_array = *array;

                size_t super_index = gridpoint[dim_index - 1];

                vector<T> *the_array = (vector<T> *) super_array->at(super_index);

                delete the_array;

                super_array->at(super_index) = NULL;
            }

        }

        void
        copy_recursive(const MultiArrayRecursive<T> *otherIndex,
                       size_t dim_index,
                       vector<int> &gridpoint) {
            size_t dimSize = this->m_dims[dim_index];

            if (dim_index < (this->m_dims.size() - 1)) {
                for (size_t index = 0; index < dimSize; index++) {
                    gridpoint[dim_index] = index;

                    copy_recursive(otherIndex, dim_index + 1, gridpoint);
                }
            } else {
                vector<int> gIter = gridpoint;

                for (size_t i = 0; i < dimSize; i++) {
                    gIter[dim_index] = i;

                    this->set(gIter, otherIndex->get(gIter));
                }
            }
        }

        void count_recursive(T value,
                             size_t &count,
                             size_t dim_index,
                             vector<int> &gp) {
            size_t dimSize = this->m_dims[dim_index];

            if (dim_index < (this->m_dims.size() - 1)) {
                for (size_t index = 0; index < dimSize; index++) {
                    gp[dim_index] = index;

                    count_recursive(value, count, dim_index + 1, gp);
                }
            } else {
                vector<int> gIter = gp;

                for (size_t i = 0; i < dimSize; i++) {
                    gIter[dim_index] = i;

                    T val = this->get(gIter);

                    if (val == value) {
                        count++;
                    }
                }
            }
        }

#pragma mark -
#pragma mark Constructors/Destructors

    public:

        MultiArrayRecursive() : MultiArray<T>() {
        };

        MultiArrayRecursive(const vector<size_t> &dims) : MultiArray<T>(dims) {
            vector<int> gp(dims.size(), 0);
            this->construct_recursive(0, &m_data, gp, NULL, NULL);
        };

        MultiArrayRecursive(const vector<size_t> &dims,
                            const T &default_value)
                : MultiArray<T>(dims, default_value) {
            vector<int> gp(dims.size(), 0);
            this->construct_recursive(0, &m_data, gp, &default_value, NULL);
        }

        MultiArrayRecursive(const MultiArrayRecursive<T> &other) {
            this->m_dims = other.get_dimensions();
            vector<int> gp(this->m_dims.size(), 0);
            this->copy_recursive(other, 0, gp);
        }

        MultiArrayRecursive<T>
        operator=(const MultiArrayRecursive<T> &other) {
            return MultiArrayRecursive(other);
        }

        /** Destructor
         */
        virtual ~MultiArrayRecursive() {
            vector<int> gp(this->m_dims.size(), 0);
            this->destroy_recursive(0, &m_data, gp);
        }

#pragma mark -
#pragma mark Accessors

        T get(const vector<int> &gp) const {
            vector<void *> *array = m_data;

            T result = 0;

            for (size_t dim_index = 0; dim_index < gp.size(); dim_index++) {
                if (dim_index < gp.size() - 1) {
                    array = (vector<void *> *) array->at(gp[dim_index]);
                } else {
                    vector<T> *points = (vector<T> *) array;

                    // Apparently the reverse lookup can hit indexes too
                    // high on occasion. This is a problem with the coordinate
                    // system class and needs to be fixed eventually

                    size_t index = gp[dim_index];

                    assert(index < points->size());

                    result = points->at(index);
                }
            }

            return result;
        }

        void
        set(const vector<int> &gp, const T &value) {
            vector<void *> *array = m_data;

            for (size_t dim_index = 0; dim_index < gp.size(); dim_index++) {
                if (dim_index < gp.size() - 1) {
                    array = (vector<void *> *) array->at(gp[dim_index]);
                } else {
                    vector<T> *points = (vector<T> *) array;

                    size_t index = gp[dim_index];

                    if (index < points->size()) {
                        points->at(index) = value;
                    }
                }
            }
        }

#pragma mark -
#pragma mark Stuff

        size_t
        count_value(const T &value) {
            size_t count = 0;
            vector<int> gp(this->m_dims.size(), 0);
            this->count_recursive(value, count, 0, gp);
            return count;
        }

        void resize(vector<size_t> dimensions) {
            throw "not implemented";
        }

        void populate_array(const T &value) {
            throw "not implemented";
        }

        void copy_from(const MultiArray <T> *other) {
            throw "not implemented";
        }

    };

}

#endif
