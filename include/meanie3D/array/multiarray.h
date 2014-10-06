#ifndef M3D_MULTI_ARRAY_H
#define M3D_MULTI_ARRAY_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils/vector_utils.h>

#include <vector>

namespace m3D { 

    template <typename T>
    class MultiArray
    {

#pragma mark -
#pragma mark Attributes

    protected:

        vector<size_t>  m_dims;

#pragma mark -
#pragma mark Constructors/Destructors

    public:

        class ForEachFunctor
        {
            public:
            virtual void operator()(const vector<int> &index, const T value) = 0;
        };

        /** Constructs an empty mutli-dimensional array
         */
        MultiArray()
            : m_dims(vector<size_t>()) {};

        /** Constructs an empty mutli-dimensional array
         * with as many dimensions as the dimension vector
         * indicates. Example: [4,5] constructs a 4x5 array.
         * @param dimensions
         */
        MultiArray(const vector<size_t> &dims)
            : m_dims(dims) {};

        /** Constructs an empty mutli-dimensional array
         * with as many dimensions as the dimension vector
         * indicates and populates it with a default value.
         * @param dimensions
         * @param default value
         */
        MultiArray(const vector<size_t> &dims, T default_value)
            : m_dims(dims)
        {};

        /** Destructor
         */
        virtual ~MultiArray() {};

        /** Copy constructor.
         * @param other index
         */
        MultiArray(const MultiArray<T> &other)
            : m_dims(other.get_dimensions())
        {
            this->copy_from(&other);
        };

        /** Copy constructor on pointer 
         */
        MultiArray(const MultiArray<T> *other)
        : m_dims(other->get_dimensions())
        {
            this->copy_from(other);
        };

        /** Copy operator
         * @param other index
         */
        MultiArray<T> *
        operator = (const MultiArray* other)
        {
            this->m_dims = other->get_dimensions();
            this->copy_from(other);
        }

#pragma mark -
#pragma mark Accessors

        /** Gets value at given gridpoint
         * @param grid point
         */
        virtual
        T get(const vector<int> &index) const = 0;

        /** Sets a value in the index at a given grid point.
         * @param grid point
         * @param value
         */
        virtual
        void set(const vector<int> &index, const T &value) = 0;

        /** @return const reference to the dimension vector
         * this array was build on
         */
        const vector<size_t> & get_dimensions() const { return m_dims; };

        /** @return number of points in this array
         */
        const size_t size() const {
            size_t num_values = 1;
            for (int i=0; i<m_dims.size(); i++) num_values *= m_dims[i];
            return num_values;
        }

        /** @return the array's dimensionality
         */
        const size_t rank() const { return m_dims.size(); }

#pragma mark -
#pragma mark Stuff

        /** Resizes the array. All data will be lost.
         */
        virtual
        void resize(vector<size_t> dimensions) = 0;

        /** Fills the whole array with the given value.
         * @param value
         */
        virtual
        void populate_array(const T& value) = 0;

        /** Copy the entire data over from another array
         * @param other multiarray
         */
        virtual
        void copy_from(const MultiArray<T> *other) = 0;

        /** Iterates over the array and counts the number
         * of occurences of the given value
         * @param value
         * @return number of occurences
         */
        virtual
        size_t count_value(const T &value) = 0;
        
#pragma mark -
#pragma mark I want variants ...

        /**
         * @param fileName
         * @param variableName
         * TODO: implement generically using for_each 
         */
        virtual
        void write(const std::string &fileName,const std::string &variableName) const = 0;
        
#pragma mark -
#pragma mark For each

        /** Iterates over each point in the array and calls the given function.
         * @param function to call
         */
        void for_each(ForEachFunctor *f)
        {
            vector<int> index(this->get_dimensions().size(),0);
            this->for_each_recursive(f, 0, index);
        }

    private:

        void for_each_recursive(ForEachFunctor *f, size_t dim_index, vector<int> &index)
        {
            size_t dimSize = this->get_dimensions()[dim_index];

            if (dim_index < (this->get_dimensions().size()-1) )
            {
                for ( size_t i = 0; i < dimSize; i++ )
                {
                    index[dim_index] = i;

                    for_each_recursive(f,dim_index+1,index);
                }
            }
            else
            {
                vector<int> gIter = index;

                for (size_t i=0; i<dimSize; i++)
                {
                    gIter[dim_index] = i;

                    f->operator()(gIter, this->get(gIter));
                }
            }
        }

#pragma mark -
#pragma mark Range search

    private:

        void
        add_values_around_recursive(const vector<int> lower_bounds,
                                    const vector<int> upper_bounds,
                                    vector<int> &gridpoint,
                                    size_t dim_index,
                                    vector<T> &values)
        {
            if (dim_index < (this->get_dimensions().size()-1) )
            {
                for ( size_t i = lower_bounds[dim_index]; i <= upper_bounds[dim_index]; i++ )
                {
                    gridpoint[dim_index] = i;

                    add_values_around_recursive(lower_bounds,upper_bounds,gridpoint,dim_index+1,values);
                }
            }
            else
            {
                vector<int> gIter = gridpoint;

                for ( size_t i = lower_bounds[dim_index]; i <= upper_bounds[dim_index]; i++ )
                {
                    gIter[dim_index] = i;

                    values.push_back(this->get(gIter));
                }
            }
        }

    public:

        void
        values_around(const vector<int> gridpoint, const vector<int> &bandwidth, vector<T> &result)
        {
            // set search ranges and make sure they are within bounds
            
            using namespace utils::vectors;

            vector<int> lower_bounds = gridpoint - bandwidth;
            vector<int> upper_bounds = gridpoint + bandwidth;

            for (size_t dim_index=0; dim_index < (this->get_dimensions().size()-1); dim_index++ )
            {
                if (lower_bounds[dim_index] < 0)
                {
                    lower_bounds[dim_index] = 0;
                }
                if (upper_bounds[dim_index] > (this->get_dimensions()[dim_index]-1))
                {
                    upper_bounds[dim_index] = (this->get_dimensions()[dim_index]-1);
                }
            }

            vector<int> gp = lower_bounds;

            this->add_values_around_recursive(lower_bounds,upper_bounds,gp,0,result);
        }
    };
    
}

#endif
