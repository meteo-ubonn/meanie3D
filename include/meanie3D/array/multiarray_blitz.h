#ifndef M3D_MULTIARRAY_BLITZ_H
#define M3D_MULTIARRAY_BLITZ_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <meanie3D/array/multiarray.h>

#include <blitz/array.h>

#include <vector>
#include <netcdf>

namespace m3D { 

#pragma mark -
#pragma mark Helper macros

/* Indexing */

#define A1(dim) dim[0]
#define A2(dim) dim[0],dim[1]
#define A3(dim) dim[0],dim[1],dim[2]
#define A4(dim) dim[0],dim[1],dim[2],dim[3]
#define A5(dim) dim[0],dim[1],dim[2],dim[3],dim[4]

/* Filling */

#define FILL_A1(dim,array,value) \
for(int i1=0;i1<dim[0];i1++) \
    array(i1)=value

#define FILL_A2(dim,array,value) \
for(int i1=0;i1<dim[0];i1++) \
    for(int i2=0;i2<dim[1];i2++) \
        array(i1,i2)=value

#define FILL_A3(dim,array,value) \
for(int i1=0;i1<dim[0];i1++) \
    for(int i2=0;i2<dim[1];i2++) \
        for(int i3=0;i3<dim[2];i3++) \
            array(i1,i2,i3)=value

#define FILL_A4(dim,array,value) \
for(int i1=0;i1<dim[0];i1++) \
    for(int i2=0;i2<dim[1];i2++) \
        for(int i3=0;i3<dim[2];i3++) \
            for(int i4=0;i4<dim[3];i4++) \
                array(i1,i2,i3,i4)=value

#define FILL_A5(dim,array,value) \
for(int i1=0;i1<dim[0];i1++) \
    for(int i2=0;i2<dim[1];i2++) \
        for(int i3=0;i3<dim[2];i3++) \
            for(int i4=0;i4<dim[3];i4++) \
                for(int i5=0;i5<dim[4];i5++) \
                    array(i1,i2,i3,i4,i5)=value

    /* Copying */

#define COPY_A1(dim,array,value) \
for(int i1=0;i1<dim[0];i1++){ \
    vector<int> idx(1); \
    idx[0]=i1; \
    m_a1(i1)=array->get(idx);}

#define COPY_A2(dim,array,value) \
for(int i1=0;i1<dim[0];i1++){ \
    for(int i2=0;i2<dim[1];i2++){ \
        vector<int> idx(2); \
        idx[0]=i1; idx[1]=i2; \
        m_a2(i1,i2)=array->get(idx);}}

#define COPY_A3(dim,array,value) \
for(int i1=0;i1<dim[0];i1++){ \
    for(int i2=0;i2<dim[1];i2++){ \
        for(int i3=0;i3<dim[2];i3++){ \
            vector<int> idx(3); \
            idx[0]=i1; idx[1]=i2; idx[2]=i3; \
            m_a3(i1,i2,i3)=array->get(idx);}}}

#define COPY_A4(dim,array,value) \
for(int i1=0;i1<dim[0];i1++){ \
    for(int i2=0;i2<dim[1];i2++){ \
        for(int i3=0;i3<dim[2];i3++){ \
            for(int i4=0;i4<dim[3];i4++){ \
                vector<int> idx(4); \
                idx[0]=i1; idx[1]=i2; idx[2]=i3; idx[3]=i4; \
                m_a4(i1,i2,i3,i4)=array->get(idx);}}}}

#define COPY_A5(dim,array,value) \
for(int i1=0;i1<dim[0];i1++){ \
    for(int i2=0;i2<dim[1];i2++){ \
        for(int i3=0;i3<dim[2];i3++){ \
            for(int i4=0;i4<dim[3];i4++){ \
                for(int i5=0;i5<dim[4];i5++){ \
                    vector<int> idx(5); \
                    idx[0]=i1; idx[1]=i2; idx[2]=i3; idx[3]=i4; idx[4]=i5; \
                    m_a5(i1,i2,i3,i4,i5)=array->get(idx);}}}}}

    template <typename T>
    class MultiArrayBlitz : public MultiArray<T>
    {

#pragma mark -
#pragma mark Attributes

    private:

        // only 1D .. 10D is currently supported.
        // TODO: think about employing template
        // metaprogramming for this

        blitz::Array<T,1> m_a1;
        blitz::Array<T,2> m_a2;
        blitz::Array<T,3> m_a3;
        blitz::Array<T,4> m_a4;
        blitz::Array<T,5> m_a5;

        void destroy();

#pragma mark -
#pragma mark Constructors/Destructors

    public:

        MultiArrayBlitz() : MultiArray<T>() {};

        MultiArrayBlitz(const vector<size_t> &dims)
            : MultiArray<T>(dims)
        {
            switch (this->m_dims.size())
            {
                case 1: m_a1.resize(A1(this->m_dims)); break;
                case 2: m_a2.resize(A2(this->m_dims)); break;
                case 3: m_a3.resize(A3(this->m_dims)); break;
                case 4: m_a4.resize(A4(this->m_dims)); break;
                case 5: m_a5.resize(A5(this->m_dims)); break;
                default: throw std::out_of_range("only dimensions 1 .. 5 are supported");
            }
        };

        MultiArrayBlitz(const vector<size_t> &dims, T default_value)
            : MultiArray<T>(dims,default_value)
        {
            switch (this->m_dims.size())
            {
                case 1: {
                    m_a1.resize(A1(this->m_dims));
                    FILL_A1(this->m_dims, m_a1, default_value);
                } break;

                case 2: {
                    m_a2.resize(A2(this->m_dims));
                    FILL_A2(this->m_dims, m_a2, default_value);
                } break;

                case 3: {
                    m_a3.resize(A3(this->m_dims));
                    FILL_A3(this->m_dims, m_a3, default_value);
                } break;

                case 4: {
                    m_a4.resize(A4(this->m_dims));
                    FILL_A4(this->m_dims, m_a4, default_value);
                } break;

                case 5: {
                    m_a5.resize(A5(this->m_dims));
                    FILL_A5(this->m_dims, m_a5, default_value);
                } break;

                default:
                    throw std::out_of_range("only 5 dimensions are currently supported");
            }
        };

        /** Copy constructor
         */
        MultiArrayBlitz(const MultiArrayBlitz<T> &other) : MultiArray<T>(other)
        {
            switch (this->m_dims.size())
            {
                case 1: m_a1 = blitz::Array<T,1>(other.m_a1); break;
                case 2: m_a2 = blitz::Array<T,2>(other.m_a2); break;
                case 3: m_a3 = blitz::Array<T,3>(other.m_a3); break;
                case 4: m_a4 = blitz::Array<T,4>(other.m_a4); break;
                case 5: m_a5 = blitz::Array<T,5>(other.m_a5);
                default: throw std::out_of_range("only 5 dimensions are currently supported");
            }
        }

        /** Copy constructor from pointer 
         */
        MultiArrayBlitz(const MultiArrayBlitz<T> *other) : MultiArray<T>(other)
        {
            switch (this->m_dims.size())
            {
                case 1: m_a1 = blitz::Array<T,1>(other->m_a1); break;
                case 2: m_a2 = blitz::Array<T,2>(other->m_a2); break;
                case 3: m_a3 = blitz::Array<T,3>(other->m_a3); break;
                case 4: m_a4 = blitz::Array<T,4>(other->m_a4); break;
                case 5: m_a5 = blitz::Array<T,5>(other->m_a5);
                default: throw std::out_of_range("only 5 dimensions are currently supported");
            }
        }

        MultiArrayBlitz<T>
        operator = (const MultiArrayBlitz& other)
        {
            return MultiArrayBlitz<T>(other);
        }

        /** Destructor
         */
        ~MultiArrayBlitz() {};

#pragma mark -
#pragma mark Accessors

        T get(const vector<int> &index) const
        {
            T result;
            switch (this->m_dims.size())
            {
                case 1: result = m_a1(A1(index)); break;
                case 2: result = m_a2(A2(index)); break;
                case 3: result = m_a3(A3(index)); break;
                case 4: result = m_a4(A4(index)); break;
                case 5: result = m_a5(A5(index)); break;
                default: throw std::out_of_range("only 5 dimensions are currently supported");
            }

            return result;
        }

        void
        set(const vector<int> &index, const T &value)
        {

            switch (this->m_dims.size())
            {
                case 1: m_a1(A1(index)) = value; break;
                case 2: m_a2(A2(index)) = value; break;
                case 3: m_a3(A3(index)) = value; break;
                case 4: m_a4(A4(index)) = value; break;
                case 5: m_a5(A5(index)) = value; break;
                default: throw std::out_of_range("only 5 dimensions are currently supported");
            }
        }

#pragma mark -
#pragma mark Stuff

        void resize(vector<size_t> dimensions)
        {
            if (dimensions.size() != this->m_dims.size())
            {
                // release old data

                vector<int> empty(this->m_dims.size(),0);

                switch (dimensions.size())
                {
                    case 1: m_a1.resize(0); break;
                    case 2: m_a2.resize(0,0); break;
                    case 3: m_a3.resize(0,0,0); break;
                    case 4: m_a4.resize(0,0,0,0); break;
                    case 5: m_a5.resize(0,0,0,0,0); break;
                    default: throw std::out_of_range("only 5 dimensions are currently supported");
                }
            }

            this->m_dims = dimensions;

            // Re-allocate

            switch (this->m_dims.size())
            {
                case 1: m_a1.resize(A1(this->m_dims)); break;
                case 2: m_a2.resize(A2(this->m_dims)); break;
                case 3: m_a3.resize(A3(this->m_dims)); break;
                case 4: m_a4.resize(A4(this->m_dims)); break;
                case 5: m_a5.resize(A5(this->m_dims)); break;
                default: throw std::out_of_range("only 5 dimensions are currently supported");
            }
        }

        void populate_array(const T& value)
        {
            switch (this->m_dims.size())
            {
                case 1: FILL_A1(this->m_dims, m_a1, value); break;
                case 2: FILL_A2(this->m_dims, m_a2, value); break;
                case 3: FILL_A3(this->m_dims, m_a3, value); break;
                case 4: FILL_A4(this->m_dims, m_a4, value); break;
                case 5: FILL_A5(this->m_dims, m_a5, value); break;
                default: throw std::out_of_range("only 5 dimensions are currently supported");
            }
        }

        void copy_from(const MultiArray<T> *other)
        {
            assert(this->m_dims == other->get_dimensions());

            switch (this->m_dims.size())
            {
                case 1: COPY_A1(this->m_dims, other, default_value); break;
                case 2: COPY_A2(this->m_dims, other, default_value); break;
                case 3: COPY_A3(this->m_dims, other, default_value); break;
                case 4: COPY_A4(this->m_dims, other, default_value); break;
                case 5: COPY_A5(this->m_dims, other, default_value); break;
                default: throw std::out_of_range("only 5 dimensions are currently supported");
            }
        }

        /** Iterates over the array and counts the number
         * of occurences of the given value
         * @param value
         * @return number of occurences
         */
        size_t count_value(const T &value)
        {
            throw "not implemented";
        }

    };

}

#endif
