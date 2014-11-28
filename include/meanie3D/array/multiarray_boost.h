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


#ifndef M3D_MULTIARRAY_BOOST_H
#define M3D_MULTIARRAY_BOOST_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <boost/multi_array.hpp>

#include <vector>

namespace m3D { 

#pragma mark -
#pragma mark Helper macros

/* Indexing */

#define BOOST_A1(e,dim) e[dim[0]]
#define BOOST_A2(e,dim) e[dim[0]][dim[1]]
#define BOOST_A3(e,dim) e[dim[0]][dim[1]][dim[2]]
#define BOOST_A4(e,dim) e[dim[0]][dim[1]][dim[2]][dim[3]]
#define BOOST_A5(e,dim) e[dim[0]][dim[1]][dim[2]][dim[3]][dim[4]]

#define BOOST_ID1(a,idx) a[idx[0]]
#define BOOST_ID2(a,idx) a[idx[0]][idx[1]]
#define BOOST_ID3(a,idx) a[idx[0]][idx[1]][idx[2]]
#define BOOST_ID4(a,idx) a[idx[0]][idx[1]][idx[2]][idx[3]]
#define BOOST_ID5(a,idx) a[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]

/* Filling */

#define BOOST_FILL_A1(dim,array,value) \
for(int i1=0;i1<dim[0];i1++) \
    array[i1]=value

#define BOOST_FILL_A2(dim,array,value) \
for(int i1=0;i1<dim[0];i1++) \
    for(int i2=0;i2<dim[1];i2++) \
        array[i1][i2]=value

#define BOOST_FILL_A3(dim,array,value) \
for(int i1=0;i1<dim[0];i1++) \
    for(int i2=0;i2<dim[1];i2++) \
        for(int i3=0;i3<dim[2];i3++) \
            array[i1][i2][i3]=value

#define BOOST_FILL_A4(dim,array,value) \
for(int i1=0;i1<dim[0];i1++) \
    for(int i2=0;i2<dim[1];i2++) \
        for(int i3=0;i3<dim[2];i3++) \
            for(int i4=0;i4<dim[3];i4++) \
                array[i1][i2][i3][i4]=value

#define BOOST_FILL_A5(dim,array,value) \
for(int i1=0;i1<dim[0];i1++) \
    for(int i2=0;i2<dim[1];i2++) \
        for(int i3=0;i3<dim[2];i3++) \
            for(int i4=0;i4<dim[3];i4++) \
                for(int i5=0;i5<dim[4];i5++) \
                    array[i1][i2][i3][i4][i5]=value


/* Copying */

#define BOOST_COPY_A1(dim,array,value) \
for(int i1=0;i1<dim[0];i1++){ \
    vector<int> idx(1); \
    idx[0]=i1; \
    m_a1[i1]=array->get(idx);}

#define BOOST_COPY_A2(dim,array,value) \
for(int i1=0;i1<dim[0];i1++){ \
    for(int i2=0;i2<dim[1];i2++){ \
        vector<int> idx(2); \
        idx[0]=i1; idx[1]=i2; \
        m_a2[i1][i2]=array->get(idx);}}

#define BOOST_COPY_A3(dim,array,value) \
for(int i1=0;i1<dim[0];i1++){ \
    for(int i2=0;i2<dim[1];i2++){ \
        for(int i3=0;i3<dim[2];i3++){ \
            vector<int> idx(3); \
            idx[0]=i1; idx[1]=i2; idx[2]=i3; \
            m_a3[i1][i2][i3]=array->get(idx);}}}

#define BOOST_COPY_A4(dim,array,value) \
for(int i1=0;i1<dim[0];i1++){ \
    for(int i2=0;i2<dim[1];i2++){ \
        for(int i3=0;i3<dim[2];i3++){ \
            for(int i4=0;i4<dim[3];i4++){ \
                vector<int> idx(4); \
                idx[0]=i1; idx[1]=i2; idx[2]=i3; idx[3]=i4; \
                m_a4[i1][i2][i3][i4]=array->get(idx);}}}}

#define BOOST_COPY_A5(dim,array,value) \
for(int i1=0;i1<dim[0];i1++){ \
    for(int i2=0;i2<dim[1];i2++){ \
        for(int i3=0;i3<dim[2];i3++){ \
            for(int i4=0;i4<dim[3];i4++){ \
                for(int i5=0;i5<dim[4];i5++){ \
                    vector<int> idx(5); \
                    idx[0]=i1; idx[1]=i2; idx[2]=i3; idx[3]=i4; idx[4]=i5; \
                    m_a5[i1][i2][i3][i4][i5]=array->get(idx);}}}}}

    template <typename T>
    class MultiArrayBoost : public MultiArray<T>
    {

#pragma mark -
#pragma mark Attributes

    private:

        // only 1D .. 10D is currently supported.
        // TODO: think about employing template
        // metaprogramming for this

        boost::multi_array<T,1> m_a1;
        boost::multi_array<T,2> m_a2;
        boost::multi_array<T,3> m_a3;
        boost::multi_array<T,4> m_a4;
        boost::multi_array<T,5> m_a5;

    private:

        void destroy();

        void write(const std::string &fileName,
                const std::string &variableName) const 
        {
            throw "not_implemented";
        }

#pragma mark -
#pragma mark Constructors/Destructors

    public:

        MultiArrayBoost() : MultiArray<T>() {};

        MultiArrayBoost(const vector<size_t> &dims)
            : MultiArray<T>(dims)
        {
            switch (this->m_dims.size())
            {
                case 1:
                {
                    typename boost::multi_array<T,1>::extent_gen e;
                    m_a1.resize(e[this->m_dims[0]]);
                } break;

                case 2:
                {
                    typename boost::multi_array<T,2>::extent_gen e;
                    m_a2.resize(e[this->m_dims[0]][this->m_dims[1]]);
                } break;

                case 3:
                {
                    typename boost::multi_array<T,1>::extent_gen e;
                    m_a3.resize(e[this->m_dims[0]][this->m_dims[1]][this->m_dims[2]]);
                } break;

                case 4:
                {
                    typename boost::multi_array<T,1>::extent_gen e;
                    m_a4.resize(e[this->m_dims[0]][this->m_dims[1]][this->m_dims[2]][this->m_dims[3]]);
                } break;

                case 5:
                {
                    typename boost::multi_array<T,1>::extent_gen e;
                    m_a5.resize(e[this->m_dims[0]][this->m_dims[1]][this->m_dims[2]][this->m_dims[3]][this->m_dims[4]]);
                } break;

                default: throw std::out_of_range("only dimensions 1 .. 5 are supported");
            }
        };

        MultiArrayBoost(const vector<size_t> &dims, T default_value)
            : MultiArray<T>(dims,default_value)
        {
            switch (this->m_dims.size())
            {
                case 1:
                {
                    typename boost::multi_array<T,1>::extent_gen extend;
                    m_a1.resize(BOOST_A1(extend,this->m_dims));
                    BOOST_FILL_A1(this->m_dims,m_a1,default_value);
                } break;

                case 2:
                {
                    typename boost::multi_array<T,2>::extent_gen extend;
                    m_a2.resize(BOOST_A2(extend,this->m_dims));
                    BOOST_FILL_A2(this->m_dims,m_a2,default_value);
                } break;

                case 3:
                {
                    typename boost::multi_array<T,1>::extent_gen extend;
                    m_a3.resize(BOOST_A3(extend,this->m_dims));
                    BOOST_FILL_A3(this->m_dims,m_a3,default_value);
                } break;

                case 4:
                {
                    typename boost::multi_array<T,1>::extent_gen extend;
                    m_a4.resize(BOOST_A4(extend,this->m_dims));
                    BOOST_FILL_A4(this->m_dims,m_a4,default_value);
                } break;

                case 5:
                {
                    typename boost::multi_array<T,1>::extent_gen extend;
                    m_a5.resize(BOOST_A5(extend,this->m_dims));
                    BOOST_FILL_A5(this->m_dims,m_a5,default_value);
                } break;

                default: throw std::out_of_range("only dimensions 1 .. 5 are supported");
            }

        };

        MultiArrayBoost(const MultiArrayBlitz<T> &other)
        {
            switch (this->m_dims.size())
            {
                case 1: m_a1 = boost::multi_array<T,1>(other.m_a1); break;
                case 2: m_a2 = boost::multi_array<T,1>(other.m_a2); break;
                case 3: m_a3 = boost::multi_array<T,1>(other.m_a3); break;
                case 4: m_a4 = boost::multi_array<T,1>(other.m_a4); break;
                case 5: m_a5 = boost::multi_array<T,1>(other.m_a5); break;
                default: throw std::out_of_range("only 5 dimensions are currently supported");
            }
        }

        MultiArrayBoost<T>
        operator = (const MultiArrayBoost& other)
        {
            return MultiArrayBoost<T>(other);
        }

        /** Destructor
         */
        ~MultiArrayBoost() {};

#pragma mark -
#pragma mark Accessors

        T get(const vector<int> &index) const
        {
            switch (this->m_dims.size())
            {
                case 1: return BOOST_ID1(m_a1,index); break;
                case 2: return BOOST_ID2(m_a2,index); break;
                case 3: return BOOST_ID3(m_a3,index); break;
                case 4: return BOOST_ID4(m_a4,index); break;
                case 5: return BOOST_ID5(m_a5,index); break;
                default: throw std::out_of_range("only 5 dimensions are currently supported");
            }
        }

        void
        set(const vector<int> &index, const T &value)
        {
            switch (this->m_dims.size())
            {
                case 1: BOOST_ID1(m_a1,index) = value; break;
                case 2: BOOST_ID2(m_a2,index) = value; break;
                case 3: BOOST_ID3(m_a3,index) = value; break;
                case 4: BOOST_ID4(m_a4,index) = value; break;
                case 5: BOOST_ID5(m_a5,index) = value; break;
                default: throw std::out_of_range("only 5 dimensions are currently supported");
            }
        }

#pragma mark -
#pragma mark Stuff

        void resize(vector<size_t> dimensions)
        {
            this->m_dims = dimensions;

            switch (dimensions.size())
            {
                case 1:
                {
                    typename boost::multi_array<T,1>::extent_gen extend;
                    m_a1.resize(BOOST_A1(extend,this->m_dims));
                } break;

                case 2:
                {
                    typename boost::multi_array<T,2>::extent_gen extend;
                    m_a2.resize(BOOST_A2(extend,this->m_dims));
                } break;

                case 3:
                {
                    typename boost::multi_array<T,1>::extent_gen extend;
                    m_a3.resize(BOOST_A3(extend,this->m_dims));
                } break;

                case 4:
                {
                    typename boost::multi_array<T,1>::extent_gen extend;
                    m_a4.resize(BOOST_A4(extend,this->m_dims));
                } break;

                case 5:
                {
                    typename boost::multi_array<T,1>::extent_gen extend;
                    m_a5.resize(BOOST_A5(extend,this->m_dims));
                } break;

                default: throw std::out_of_range("only dimensions 1 .. 5 are supported");
            }
        }

        void populate_array(const T& value)
        {
            switch (this->m_dims.size())
            {
                case 1: BOOST_FILL_A1(this->m_dims, m_a1, value); break;
                case 2: BOOST_FILL_A2(this->m_dims, m_a2, value); break;
                case 3: BOOST_FILL_A3(this->m_dims, m_a3, value); break;
                case 4: BOOST_FILL_A4(this->m_dims, m_a4, value); break;
                case 5: BOOST_FILL_A5(this->m_dims, m_a5, value); break;
                default: throw std::out_of_range("only 5 dimensions are currently supported");
            }
        }

        void copy_from(const MultiArray<T> *other)
        {
            assert(this->m_dims == other->get_dimensions());

            switch (this->m_dims.size())
            {
                case 1: BOOST_COPY_A1(this->m_dims, other, default_value); break;
                case 2: BOOST_COPY_A2(this->m_dims, other, default_value); break;
                case 3: BOOST_COPY_A3(this->m_dims, other, default_value); break;
                case 4: BOOST_COPY_A4(this->m_dims, other, default_value); break;
                case 5: BOOST_COPY_A5(this->m_dims, other, default_value); break;
                default: throw std::out_of_range("only 5 dimensions are currently supported");
            }
        }

        size_t count_value(const T &value)
        {
            throw "not implemented";
        }

    };
}

#endif
