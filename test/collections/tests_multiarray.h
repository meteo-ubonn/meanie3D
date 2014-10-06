#ifndef cf_algorithms_blitz_idx_testcases_h
#define cf_algorithms_blitz_idx_testcases_h

#include <meanie3D/array/multiarray_blitz.h>

#include <gtest/gtest.h>

using namespace m3D;
using namespace testing;

#define TEST_A1(dim,array,idx,value) \
    idx.resize(1); \
    for(int i1=0;i1<dim[0];i1++) { idx[0]=i1; \
        EXPECT_EQ(value, array.get(idx)); }

#define TEST_A2(dim,array,idx,value) \
    idx.resize(2); \
    for(int i1=0;i1<dim[0];i1++) { idx[0]=i1; \
        for(int i2=0;i2<dim[1];i2++) { idx[1]=i2; \
            EXPECT_EQ(value, array.get(idx)); }}

#define TEST_A3(dim,array,idx,value) \
    idx.resize(3); \
    for(int i1=0;i1<dim[0];i1++) { idx[0]=i1; \
        for(int i2=0;i2<dim[1];i2++) { idx[1]=i2; \
            for(int i3=0;i3<dim[2];i3++) { idx[2]=i3; \
                EXPECT_EQ(value, array.get(idx)); }}}

#define TEST_A4(dim,array,idx,value) \
    idx.resize(4); \
    for(int i1=0;i1<dim[0];i1++) { idx[0]=i1; \
        for(int i2=0;i2<dim[1];i2++) { idx[1]=i2; \
            for(int i3=0;i3<dim[2];i3++) { idx[2]=i3; \
                for(int i4=0;i4<dim[3];i4++) { idx[3]=i4; \
                    EXPECT_EQ(value, array.get(idx)); }}}}

#define TEST_A5(dim,array,idx,value) \
    idx.resize(5); \
    for(int i1=0;i1<dim[0];i1++) { idx[0]=i1; \
        for(int i2=0;i2<dim[1];i2++) { idx[1]=i2; \
            for(int i3=0;i3<dim[2];i3++) { idx[2]=i3; \
                for(int i4=0;i4<dim[3];i4++) { idx[3]=i4; \
                    for(int i5=0;i5<dim[4];i5++) { idx[4]=i5; \
                        EXPECT_EQ(value, array.get(idx)); }}}}}

#define SET_A1(dim,array,idx,value) \
    idx.resize(1); \
    for(int i1=0;i1<dim[0];i1++) { idx[0]=i1; \
        array.set(idx,value); }

#define SET_A2(dim,array,idx,value) \
    idx.resize(2); \
    for(int i1=0;i1<dim[0];i1++) { idx[0]=i1; \
        for(int i2=0;i2<dim[1];i2++) { idx[1]=i2; \
            array.set(idx,value); }}


#define SET_A3(dim,array,idx,value) \
    idx.resize(3); \
    for(int i1=0;i1<dim[0];i1++) { idx[0]=i1; \
        for(int i2=0;i2<dim[1];i2++) { idx[1]=i2; \
            for(int i3=0;i3<dim[2];i3++) { idx[2]=i3; \
                array.set(idx,value); }}}

#define SET_A4(dim,array,idx,value) \
    idx.resize(4); \
    for(int i1=0;i1<dim[0];i1++) { idx[0]=i1; \
        for(int i2=0;i2<dim[1];i2++) { idx[1]=i2; \
            for(int i3=0;i3<dim[2];i3++) { idx[2]=i3; \
                for(int i4=0;i4<dim[3];i4++) { idx[3]=i4; \
                    array.set(idx,value); }}}}

#define SET_A5(dim,array,idx,value) \
    idx.resize(5); \
    for(int i1=0;i1<dim[0];i1++) { idx[0]=i1; \
        for(int i2=0;i2<dim[1];i2++) { idx[1]=i2; \
            for(int i3=0;i3<dim[2];i3++) { idx[2]=i3; \
                for(int i4=0;i4<dim[3];i4++) { idx[3]=i4; \
                    for(int i5=0;i5<dim[4];i5++) { idx[4]=i5; \
                        array.set(idx,value); }}}}}


#pragma mark -
#pragma mark Blitz Multi-Array

template <typename T>
class MultiArrayBlitzTest : public testing::Test {};

typedef testing::Types< float, double > VectorDataTypes;

TYPED_TEST_CASE( MultiArrayBlitzTest, VectorDataTypes );

TYPED_TEST( MultiArrayBlitzTest, VectorDataTypes )
{
    vector<size_t>  dims;
    vector<int>     index;

    // 1D
//    dims.resize(1);
//    dims[0]=5;
//    
//    // Test default constructor
//    
//    MultiArrayBlitz<TypeParam> a11(dims);
//
//    // Test set / get
//
//    SET_A1(dims,a11,index,120);
//    TEST_A1(dims,a11,index,120);
//    
//    // test populate
//    a11.populate_array(125);
//    TEST_A1(dims,a11,index,125);
//
//    // Test constructor with default value
//    
//    MultiArrayBlitz<TypeParam> a12(dims,100);
//    TEST_A1(dims,a12,index,100);

    // 2D
    
    dims.resize(2);
    dims[0] = 1000;
    dims[1] = 1000;

    MultiArrayBlitz<TypeParam> a21(dims);
    
    // Test set / get
    
    SET_A2(dims,a21,index,200);
    TEST_A2(dims,a21,index,200);
    
    // Test constructor with default value
    
    MultiArrayBlitz<TypeParam> a22(dims,250);
    TEST_A2(dims,a22,index,250);

    // 3D
    
    dims.resize(3);
    dims[0] = 1000;
    dims[1] = 1000;
    dims[2] = 100;
    
    MultiArrayBlitz<TypeParam> a31(dims);
    
    // Test set / get
    
    SET_A3(dims,a31,index,300);
    TEST_A3(dims,a31,index,300);
    
    // Test constructor with default value
    
    MultiArrayBlitz<TypeParam> a32(dims,350);
    TEST_A3(dims,a32,index,350);

    // 4D
    
    dims.resize(4);
    dims[0] = 10;
    dims[1] = 10;
    dims[2] = 10;
    dims[3] = 10;
    
    MultiArrayBlitz<TypeParam> a41(dims);
    
    // Test set / get
    
    SET_A4(dims,a41,index,400);
    TEST_A4(dims,a41,index,400);
    
    // Test constructor with default value
    
    MultiArrayBlitz<TypeParam> a42(dims,450);
    TEST_A4(dims,a42,index,450);

    // 4D
    
    dims.resize(5);
    dims[0] = 10;
    dims[1] = 10;
    dims[2] = 10;
    dims[3] = 10;
    dims[4] = 10;
    
    MultiArrayBlitz<TypeParam> a51(dims);
    
    // Test set / get
    
    SET_A5(dims,a51,index,500);
    TEST_A5(dims,a51,index,500);
    
    // Test constructor with default value
    
    MultiArrayBlitz<TypeParam> a52(dims,550);
    TEST_A5(dims,a52,index,550);
}

#pragma mark -
#pragma mark Recursive Multi-Array

template <typename T>
class MultiArrayRecursiveTest : public testing::Test {};

TYPED_TEST_CASE( MultiArrayRecursiveTest, VectorDataTypes );

TYPED_TEST( MultiArrayRecursiveTest, VectorDataTypes )
{
    vector<size_t>  dims;
    vector<int>     index;
    
    // 1D
//    dims.resize(1);
//    dims[0]=5;
//    
//    // Test default constructor
//    
//    MultiArrayRecursive<TypeParam> a11(dims);
//    
//    // Test set / get
//    
//    SET_A1(dims,a11,index,120);
//    TEST_A1(dims,a11,index,120);
//    
//    // Test constructor with default value
//    
//    MultiArrayRecursive<TypeParam> a12(dims,100);
//    TEST_A1(dims,a12,index,100);
    
    // 2D
    
    dims.resize(2);
    dims[0] = 1000;
    dims[1] = 1000;
    
    MultiArrayRecursive<TypeParam> a21(dims);
    
    // Test set / get
    
    SET_A2(dims,a21,index,200);
    TEST_A2(dims,a21,index,200);
    
    // Test constructor with default value
    
    MultiArrayRecursive<TypeParam> a22(dims,250);
    TEST_A2(dims,a22,index,250);
    
    // 3D
    
    dims.resize(3);
    dims[0] = 1000;
    dims[1] = 1000;
    dims[2] = 100;
    
    MultiArrayRecursive<TypeParam> a31(dims);
    
    // Test set / get
    
    SET_A3(dims,a31,index,300);
    TEST_A3(dims,a31,index,300);
    
    // Test constructor with default value
    
    MultiArrayRecursive<TypeParam> a32(dims,350);
    TEST_A3(dims,a32,index,350);
    
    // 4D
    
    dims.resize(4);
    dims[0] = 10;
    dims[1] = 10;
    dims[2] = 10;
    dims[3] = 10;
    
    MultiArrayRecursive<TypeParam> a41(dims);
    
    // Test set / get
    
    SET_A4(dims,a41,index,400);
    TEST_A4(dims,a41,index,400);
    
    // Test constructor with default value
    
    MultiArrayRecursive<TypeParam> a42(dims,450);
    TEST_A4(dims,a42,index,450);
    
    // 4D
    
    dims.resize(5);
    dims[0] = 10;
    dims[1] = 10;
    dims[2] = 10;
    dims[3] = 10;
    dims[4] = 10;
    
    MultiArrayRecursive<TypeParam> a51(dims);
    
    // Test set / get
    
    SET_A5(dims,a51,index,500);
    TEST_A5(dims,a51,index,500);
    
    // Test constructor with default value
    
    MultiArrayRecursive<TypeParam> a52(dims,550);
    TEST_A5(dims,a52,index,550);
};

#pragma mark -
#pragma mark Blitz Multi-Array

template <typename T>
class MultiArrayBoostTest : public testing::Test {};

typedef testing::Types< float, double > VectorDataTypes;

TYPED_TEST_CASE( MultiArrayBoostTest, VectorDataTypes );

TYPED_TEST( MultiArrayBoostTest, VectorDataTypes )
{
    vector<size_t>  dims;
    vector<int>     index;
    
    // 1D
    //    dims.resize(1);
    //    dims[0]=5;
    //
    //    // Test default constructor
    //
    //    MultiArrayBoost<TypeParam> a11(dims);
    //
    //    // Test set / get
    //
    //    BLITZ_SET_A1(dims,a11,index,120);
    //    BLITZ_TEST_A1(dims,a11,index,120);
    //
    //    // test populate
    //    a11.populate_array(125);
    //    BLITZ_TEST_A1(dims,a11,index,125);
    //
    //    // Test constructor with default value
    //
    //    MultiArrayBoost<TypeParam> a12(dims,100);
    //    BLITZ_TEST_A1(dims,a12,index,100);
    
    // 2D
    
    dims.resize(2);
    dims[0] = 1000;
    dims[1] = 1000;
    
    MultiArrayBoost<TypeParam> a21(dims);
    
    // Test set / get
    
    SET_A2(dims,a21,index,200);
    TEST_A2(dims,a21,index,200);
    
    // Test constructor with default value
    
    MultiArrayBoost<TypeParam> a22(dims,250);
    TEST_A2(dims,a22,index,250);
    
    // 3D
    
    dims.resize(3);
    dims[0] = 1000;
    dims[1] = 1000;
    dims[2] = 100;
    
    MultiArrayBoost<TypeParam> a31(dims);
    
    // Test set / get
    
    SET_A3(dims,a31,index,300);
    TEST_A3(dims,a31,index,300);
    
    // Test constructor with default value
    
    MultiArrayBoost<TypeParam> a32(dims,350);
    TEST_A3(dims,a32,index,350);
    
    // 4D
    
    dims.resize(4);
    dims[0] = 10;
    dims[1] = 10;
    dims[2] = 10;
    dims[3] = 10;
    
    MultiArrayBoost<TypeParam> a41(dims);
    
    // Test set / get
    
    SET_A4(dims,a41,index,400);
    TEST_A4(dims,a41,index,400);
    
    // Test constructor with default value
    
    MultiArrayBoost<TypeParam> a42(dims,450);
    TEST_A4(dims,a42,index,450);
    
    // 4D
    
    dims.resize(5);
    dims[0] = 10;
    dims[1] = 10;
    dims[2] = 10;
    dims[3] = 10;
    dims[4] = 10;
    
    MultiArrayBoost<TypeParam> a51(dims);
    
    // Test set / get
    
    SET_A5(dims,a51,index,500);
    TEST_A5(dims,a51,index,500);
    
    // Test constructor with default value
    
    MultiArrayBoost<TypeParam> a52(dims,550);
    TEST_A5(dims,a52,index,550);
};

#endif

