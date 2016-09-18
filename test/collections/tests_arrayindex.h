#ifndef M3D_ARRAY_INDEX_TEST_H
#define M3D_ARRAY_INDEX_TEST_H

#include <meanie3D/array/array_index.h>
#include <meanie3D/featurespace.h>

#include <gtest/gtest.h>
#include <vector>

using namespace m3D;
using namespace testing;

#pragma mark -
#pragma mark Blitz Multi-Array

typedef testing::Types<float, double> VectorDataTypes;

// ONE DEE

template<typename T>
class ArrayIndexTest1D : public testing::Test
{
};

TYPED_TEST_CASE(ArrayIndexTest1D, VectorDataTypes);

TYPED_TEST(ArrayIndexTest1D, VectorDataTypes) {
    PointFactory<TypeParam>::set_instance(new PointDefaultFactory<TypeParam>());

    typename Point<TypeParam>::list points;

    vector<size_t> dimensions(1);
    dimensions[0] = 10;

    // construct point list

    vector<int> g(dimensions.size(), 0);
    vector<TypeParam> c(dimensions.size(), 0);

    for (int iz = 0; iz < dimensions[0]; iz++) {
        g[0] = iz;
        c[0] = 0.1 * iz;

        typename Point<TypeParam>::ptr p = PointFactory<TypeParam>::get_instance()->create(g, c, c);
        points.push_back(p);
    }

    // Construct an non-copying ArrayIndex

    ArrayIndex<TypeParam> index(dimensions, points, false);

    // Retrieve all points and compare

    for (size_t pi = 0; pi < points.size(); pi++) {
        typename Point<TypeParam>::ptr a = points.at(pi);
        typename Point<TypeParam>::ptr b = index.get(a->gridpoint);

        EXPECT_EQ(a->gridpoint, b->gridpoint);
        EXPECT_EQ(a->coordinate, b->coordinate);
        EXPECT_EQ(a->values, b->values);
    }


    // clean up

    while (!points.empty()) {
        typename Point<TypeParam>::ptr a = points.back();
        points.pop_back();
        delete a;
    }
}

// TWO DEE

template<typename T>
class ArrayIndexTest2D : public testing::Test
{
};

TYPED_TEST_CASE(ArrayIndexTest2D, VectorDataTypes);

TYPED_TEST(ArrayIndexTest2D, VectorDataTypes) {
    PointFactory<TypeParam>::set_instance(new PointDefaultFactory<TypeParam>());

    typename Point<TypeParam>::list points;

    vector<size_t> dimensions(2);
    dimensions[0] = 10;
    dimensions[1] = 10;

    // construct point list

    vector<int> g(dimensions.size(), 0);
    vector<TypeParam> c(dimensions.size(), 0);

    for (int iz = 0; iz < dimensions[0]; iz++) {
        g[0] = iz;
        c[0] = 0.1 * iz;

        for (int iy = 0; iy < dimensions[1]; iy++) {
            g[1] = iy;
            c[1] = 0.1 * iy;

            typename Point<TypeParam>::ptr p = PointFactory<TypeParam>::get_instance()->create(g, c, c);
            points.push_back(p);
        }
    }

    // Construct an non-copying ArrayIndex

    ArrayIndex<TypeParam> index(dimensions, points, false);

    // Retrieve all points and compare

    for (size_t pi = 0; pi < points.size(); pi++) {
        typename Point<TypeParam>::ptr a = points.at(pi);
        typename Point<TypeParam>::ptr b = index.get(a->gridpoint);

        EXPECT_EQ(a->gridpoint, b->gridpoint);
        EXPECT_EQ(a->coordinate, b->coordinate);
        EXPECT_EQ(a->values, b->values);
    }


    // clean up

    while (!points.empty()) {
        typename Point<TypeParam>::ptr a = points.back();
        points.pop_back();
        delete a;
    }
}

// THREE DEE

template<typename T>
class ArrayIndexTest3D : public testing::Test
{
};

TYPED_TEST_CASE(ArrayIndexTest3D, VectorDataTypes);

TYPED_TEST(ArrayIndexTest3D, VectorDataTypes) {
    PointFactory<TypeParam>::set_instance(new PointDefaultFactory<TypeParam>());

    typename Point<TypeParam>::list points;

    vector<size_t> dimensions(3);
    dimensions[0] = 10;
    dimensions[1] = 10;
    dimensions[2] = 10;

    // construct point list

    vector<int> g(dimensions.size(), 0);
    vector<TypeParam> c(dimensions.size(), 0);

    for (int iz = 0; iz < dimensions[0]; iz++) {
        g[0] = iz;
        c[0] = 0.1 * iz;

        for (int iy = 0; iy < dimensions[1]; iy++) {
            g[1] = iy;
            c[1] = 0.1 * iy;

            for (int ix = 0; ix < dimensions[2]; ix++) {
                g[2] = ix;
                c[2] = 0.1 * ix;

                typename Point<TypeParam>::ptr p = PointFactory<TypeParam>::get_instance()->create(g, c, c);
                points.push_back(p);
            }
        }
    }

    // Construct an non-copying ArrayIndex

    ArrayIndex<TypeParam> index(dimensions, points, false);

    // Retrieve all points and compare

    for (size_t pi = 0; pi < points.size(); pi++) {
        typename Point<TypeParam>::ptr a = points.at(pi);
        typename Point<TypeParam>::ptr b = index.get(a->gridpoint);

        EXPECT_EQ(a->gridpoint, b->gridpoint);
        EXPECT_EQ(a->coordinate, b->coordinate);
        EXPECT_EQ(a->values, b->values);
    }


    // clean up

    while (!points.empty()) {
        typename Point<TypeParam>::ptr a = points.back();
        points.pop_back();
        delete a;
    }

}

#endif

