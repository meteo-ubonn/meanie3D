#ifndef M3D_VECTOR_TEST_H
#define M3D_VECTOR_TEST_H

//
//  kernel_testcases_h
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#include <meanie3D/utils/vector_utils.h>

#include <gtest/gtest.h>
#include <limits>
#include <vector>

using namespace std;
using namespace testing;
using namespace m3D::utils::vectors;

template <typename T>
class VectorTest : public testing::Test
{
public:
};

typedef testing::Types< float, double > VectorDataTypes;

TYPED_TEST_CASE(VectorTest, VectorDataTypes);

TYPED_TEST(VectorTest, VectorDataTypes)
{
    // T vector_norm( vector<T> *v )
    vector<TypeParam> v1, v2, v3;
    TypeParam s1;
    TypeParam sqr5 = sqrt(5);

    // vector_norm

    v1 = vector<TypeParam>(5, 1);
    s1 = vector_norm(&v1);
    EXPECT_EQ(sqr5, s1);

    // T vector_norm( vector<T> &v )

    v1 = vector<TypeParam>(5, 1);
    s1 = vector_norm(v1);
    EXPECT_EQ(sqr5, s1);

    v1 = vector<TypeParam>(5, 1);
    v2 = vector<TypeParam>(5, -1);
    v1 += v2;
    EXPECT_TRUE(vector_is_null(v1));

    // operator v + v
    v1 = vector<TypeParam>(5, 0.5);
    v2 = vector<TypeParam>(5, 0.5);
    v3 = v1 + v2;
    EXPECT_EQ(sqr5, vector_norm(v3));

    // operator v += v
    v1 = vector<TypeParam>(5, 0.5);
    v2 = vector<TypeParam>(5, 0.5);
    v1 += v2;
    EXPECT_EQ(sqr5, vector_norm(v1));

    // operator v - v
    v1 = vector<TypeParam>(5, 2);
    v2 = vector<TypeParam>(5, 1);
    v3 = v1 - v2;
    EXPECT_EQ(sqr5, vector_norm(v3));

    // operator v -= v
    v1 = vector<TypeParam>(5, 2);
    v2 = vector<TypeParam>(5, 1);
    v1 -= v2;
    EXPECT_EQ(sqr5, vector_norm(v1));

    // operator s * v
    s1 = 0.5;
    v1 = vector<TypeParam>(5, 2);
    v2 = s1 * v1;
    EXPECT_EQ(sqr5, vector_norm(v2));

    // operator v * s
    s1 = 0.5;
    v1 = vector<TypeParam>(5, 2);
    v2 = v1 * s1;
    EXPECT_EQ(sqr5, vector_norm(v2));

    // operator v *= s
    s1 = 0.5;
    v1 = vector<TypeParam>(5, 2);
    v1 *= s1;
    EXPECT_EQ(sqr5, vector_norm(v1));

    // operator v /= s
    s1 = 2.0;
    v1 = vector<TypeParam>(5, 2);
    v1 /= s1;
    EXPECT_EQ(sqr5, vector_norm(v1));

    // operaror v / s
    s1 = 2.0;
    v1 = vector<TypeParam>(5, 2);
    v2 = v1 / s1;
    EXPECT_EQ(sqr5, vector_norm(v2));

    // operator v * v (scalar product)
    v1 = vector<TypeParam>(5, 2);
    v2 = vector<TypeParam>(5, 0.5);
    s1 = v1 * v2;
    EXPECT_EQ(5, s1);

    // vector_is_null

    v1 = vector<TypeParam>(5, 1);
    EXPECT_FALSE(vector_is_null(&v1));

    v1 = vector<TypeParam>(5, 0);
    EXPECT_TRUE(vector_is_null(&v1));

    v1 = vector<TypeParam>(5, 1);
    EXPECT_FALSE(vector_is_null(v1));

    v1 = vector<TypeParam>(5, 0);
    EXPECT_TRUE(vector_is_null(v1));

    // within_range

    v1 = vector<TypeParam>(5, 2);
    v2 = vector<TypeParam>(5, 1);
    v3 = vector<TypeParam>(5, 1);
    EXPECT_TRUE(within_range(v1, v2, v3));

    v1 = vector<TypeParam>(5, 3);
    v2 = vector<TypeParam>(5, 1);
    v3 = vector<TypeParam>(5, 1);
    EXPECT_FALSE(within_range(v1, v2, v3));

    // closer_than

    v1 = vector<TypeParam>(5, 1.99999);
    v2 = vector<TypeParam>(5, 1);
    v3 = vector<TypeParam>(5, 1);
    EXPECT_TRUE(closer_than(v1, v2, v3));

    v1 = vector<TypeParam>(5, 2);
    v2 = vector<TypeParam>(5, 1);
    v3 = vector<TypeParam>(5, 1);
    EXPECT_FALSE(closer_than(v1, v2, v3));

    // vector_angle

    v1 = vector<TypeParam>(5, 1);
    v2 = vector<TypeParam>(5, 1);
    s1 = vector_angle(v1, v2);
    EXPECT_NEAR(0.0, s1, 10e-8);

    v1 = vector<TypeParam>(2);
    v1[0] = 1.0;
    v1[1] = 0.0;
    v2 = vector<TypeParam>(2);
    v2[0] = 0.0;
    v2[1] = 1.0;
    s1 = vector_angle(v1, v2);
    EXPECT_EQ((TypeParam) M_PI_2, s1);

    // mahalabonis_distance_sqr
    v1 = vector<TypeParam>(5, 2);
    v2 = vector<TypeParam>(5, 1);
    v3 = vector<TypeParam>(5, 1);
    s1 = mahalabonis_distance_sqr(v1, v2, v3);
    EXPECT_EQ(5, s1);

    // vector_diagonal_product

    v1 = vector<TypeParam>(5, 2);
    v2 = vector<TypeParam>(5, 0.5);
    v3 = vector_diagonal_product(v1, v2);
    EXPECT_EQ(sqr5, vector_norm(v3));

    v1 = vector<TypeParam>(3);
    v1[0] = 1;
    v1[1] = 2;
    v1[2] = 2;
    EXPECT_EQ(index_of_first<TypeParam>(v1, 1), 0);
    EXPECT_EQ(index_of_first<TypeParam>(v1, 2), 1);
    EXPECT_EQ(index_of_first<TypeParam>(v1, 5), -1);

    // to_string / from_string

    v1 = vector<TypeParam>(3);
    v1[0] = 1;
    v1[1] = 2;
    v1[2] = 3;

    EXPECT_EQ("(1,2,3)", to_string<TypeParam>(v1));
    EXPECT_EQ(v1, from_string<TypeParam>("(1,2,3)"));
}

TEST(VectorTest, Double)
{
    typedef double TypeParam;
}

#endif
