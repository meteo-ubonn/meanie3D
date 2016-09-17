#ifndef M3D_SET_TEST_H
#define M3D_SET_TEST_H

//
//  kernel_testcases_h
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#include <meanie3D/utils/set_utils.h>

#include <gtest/gtest.h>
#include <set>

using namespace m3D::utils;
using namespace testing;

template<typename T>
class SetTest : public testing::Test
{
public:
};

typedef testing::Types<float, double, m3D::id_t> SetDataTypes;

TYPED_TEST_CASE(SetTest, SetDataTypes);

TYPED_TEST(SetTest, SetDataTypes) {
    // to_string / from_string

    set <TypeParam> s;
    s.insert(1);
    s.insert(2);
    s.insert(3);
    s.insert(4);

    const char *expected_serialized = "{1,2,3,4}";

    std::string serialized = sets::to_string<TypeParam>(s);

    EXPECT_STREQ(expected_serialized, serialized.c_str());

    EXPECT_EQ(s, sets::from_string<TypeParam>(serialized));
}

TEST(SetTest, Float) {
    typedef float TypeParam;
}

TEST(SetTest, Double) {
    typedef double TypeParam;
}

TEST(SetTest, ID_T) {
    typedef m3D::id_t TypeParam;
}

#endif
