#ifndef M3D_MAP_TEST_H
#define M3D_MAP_TEST_H

//
//  kernel_testcases_h
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#include <meanie3D/utils/map_utils.h>

#include <gtest/gtest.h>
#include <map>

using namespace m3D;
using namespace m3D::utils;
using namespace testing;

class IDMapTest : public testing::Test
{
public:
};

TEST(IDMapTest, IDMapTest) {
    m3D::id_map_t m;

    m3D::id_set_t s1;
    s1.insert(1);
    s1.insert(2);

    m3D::id_set_t s2;
    s2.insert(3);
    s2.insert(4);

    m[(m3D::id_t) 0] = s1;
    m[(m3D::id_t) 1] = s2;

    const char *expected_str = "[0:{1,2};1:{3,4};]";

    std::string serialized = maps::id_map_to_string(m);
    EXPECT_STREQ(expected_str, serialized.c_str());

    m3D::id_map_t deserialized = maps::id_map_from_string(serialized);
    EXPECT_EQ(m, deserialized);
}

#endif
