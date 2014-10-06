//
//  kernel_test.cpp
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#define GTEST_HAS_TR1_TUPLE 0

#include <gtest/gtest.h>
#include <meanie3D/meanie3D.h>

#include "tests_vector.h"
#include "test_linear_mapping.h"
#include "tests_map.h"
#include "tests_set.h"
#include "tests_arrayindex.h"
#include "tests_multiarray.h"

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    
    return RUN_ALL_TESTS();
}