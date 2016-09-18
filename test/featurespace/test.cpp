//
//  featurespace_test.cpp
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

//#define GTEST_GCC_VER_ 40303
//#define GTEST_HAS_RTTI 1
//#define GTEST_USE_OWN_TR1_TUPLE 1
//#define GTEST_HAS_TR1_TUPLE 0

#include <gtest/gtest.h>
#include "testcases.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
