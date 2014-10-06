#ifndef _M3D_TESTS_FS_TESTCASES_H
#define _M3D_TESTS_FS_TESTCASES_H

#include <meanie3D/meanie3D.h>

#include <gtest/gtest.h>

#include <string>
#include <iostream>
#include <fstream>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#pragma mark - 
#pragma mark Switch individual tests on/off here

#define TEST_PRINT_PRECISION 3

#define RUN_2D 1
#define RUN_3D 1

#define RUN_CIRCULAR_PATTERN 1
#define RUN_UNWEIGHED_SAMPLE 1
#define RUN_WEIGHED_SAMPLE 1
#define RUN_ITERATION 1
#define RUN_CLUSTERING 1

#pragma mark -
#pragma mark Data Types 

using ::testing::Types;
using namespace std;

typedef Types< float, double > DataTypes;

#pragma mark -
#pragma mark Clustering 

#if RUN_CLUSTERING
#include "clustering.h"
#endif

#endif
