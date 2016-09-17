#ifndef M3D_TESTS_FS_TESTCASES_H
#define M3D_TESTS_FS_TESTCASES_H

#include <meanie3D/meanie3D.h>

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

#pragma mark -
#pragma mark Data Types 

using namespace testing;
using namespace std;

typedef Types<float, double> DataTypes;

#pragma mark -
#pragma mark Circular Pattern Tests

#if RUN_CIRCULAR_PATTERN

#include "circular_pattern.h"

#endif

#pragma mark -
#pragma mark Unweighed uniform distribution sample (no weight function)

#if RUN_UNWEIGHED_SAMPLE

#include "uniform.h"

#endif

#pragma mark -
#pragma mark Weighed uniform distribution sample

#if RUN_WEIGHED_SAMPLE

#include "weighed.h"

#endif

#pragma mark -
#pragma mark Mean-Shift Iterations

#if RUN_ITERATION

#include "iteration.h"

#endif

#endif
