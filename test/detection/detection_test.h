#ifndef cf_algorithms_featurespace_test_h
#define cf_algorithms_featurespace_test_h

//
//  featurespace_test.h
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#include <gtest/gtest.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <set>
#include <map>
#include <vector>
#include <netcdf>

#include <boost/program_options.hpp>

#pragma mark Constants

// Set some constants

/** Number of dimensions and dimension variables in the testfile (max 4) 
 */
const size_t num_dimensions = 3;

/** Number of non-dimension variables in the testfile 
 */
const size_t num_vars = 1;

/** Dimensionality of feature space 
 */
const size_t FS_DIM = num_dimensions + num_vars;

/** Name of the file used for the test data 
 */
const char *SAMPLETEST_FILENAME = "sampletest.nc";

/** Number of grid points in each direction 
 */
size_t NUM_GRIDPOINTS = 100;

/** 'Sharpness' of the test data ellipsoids 
 */
float ELLIPSE_FUZZINESS = 0.05;

/** 'Exactness' of the meanshift result vector 
 */
float MEANSHIFT_FUZZINESS = 1.5;

/** Data type for the feature space template instantiation 
 */
#define FEATURESPACE_TYPE float


class FeatureSpaceTest : public testing::Test
{
public:
    
    // You can remove any or all of the following functions if its body
    // is empty.
    
    FeatureSpaceTest();

    virtual ~FeatureSpaceTest();
        // You can do clean-up work that doesn't throw exceptions here.
    
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:
    
    virtual void SetUp();
        // Code here will be called immediately after the constructor (right
        // before each test).
    
    virtual void TearDown();
        // Code here will be called immediately after each test (right
        // before the destructor).
    

    size_t pointCount;            // number of points in the data set of the current run
    
    size_t totalPointCount;       // overall number of data points in the file
    
    netCDF::NcFile *file;

    // Stores dimension variable data
    typename std::map< netCDF::NcVar * , float * > dimensionVars;
    
};

#endif
