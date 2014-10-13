#ifndef M3D_TEST_FSBASE_H
#define M3D_TEST_FSBASE_H

//
//  circular_pattern.h
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#include <meanie3D/meanie3D.h>

#include <gtest/gtest.h>
#include <boost/lexical_cast.hpp>

#include <netcdf>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <typeinfo>

#include "settings.h"

#pragma mark -
#pragma mark Declarations

#define INFO_ENABLED false
#define INFO if (INFO_ENABLED) cout << "INFO:" 

using namespace std;
using namespace testing;
using namespace netCDF;
using namespace m3D;

template <class T>
class FSTestBase : public Test 
{
protected:

    NcFile              *m_file;
    std::string         m_filename;
    vector<string>      m_variable_names;

    DataStore<T>        *m_data_store;
    CoordinateSystem<T> *m_coordinate_system;
    FeatureSpace<T>     *m_featureSpace;
    PointIndex<T>       *m_featureSpaceIndex;
    
    FSTestSettings      *m_settings;

    size_t              m_pointCount;
    size_t              m_totalPointCount;

    /** Returns variables names for dimension variables to be 
     * created as part of the test
     * procedure. 0=x, 1=y, 2=z, 3=t.
     * 
     * @param dimensionIndex
     * @return 
     */
    const char *dimensionName(size_t dimensionIndex);

    /** Generates a number of dimensions in the given NetCDF 
     * file, using the values in the current settings object. 
     * Also creates the dimension variables and buffers
     * their values for quick access in m_dimension_data
     */
    void generate_dimensions();

    /** Auto-generates a filename for the netcdf from the 
     * currently running testcase info
     * 
     * @return 
     */
    string filename_from_current_testcase();

public:
    
    static const T FILL_VALUE;

    /**
     * 
     */
    FSTestBase();

    /**
     * 
     */
    ~FSTestBase();

    /**
     * 
     */
    virtual void SetUp();

    /**
     * 
     */
    virtual void TearDown();

    /**
     * 
     * @return 
     */
    NcFile *file();
    
    /**
     */
    void reopen_file_for_reading();

    /**
     * 
     * @return 
     */
    CoordinateSystem<T> *coordinate_system();
    
    /** Helper method. Adds a variable with the given name to the file. 
     * The variable has the given name. The valid_min and valid_max 
     * attributes can also be set (default 0..1)
     * 
     * @param name
     * @param valid_min
     * @param valid_max
     * @return the created variable
     */
    NcVar add_variable(string name, T valid_min = 0.0, T valid_max = 1.0);
    
    /**
     * 
     */
    void generate_featurespace();
};

#include "testcase_base_impl.h"

#endif
