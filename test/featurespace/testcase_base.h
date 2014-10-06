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

using namespace std;
using namespace testing;
using namespace netCDF;
using namespace m3D;

template <class T> 
class FSTestBase : public Test
{
    
protected:
    
    //
    // Protected member variables
    //
    
    // Common settings, like grid resolution etc.
    FSTestSettings  *m_settings;
    
    // NetCDF related
    
    std::string             m_filename;
    
    NcFile                  *m_file;
    
    CoordinateSystem<T>     *m_coordinate_system;
    
//    vector<NcDim *> m_dimensions;
//
//    vector<NcVar *> m_dimension_variables;
//    
//    // Stores dimension variable data
//    map< NcVar* , float* > m_dimension_data;

    vector<NcVar>           m_variables;

    size_t                  m_pointCount;            
    
    size_t                  m_totalPointCount;       
    
    FeatureSpace<T>         *m_featureSpace;
    
    PointIndex<T>    *m_featureSpaceIndex;
    
    //
    // Protected methods
    //
    
    /** Returns variables names for dimension variables to be created as part of the test
     * procedure. 0=x, 1=y, 2=z, 3=t.
     */
    const char *dimensionName( size_t dimensionIndex );

    /** Generates a number of dimensions in the given NetCDF file, using the values
     * in the current settings object. Also creates the dimension variables and buffers
     * their values for quick access in m_dimension_data
     */
    void generate_dimensions();
    
    /** Auto-generates a filename for the netcdf from the currently running 
     * testcase info
     */
    string filename_from_current_testcase();
    
public:

    FSTestBase();
    
    ~FSTestBase();
    
    virtual void SetUp();
    
    virtual void TearDown();
    
    NcFile * file() { return m_file; };
    
    vector<NcVar *> variables() { return m_variables; };
    
    CoordinateSystem<T> *coordinate_system() { return m_coordinate_system; };
    
    /** Helper method. Adds a variable with the given name to the file. The variable
     * has the given name. The valid_min and valid_max attributes can also be set (default 0..1)
     */
    NcVar add_variable( string name, T valid_min = 0.0, T valid_max = 1.0 )
    {
        // test case 1 : ellipsis for unweighed sample mean
        
        std::string type_name = typeid(T).name();
        
        if ( ! ( type_name == "f" || type_name == "d") )
        {
            cerr << "Only float and double are supported at this time." << endl;
            exit(-1);
        }
        
        NcType type = ( type_name == "f" ) ? ncFloat.getTypeClass() : ncDouble.getTypeClass();
        
        NcVar var = file()->addVar( name, type, this->coordinate_system()->dimensions() );
        
        var.putAtt( "valid_min", type, valid_min );
        
        var.putAtt( "valid_max", type, valid_max );
        
        this->m_variables.push_back( var );
        
        return var;
    }
    
    void generate_featurespace();
};

#include "testcase_base_impl.h"

#endif
