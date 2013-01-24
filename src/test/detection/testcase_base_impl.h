#ifndef _M3D_TEST_FSBASE_IMPL_H_
#define _M3D_TEST_FSBASE_IMPL_H_

#include <iostream>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <cf-algorithms/utils.h>

#pragma mark -
#pragma mark Test case 1 - ellipsoid sample means - data generation

using std::vector;
using std::cout;
using std::endl;
using namespace cfa::utils::coords;

template<class T>
FSTestBase<T>::FSTestBase() : m_settings(NULL), m_coordinate_system(NULL), m_totalPointCount(0), m_featureSpace(NULL), m_featureSpaceIndex(NULL) {}

template<class T>
FSTestBase<T>::~FSTestBase()
{
    if (m_featureSpaceIndex) {
        delete m_featureSpaceIndex;
        m_featureSpaceIndex = NULL;
    }

    if (m_settings) {
        delete m_settings;
        m_settings = NULL;
    }
    
    if (m_coordinate_system) {
        delete m_coordinate_system;
        m_coordinate_system = NULL;
    }
    
    if (m_featureSpace) {
        delete m_featureSpace;
        m_featureSpace = NULL;
    }
}

template<class T>
const char * FSTestBase<T>::dimensionName( size_t dimensionIndex )
{
    assert( dimensionIndex < 4 );
    
    const char* dimensionName = NULL;
    
    switch (dimensionIndex)
    {
        case 0:
            dimensionName = "x";
            break;
            
        case 1:
            dimensionName = "y";
            break;
            
        case 2:
            dimensionName = "z";
            break;
            
        case 3:
            dimensionName = "t";
            break;
            
        default:
            break;
    }
    
    return dimensionName;
}


template<class T>
void FSTestBase<T>::generate_dimensions()
{
    // Figure out the max bandwidth
    
    // create dimensions
    
    size_t nupoints = m_settings->num_gridpoints();
    
    vector<NcDim> dimensions;
    
    for ( size_t i = 0; i < m_settings->num_dimensions(); i++ )
    {
        // add grid dimension
        
        dimensions.push_back( m_file->addDim( dimensionName(i), nupoints + 1 ) );
    }
    
    // and dimension variables
    
    vector<NcVar> dimension_variables;
    
    for ( size_t i = 0; i < m_settings->num_dimensions(); i++ )
    {
        // create dimension variable. Values range [-max_h ... to max_h]
        
        float max_h = m_settings->axis_bound_values()[i];
        
        // Create variable with same name as dimension
        
        NcVar var = m_file->addVar( dimensionName(i), ncDouble, dimensions[i] );
        
        var.putAtt("valid_min", ncDouble, -max_h );
        
        var.putAtt("valid_max", ncDouble, max_h );
        
        // Write out it's axis data
        
        float* values = (float *) malloc( ( nupoints + 1 ) * sizeof(float) );
        
        for ( size_t gridIndex = 0; gridIndex < nupoints + 1; gridIndex++ )
        {
            values[gridIndex] = - max_h + gridIndex * ( 2 * max_h / ( (float)nupoints) );
        }
        
        var.putVar( values );
        
        free( values );
        
        // Save the dimension variable
        
        dimension_variables.push_back( var );
    }
    
    m_coordinate_system = new cfa::utils::coords::CoordinateSystem<T>( dimensions, dimension_variables );
}

template<class T>
void FSTestBase<T>::SetUp()
{
    // Select the correct point factory
    cfa::meanshift::PointFactory<T>::set_instance( new m3D::M3DPointFactory<T>() );
    
    // Create NetCDF file
    
    this->m_filename = m_settings->test_filename();
    
    try
    {
        this->m_file = new NcFile( this->m_filename, NcFile::replace, NcFile::nc4 );
    }
    catch ( const netCDF::exceptions::NcException& e )
    {
        std::cerr << "ERROR:could not open " << this->m_filename << " for writing : " << e.what() << std::endl;
        
        exit(-1);
    }
}

template<class T>
void FSTestBase<T>::TearDown()
{
    if ( m_file ) {
        delete m_file;
        m_file = NULL;
    }
}

template<class T>
void FSTestBase<T>::generate_featurespace()
{
    cout << "Creating featurespace ... ";
    
    this->m_featureSpace = new FeatureSpace<T>( this->m_filename, this->coordinate_system(), this->m_variables );
    
    this->m_featureSpaceIndex = PointIndex<T>::create( this->m_featureSpace );
    
    cout << "done. (" << m_featureSpace->size() << " points)" << endl;
}

template<class T>
string FSTestBase<T>::filename_from_current_testcase()
{
    const ::testing::TestInfo* const test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    
    string filename = test_info->test_case_name() + string(".nc");
    
    boost::replace_all( filename, "/", "_" );
    
    return filename;
    
}

#endif
