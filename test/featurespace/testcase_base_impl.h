#ifndef M3D_TEST_FSBASE_IMPL_H
#define M3D_TEST_FSBASE_IMPL_H

#include <meanie3D/meanie3D.h>

#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>

using namespace std;

#pragma mark -
#pragma mark Test case 1 - ellipsoid sample means - data generation

template<class T> 
FSTestBase<T>::FSTestBase() : m_settings(NULL), m_coordinate_system(NULL),
m_totalPointCount(0), m_featureSpace(NULL), m_featureSpaceIndex(NULL) {}

template<class T> 
FSTestBase<T>::~FSTestBase()
{
    if ( m_settings != NULL )
    {
        delete m_settings;
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
    
    m_coordinate_system = new CoordinateSystem<T>( dimensions, dimension_variables );
}

template<class T> 
void FSTestBase<T>::SetUp()
{
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
    delete m_file;
    
    m_file = NULL;

    delete m_coordinate_system;
    
    m_coordinate_system = NULL;
    
    delete m_featureSpace;
    
    m_featureSpace = NULL;
    
    delete m_featureSpaceIndex;
    
    m_featureSpaceIndex = NULL;
}

template<class T>
void FSTestBase<T>::generate_featurespace() 
{
    cout << "Creating featurespace ... ";
    
    const map<int,double> lower_thresholds,upper_thresholds,fill_values;
    
    vector<string> variable_names;
    for (size_t i=0; i<this->m_variables.size(); i++)
        variable_names.push_back(this->m_variables.at(i).getName());
    
    NetCDFDataStore<T> *dataStore = new NetCDFDataStore<T>(this->m_filename,
                                                           this->coordinate_system(),
                                                           variable_names,0);
    
    this->m_featureSpace = new FeatureSpace<T>(this->coordinate_system(),
                                               dataStore,
                                               lower_thresholds,
                                               upper_thresholds,
                                               fill_values );
    
    this->m_featureSpaceIndex = PointIndex<T>::create(this->m_featureSpace->get_points(),
                                                      this->m_featureSpace->rank());
    
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
