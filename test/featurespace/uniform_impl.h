#ifndef M3D_TEST_FS_UNIFORM_IMPL_H
#define M3D_TEST_FS_UNIFORM_IMPL_H

#include <typeinfo>

template <class T> 
void 
FSUniformTest2D<T>::create_uniform_distribution_recursive( const NcVar &var,
                                                           size_t modulo,
                                                           size_t dimensionIndex, 
                                                           typename CoordinateSystem<T>::GridPoint &gridpoint )
{
    NcDim dim = var.getDim(dimensionIndex);
    
    if ( dimensionIndex < this->file()->getDimCount() - 1 )
    {
        for ( int index = 0; index < dim.getSize(); index++ )
        {
            gridpoint[dimensionIndex] = index;
            
            if ( gridpoint[dimensionIndex] % modulo == 0 )
            {
                create_uniform_distribution_recursive( var, modulo, dimensionIndex+1, gridpoint );
            }
        }
    }
    else
    {
        for ( int index = 0; index < dim.getSize(); index++ )
        {
            gridpoint[dimensionIndex] = index;
            
            if ( gridpoint[dimensionIndex] % modulo == 0 )
            {
                this->m_pointCount++;
                
                T value = (T) FS_VALUE_MAX;
                
                vector<size_t> gp(gridpoint.begin(), gridpoint.end());
                
                var.putVar( gp, value );
            }
        }
    }
}

/** Creates a distribution of points every modulo grid points, 
 * where all values are 1.0
 */
template <class T> 
void 
FSUniformTest2D<T>::create_uniform_distribution( const NcVar &var, size_t modulo  )
{
    // homogenous point distribution in N-D, every modulo grid point is used
    
    typename CoordinateSystem<T>::GridPoint gridpoint( this->file()->getDimCount(), 0 );
    
    this->m_pointCount = 0;
    
    INFO << "Creating uniform distribution with data every " << modulo << " nodes ...";
    
    create_uniform_distribution_recursive( var, modulo, 0, gridpoint );
    
    if (INFO_ENABLED) cout << "done. (" << this->m_pointCount << " points)" << endl;
    
    this->m_totalPointCount += this->m_pointCount;
}

template<class T> 
void FSUniformTest2D<T>::SetUp()
{
    const ::testing::TestInfo* const test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    
    INFO << "Setting up test " << test_info->test_case_name() << " with typeid " << typeid(T).name() << endl;
    
    this->m_settings = new FSTestSettings( 2, 1, NUMBER_OF_GRIDPOINTS, FSTestBase<T>::filename_from_current_testcase() );
    
    FSTestBase<T>::SetUp();
    
    // Set the bandwidths
    
    size_t bw_size = this->m_settings->fs_dim();
    
    this->m_bandwidths.push_back( vector<T>( bw_size, 1.0 ) );
    this->m_bandwidths.push_back( vector<T>( bw_size, 2.0 ) );
    this->m_bandwidths.push_back( vector<T>( bw_size, 3.0 ) );
    this->m_bandwidths.push_back( vector<T>( bw_size, 4.0 ) );
    
    // Figure out the bandwidth maximum
    
    vector<float> max_h( this->m_settings->fs_dim(), 0.0 );
    
    for ( size_t i = 0; i < this->m_bandwidths.size(); i++ )
    {
        vector<T> h = this->m_bandwidths.at(i);
        
        assert( h.size() == bw_size );
        
        for ( size_t j = 0; j < bw_size; j++ )
        {
            assert( h[j] > 0 );
            
            if ( h[j] > max_h[j] )
            {
                max_h[j] = h[j];
            }
        }
    }
    
    // Make the actual values 1.1 times bigger, so that in graphical
    // representations there is some boundary
    
    max_h *= 1.1f;

    this->m_settings->set_axis_bound_values( max_h );
    
    // Generate dimensions and dimension variables according to
    // the current settings
    
    try
    {
        this->generate_dimensions();
    }
    catch (const netCDF::exceptions::NcException &e)
    {
        cerr << "FATAL:error while generating dimensions: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Create a variable
    
    // test case 1 : ellipsis for unweighed sample mean
    
    NcVar var = this->add_variable( "uniform_sample", 0.0, FS_VALUE_MAX );
    
    create_uniform_distribution( var, 5 );
    
    FSTestBase<T>::generate_featurespace();
}

template<class T> 
void FSUniformTest2D<T>::TearDown()
{
    FSTestBase<T>::TearDown();
}

#pragma mark -
#pragma mark Test parameterization

template<class T> 
FSUniformTest2D<T>::FSUniformTest2D() : FSTestBase<T>()
{
}

template<class T> 
FSUniformTest3D<T>::FSUniformTest3D() : FSUniformTest2D<T>()
{
    std::string filename = FSTestBase<T>::filename_from_current_testcase();
    this->m_settings = new FSTestSettings( 3, 1, NUMBER_OF_GRIDPOINTS, filename );
}

// 2D
#if RUN_2D

TYPED_TEST_CASE( FSUniformTest2D, DataTypes );

TYPED_TEST( FSUniformTest2D, FS_UniformSampleTest_2D ) 
{
    // origin
    
    vector<TypeParam> x( this->m_settings->fs_dim(), 0 );
    x[ x.size()-1 ] = FS_VALUE_MAX;
    
    // G'is a kernel
    
    Kernel<TypeParam> *kernel = new GaussianNormalKernel<TypeParam>( 1.0f );
    
    // iterate over the bandwidths
    
    MeanshiftOperation<TypeParam> op( this->m_featureSpace, this->m_featureSpaceIndex );
    
    for ( size_t i = 0; i < this->m_bandwidths.size(); i++ )
    {
        vector<TypeParam> h = this->m_bandwidths.at(i);
        
        RangeSearchParams<TypeParam> *params = new RangeSearchParams<TypeParam>(h);
        
        vector<TypeParam> m = op.meanshift( x, params, kernel );
        
        EXPECT_NEAR( vector_norm(m), 0.0, this->coordinate_system()->resolution_norm() );
    }
    
    delete kernel;
}
#endif

// 3D
#if RUN_3D

TYPED_TEST_CASE( FSUniformTest3D, DataTypes );

TYPED_TEST( FSUniformTest3D, FS_UniformSampleTest_3D )
{
    // origin
    
    vector<TypeParam> x( this->m_settings->fs_dim(), 0 );
    x[ x.size()-1 ] = FS_VALUE_MAX;
    
    // G'is a kernel
    
    Kernel<TypeParam> *kernel = new GaussianNormalKernel<TypeParam>( 1.0f );
    
    // iterate over the bandwidths
    
    MeanshiftOperation<TypeParam> op( this->m_featureSpace, this->m_featureSpaceIndex );
    
    for ( size_t i = 0; i < this->m_bandwidths.size(); i++ )
    {
        vector<TypeParam> h = this->m_bandwidths.at(i);
        
        RangeSearchParams<TypeParam> *params = new RangeSearchParams<TypeParam>(h);
        
        vector<TypeParam> m = op.meanshift( x, params, kernel );
        
        EXPECT_NEAR( vector_norm(m), 0.0, this->coordinate_system()->resolution_norm() );
    }
    
    delete kernel;
}
#endif

#endif

