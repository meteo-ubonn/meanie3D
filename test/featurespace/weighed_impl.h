#ifndef M3D_TEST_FS_WEIGHED_IMPL_H
#define M3D_TEST_FS_WEIGHED_IMPL_H

template <class T> 
void FSWeighedTest2D<T>::create_uniform_distribution_recursive( const NcVar &var,
                                                                size_t modulo,
                                                                size_t dimensionIndex,
                                                                typename CoordinateSystem<T>::GridPoint &gridpoint )
{
    using namespace netCDF;
    
    NcDim dim = var.getDim(dimensionIndex);
    
    if ( dimensionIndex < ( var.getDimCount() - 1 ) )
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
            
            vector<size_t> count( gridpoint.size(), 1 );
            
            if ( gridpoint[dimensionIndex] % modulo == 0 )
            {
                this->m_pointCount++;
                
                typename CoordinateSystem<T>::Coordinate coordinate = this->coordinate_system()->newCoordinate();
                
                this->coordinate_system()->lookup( gridpoint, coordinate );
                
                T value = (T) m_distribution( coordinate );
                
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
void FSWeighedTest2D<T>::create_uniform_distribution( const NcVar &var, size_t modulo  )
{
    // homogenous point distribution in N-D, every modulo grid point is used
    
    typename CoordinateSystem<T>::GridPoint gridpoint( this->coordinate_system()->rank(), 0 );
    
    this->m_pointCount = 0;
    
    cout << "Creating uniform distribution with data every " << modulo << " nodes ...";
    
    create_uniform_distribution_recursive( var, modulo, 0, gridpoint );
    
    cout << "done. (" << this->m_pointCount << " points)" << endl;
    
    this->m_totalPointCount += this->m_pointCount;
}

template<class T> 
void FSWeighedTest2D<T>::SetUp()
{
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
    
    this->generate_dimensions();

    // Create a variable
    
    // test case 1 : ellipsis for unweighed sample mean
    
    NcVar var = this->add_variable( "weighed_test", 0.0, 3 * FS_VALUE_MAX );
    
    create_uniform_distribution( var, 1 );
    
    FSTestBase<T>::generate_featurespace();
}

template<class T> 
void FSWeighedTest2D<T>::TearDown()
{
    FSTestBase<T>::TearDown();
}

#pragma mark -
#pragma mark Test parameterization

template<class T> 
FSWeighedTest2D<T>::FSWeighedTest2D() : FSTestBase<T>()
{
    this->m_settings = new FSTestSettings( 2, 1, NUMBER_OF_GRIDPOINTS, FSTestBase<T>::filename_from_current_testcase() );
}

template<class T> 
FSWeighedTest3D<T>::FSWeighedTest3D() : FSWeighedTest2D<T>()
{
    this->m_settings = new FSTestSettings( 3, 1, NUMBER_OF_GRIDPOINTS, FSTestBase<T>::filename_from_current_testcase() );
}

// 2D
#if RUN_2D
TYPED_TEST_CASE( FSWeighedTest2D, DataTypes );

TYPED_TEST( FSWeighedTest2D, FS_WeighedSample_2D_Test ) 
{
    vector<TypeParam> x( this->m_settings->fs_dim(), 0 );
    x[ x.size()-1 ] = FS_VALUE_MAX;
    
    // G'is a kernel
    
    Kernel<TypeParam> *kernel = new EpanechnikovKernel<TypeParam>( 1.0f );
    
    vector< vector<TypeParam> > origins;
    
    vector< vector<TypeParam> > vectors;
    
    // iterate over the bandwidths
    
    MeanshiftOperation<TypeParam> op( this->m_featureSpace, this->m_featureSpaceIndex );
    
    for ( size_t i = 0; i < this->m_bandwidths.size(); i++ )
    {
        vector<TypeParam> h = this->m_bandwidths.at(i);
        
        RangeSearchParams<TypeParam> *params = new RangeSearchParams<TypeParam>(h);
        
        vector<TypeParam> m = op.meanshift( x, params, kernel, 0 );
        
        origins.push_back( x );
        
        vectors.push_back( m );
        
        // calculate only the spatial shift
        
        size_t dims = this->m_coordinate_system->rank();
        
        vector<TypeParam> c(dims);
        
        for ( size_t k = 0; k < dims; k++ )
        {
            c[k] = m[k];
        }
        
        // spatial shift should be within grid resolution close to zero.
        
        EXPECT_NEAR( vector_norm(c), 0.0, this->coordinate_system()->resolution_norm() );
    }
    
    delete kernel;
}
#endif

// 3D
#if RUN_3D

TYPED_TEST_CASE( FSWeighedTest3D, DataTypes );

TYPED_TEST( FSWeighedTest3D, FS_WeighedSample_3D_Test ) 
{
    vector<TypeParam> x( this->m_settings->fs_dim(), 0 );
    x[ x.size()-1 ] = FS_VALUE_MAX;
    
    // G'is a kernel
    
    Kernel<TypeParam> *kernel = new GaussianNormalKernel<TypeParam>( 1.0f );
    
    vector< vector<TypeParam> > origins;
    
    vector< vector<TypeParam> > vectors;
    
    // iterate over the bandwidths
    
    MeanshiftOperation<TypeParam> op( this->m_featureSpace, this->m_featureSpaceIndex );
    
    for ( size_t i = 0; i < this->m_bandwidths.size(); i++ )
    {
        vector<TypeParam> h = this->m_bandwidths.at(i);
        
        RangeSearchParams<TypeParam> *params = new RangeSearchParams<TypeParam>(h);
        
        vector<TypeParam> m = op.meanshift( x, params, kernel, 0 );
        
        origins.push_back( x );
        
        vectors.push_back( m );
        
        // calculate only the spatial shift
        
        size_t dims = this->m_coordinate_system->rank();
        
        vector<TypeParam> c(dims);
        
        for ( size_t k = 0; k < dims; k++ )
        {
            c[k] = m[k];
        }
        
        // spatial shift should be within grid resolution close to zero.
        
        EXPECT_NEAR( vector_norm(c), 0.0, this->coordinate_system()->resolution_norm() );
    }
    
    delete kernel;
}
#endif


#endif

