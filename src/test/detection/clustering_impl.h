#ifndef _M3D_TEST_FS_CLUSTERING_IMPL_H_
#define _M3D_TEST_FS_CLUSTERING_IMPL_H_

#include <meanie3D/meanie3D.h>

template<class T>
void FSClusteringTest2D<T>::write_cloud( const NcVar &var, vector<T> mean, vector<T> deviation )
{
    using namespace netCDF;
    using namespace m3D;
    using namespace cfa::utils::random;
    
    // start generating random points (with random values between 0 and 1)
    
    size_t numPoints = 0;
    
    // Use a simple gaussian normal function to simulate a value distribution
    
    GaussianNormal<T> gauss;
    
    T gauss_zero = gauss( vector<T>(this->coordinate_system()->size(),0) );
    
    typename CoordinateSystem<T>::GridPoint gridpoint = this->coordinate_system()->newGridPoint();
    
    do
    {
        // genrate a random coordinate
        
        for ( size_t d = 0; d < this->coordinate_system()->size(); d++ )
        {
            bool valid = false;
            
            while ( !valid )
            {
                NcDim dim = this->coordinate_system()->dimensions()[d];
                
                NcVar dimVar = this->coordinate_system()->dimension_variable( dim );
                
                T min, max;
                
                dimVar.getAtt("valid_min").getValues( &min );
                
                dimVar.getAtt("valid_max").getValues( &max );
                
                float rand = box_muller( mean[d], deviation[d] );
                
                // re-transform to grid coordinates
                
                size_t n = (size_t) round( (dim.getSize()-1) * (rand - min) / ( max - min ) );
                
                if ( n < dim.getSize() )
                {
                    gridpoint[d] = n;
                    
                    valid = true;
                }
            }
        }
        
        // value goes to
        
        vector<T> x = this->coordinate_system()->newCoordinate();
        
        this->coordinate_system()->lookup( gridpoint, x );
        
        T value = (T) (FS_VALUE_MAX * gauss( x-mean )) / gauss_zero;
        
        // cout << "x=" << x << " mean=" << mean << " g(x-m)=" << value << endl;
        
        var.putVar( gridpoint, value );
        
        numPoints++;
        
        this->m_pointCount++;
        
    } while ( numPoints < m_cloudSize );
}

template<class T>
void FSClusteringTest2D<T>::create_clouds_recursive( const NcVar &var, size_t dimensionIndex, typename CoordinateSystem<T>::GridPoint &gridpoint )
{
    using namespace netCDF;
    
    NcDim dim = var.getDim(dimensionIndex);
    
    size_t increment = this->m_division_increments[dim];
    
    if ( dimensionIndex < (var.getDimCount()-1))
    {
        for ( int index = 1; index < m_divisions; index++ )
        {
            gridpoint[dimensionIndex] = index * increment;
            
            create_clouds_recursive( var, dimensionIndex+1, gridpoint );
        }
    }
    else
    {
        for ( int index = 1; index < m_divisions; index++ )
        {
            gridpoint[dimensionIndex] = index * increment;
            
            // get the variables together and construct the cartesian coordinate
            // of the current point. If it's on the ellipse, put it in the variable
            
            vector<T> coordinate( var.getDimCount() );
            
            this->coordinate_system()->lookup( gridpoint, coordinate );
            
            cout << "\tWriting cloud at gridpoint=" <<  gridpoint << " coordinate=" << coordinate << " deviation=" << m_deviation << endl;
            
            write_cloud( var, coordinate, m_deviation );
        }
    }
}

template<class T>
void FSClusteringTest2D<T>::create_clouds( const NcVar &var )
{
    // calculate the divisions in terms of grid points
    // calculate the deviations
    
    vector<NcDim *>::iterator dim_iter;
    
    for ( size_t index = 0; index < this->coordinate_system()->size(); index++ )
    {
        NcDim dim = this->coordinate_system()->dimensions()[index];
        
        NcVar dim_var = this->coordinate_system()->dimension_variable( dim );
        
        size_t number_gridpoints = cfa::utils::netcdf::num_vals(dim_var) / m_divisions;
        
        m_division_increments[ dim ] = number_gridpoints;
        
        m_deviation.push_back( 0.4 / m_divisions );
    }
    
    typename CoordinateSystem<T>::GridPoint gridpoint = this->coordinate_system()->newGridPoint();
    
    this->m_pointCount = 0;
    
    cout << "Creating clusters at the intersection of " << m_divisions << " lines per axis ..." << endl;
    
    create_clouds_recursive( var, 0, gridpoint );
    
    cout << "done. (" << this->m_pointCount << " points)" << endl;
    
    this->m_totalPointCount += this->m_pointCount;
}

template<class T>
void FSClusteringTest2D<T>::SetUp()
{
    FSTestBase<T>::SetUp();
    
    // Set the bandwidths
    
    size_t bw_size = this->m_settings->fs_dim();
    
    for ( size_t i=0; i<bw_size; i++ )
    {
        m_fuzziness.push_back( 1.0 / m_divisions );
    }
    
    this->m_bandwidths.push_back( vector<T>( bw_size, 1.0 / m_divisions ) );
    
    // Generate dimensions and dimension variables according to
    // the current settings
    
    this->generate_dimensions();
    
    // Create a variable
    
    // test case 1 : ellipsis for unweighed sample mean
    
    NcVar var = this->add_variable( "cluster_test", 0.0, FS_VALUE_MAX );
    
    create_clouds( var );
    
    FSTestBase<T>::generate_featurespace();
}

template<class T>
void FSClusteringTest2D<T>::TearDown()
{
    FSTestBase<T>::TearDown();
}

#pragma mark -
#pragma mark Test parameterization

template<class T>
FSClusteringTest2D<T>::FSClusteringTest2D() : FSTestBase<T>()
{
    this->m_settings = new FSTestSettings( 2, 1, NUMBER_OF_GRIDPOINTS, FSTestBase<T>::filename_from_current_testcase() );
    
    this->m_divisions = 3;
    
    this->m_cloudSize = 50 * NUMBER_OF_GRIDPOINTS / m_divisions;
}

template<class T>
FSClusteringTest3D<T>::FSClusteringTest3D() : FSClusteringTest2D<T>()
{
    this->m_settings = new FSTestSettings( 3, 1, NUMBER_OF_GRIDPOINTS, FSTestBase<T>::filename_from_current_testcase() );
    
    this->m_divisions = 4;
    
    this->m_cloudSize = 100 * NUMBER_OF_GRIDPOINTS / this->m_divisions;
}



// 2D
#if RUN_2D

TYPED_TEST_CASE( FSClusteringTest2D, DataTypes );

TYPED_TEST( FSClusteringTest2D, FS_Clustering_2D_Range_Test )
{
    using cfa::utils::VisitUtils;
    using namespace cfa::utils::timer;
    using namespace cfa::meanshift;
    using namespace m3D;

    GaussianNormal<TypeParam> gauss;
    // TypeParam gauss_zero = gauss( vector<TypeParam>(this->coordinate_system()->size(),0) );

    cout << setiosflags(ios::fixed) << setprecision(TEST_PRINT_PRECISION);
    
    for ( size_t i = 0; i < this->m_bandwidths.size(); i++ )
    {
        vector<TypeParam> h = this->m_bandwidths.at(i);
        
        Kernel<TypeParam> *kernel = new GaussianNormalKernel<TypeParam>( vector_norm(h) );
        
        cout << "Clustering feature space with"
        << " epsilon=" << this->coordinate_system()->resolution_norm()
        << " max_iter=" << TERMCRIT_ITER << " ... "
        << " with bandwidths " << h  << " ... " << endl;
        
        typename Cluster<TypeParam>::list::iterator ci;
        
        start_timer();
        
        //        VariableWeighed<TypeParam> *weight = new VariableWeighed<TypeParam>( this->m_file, this->coordinate_system, this->m_variables.front() );
        
        
        // Create 'bandwidth' parameter from grid resolution and
        // the maximum value in the value range
        
        vector<TypeParam> resolution = ((TypeParam)3.0) * this->coordinate_system()->resolution();
        
        resolution.push_back( numeric_limits<TypeParam>::max());
        
        // Search parameters for neighbourhood searches in clustering
        // should be in the order of the grid resolution
        
        RangeSearchParams<TypeParam> *params = new RangeSearchParams<TypeParam>( resolution );
        
        // KNNSearchParams *params = new KNNSearchParams( 200 );
        
        ClusterOperation<TypeParam> op( this->m_featureSpace, this->m_featureSpaceIndex );
        
        ClusterList<TypeParam> clusters = op.cluster( params, kernel, 0 );
        
        delete kernel;
        
        //        delete weight;
        
        //        // Create the cluster graph
        //        this->m_featureSpace->create_cluster_graph( h, kernel, NULL );
        
        double time = stop_timer();
        cout << "done. (" << time << " seconds, " << clusters.clusters.size() << " modes found:" << endl;
        size_t cluster_number = 1;
        for ( ci = clusters.clusters.begin(); ci != clusters.clusters.end(); ci++ )
        {
            Cluster<TypeParam> *c = *ci;
            
            cout << "Cluster #" << cluster_number++ << " at " << c->mode << " (" << c->points.size() << " points.)" << endl;
        }
        
        const ::testing::TestInfo* const test_info = ::testing::UnitTest::GetInstance()->current_test_info();
        
#if WRITE_FEATURESPACE
        string fs_filename = string(test_info->test_case_name()) + "_featurespace.vtk";
        boost::replace_all( fs_filename, "/", "_" );
        cout << "Writing Featurespace to " + fs_filename << " ... ";
        VisitUtils<TypeParam>::write_featurespace_vtk( fs_filename, this->m_featureSpace );
        cout << "done." << endl;
#endif
        
#if WRITE_MEANSHIFT_VECTORS
        
        // compile and write vectors
        
        vector< vector<TypeParam> > origins, vectors;
        
        for ( size_t index = 0; index < this->m_featureSpace->points.size(); index++ )
        {
            Point<TypeParam> *p = this->m_featureSpace->points[ index ];
            vector<TypeParam> origin = p->values;
            origins.push_back( origin );
            //
            //            vector<TypeParam> shift( origin.size(), 0 );
            //
            //            if ( p->next )
            //            {
            //                vector_subtract_r( origin, p->next->values, shift );
            //            }
            //
            //            vectors.push_back( shift );
            
            vectors.push_back( p->shift );
            
        }
        
        string vectors_filename = string(test_info->test_case_name()) + "_meanshift.vtk";
        boost::replace_all( vectors_filename, "/", "_" );
        cout << "Writing meanshift vectors to " << vectors_filename << " ... ";
        cfa::utils::VisitUtils<TypeParam>::write_vectors_vtk( vectors_filename, origins, vectors );
        cout << "done." << endl;
#endif
        
#if WRITE_CLUSTERS
        cout << "Writing clusters ... ";
        m3D::utils::VisitUtils<TypeParam>::write_clusters_vtk( test_info->test_case_name(), clusters.clusters, h );
        cout << "done." << endl;
#endif
        ClusterList<TypeParam>::reset_clustering(this->m_featureSpace);
        
        EXPECT_GE( 4, clusters.clusters.size() );
    }
}
#endif

// 3D
#if RUN_3D

TYPED_TEST_CASE( FSClusteringTest3D, DataTypes );

TYPED_TEST( FSClusteringTest3D, FS_Clustering_3D_Test )
{
    using namespace m3D;
    using namespace cfa::utils::timer;
    using cfa::utils::VisitUtils;
    
    const ::testing::TestInfo* const test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    
    cout << setiosflags(ios::fixed) << setprecision(TEST_PRINT_PRECISION);
    
    for ( size_t i = 0; i < this->m_bandwidths.size(); i++ )
    {
        vector<TypeParam> h = this->m_bandwidths.at(i);
        
        Kernel<TypeParam> *kernel = new GaussianNormalKernel<TypeParam>( vector_norm(h) );
        
        cout << "Clustering feature space with"
        << " epsilon=" << this->coordinate_system()->resolution_norm()
        << " max_iter=" << TERMCRIT_ITER << " ... "
        << " with bandwidths " << h << " ... " << endl;
        
        typename Cluster<TypeParam>::list::iterator ci;
        
        start_timer();
        
        //        VariableWeighed<TypeParam> *weight = new VariableWeighed<TypeParam>( this->m_file, this->coordinate_system, this->m_variables.front() );
        ClusterOperation<TypeParam> op( this->m_featureSpace, this->m_featureSpaceIndex );
        
        // Create 'bandwidth' parameter from grid resolution and
        // the maximum value in the value range
        
        vector<TypeParam> resolution = ((TypeParam)5.0) * this->coordinate_system()->resolution();
        
        resolution.push_back( numeric_limits<TypeParam>::max());
        
        // Search parameters for neighbourhood searches in clustering
        // should be in the order of the grid resolution
        
        RangeSearchParams<TypeParam> *params = new RangeSearchParams<TypeParam>( resolution );
        
        ClusterList<TypeParam> clusters = op.cluster( params, kernel, NULL, PostAggregationMethodNone );
        
        delete kernel;
        //        delete weight;
        
        //        // Create the cluster graph
        //        this->m_featureSpace->create_cluster_graph( h, kernel, NULL );
        
        double time = stop_timer();
        cout << "done. (" << time << " seconds, " << clusters.clusters.size() << " modes found:" << endl;
        size_t cluster_number = 1;
        for ( ci = clusters.clusters.begin(); ci != clusters.clusters.end(); ci++ )
        {
            Cluster<TypeParam> *c = *ci;
            
            cout << "Cluster #" << cluster_number++ << " at " << c->mode << " (" << c->points.size() << " points.)" << endl;
        }
        
#if WRITE_FEATURESPACE
        string fs_filename = string(test_info->test_case_name()) + "_featurespace.vtk";
        boost::replace_all( fs_filename, "/", "_" );
        cout << "Writing Featurespace to " + fs_filename << " ... ";
        VisitUtils<TypeParam>::write_featurespace_vtk( fs_filename, this->m_featureSpace );
        cout << "done." << endl;
#endif
        
#if WRITE_MEANSHIFT_VECTORS
        
        // compile and write vectors
        
        vector< vector<TypeParam> > origins, vectors;
        
        for ( size_t index = 0; index < this->m_featureSpace->points.size(); index++ )
        {
            Point<TypeParam> *p = this->m_featureSpace->points[ index ];
            vector<TypeParam> origin = p->values;
            origins.push_back( origin );
            //
            //            vector<TypeParam> shift( origin.size(), 0 );
            //
            //            if ( p->next )
            //            {
            //                vector_subtract_r( origin, p->next->values, shift );
            //            }
            //
            //            vectors.push_back( shift );
            
            vectors.push_back( p->shift );
            
        }
        
        string vectors_filename = string(test_info->test_case_name()) + "_meanshift.vtk";
        boost::replace_all( vectors_filename, "/", "_" );
        cout << "Writing meanshift vectors to " << vectors_filename << " ... ";
        VisitUtils<TypeParam>::write_vectors_vtk( vectors_filename, origins, vectors );
        VisitUtils<TypeParam>::write_vectors_vtk( "blahblah.vtk", origins, vectors );
        cout << "done." << endl;
#endif
        
#if WRITE_CLUSTERS
        cout << "Writing clusters ... ";
        m3D::utils::VisitUtils<TypeParam>::write_clusters_vtk( test_info->test_case_name(), clusters.clusters, h );
        cout << "done." << endl;
#endif
        EXPECT_GE( 27, clusters.clusters.size() );
    }
}

#endif


#endif

