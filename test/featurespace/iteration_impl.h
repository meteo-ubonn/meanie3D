#ifndef M3D_TEST_FS_ITERATION_IMPL_H
#define M3D_TEST_FS_ITERATION_IMPL_H

#include <cf-algorithms/utils/rand_utils.h>
#include <cf-algorithms/utils/array_utils.h>

template<class T> 
void FSIterationTest2D<T>::write_cloud( const NcVar &variable )
{
    using cfa::utils::random::box_muller;
    
    // allocate a cursor
    
    vector<size_t> cursor( this->coordinate_system()->size(), 0 );
    
    // start generating random points (with random values between 0 and 1)
    
    size_t numPoints = 0;
    
    do 
    {
        // genrate a random coordinate
        
        for ( size_t d = 0; d < this->coordinate_system()->size(); d++ )
        {
            bool valid = false;
            
            while ( !valid )
            {
                NcDim dim = this->coordinate_system()->dimensions()[d];
                
                NcVar var = this->m_file->getVar( dim.getName() );
                
                float rand = box_muller( m_mean->at(d), m_deviation->at(d) );
                
                // re-transform to grid coordinates
                
                float min,max;
                
                var.getAtt("valid_min").getValues( &min );
                var.getAtt("valid_max").getValues( &max );
                
                long n = (long)round( (dim.getSize()-1)*(rand - min) / ( max - min ) );
                
                if ( n >=0 && n < dim.getSize() )
                {
                    cursor[d] = n;
                    
                    valid = true;
                }
            }
        }
        
        // generate a random value
        
        T value = (T) FS_VALUE_MAX;
        
        variable.putVar( cursor, value );
        
        numPoints++;
        
    } while ( numPoints < m_cloudSize );
}

template<class T> 
void FSIterationTest2D<T>::SetUp()
{
    FSTestBase<T>::SetUp();
    
    // Generate dimensions and dimension variables according to
    // the current settings
    
    this->generate_dimensions();
    
    // Create a variable
    
    // test case 1 : ellipsis for unweighed sample mean
    
    NcVar var = this->add_variable( "iteration_test", 0.0, FS_VALUE_MAX );    

    write_cloud( var );
    
    FSTestBase<T>::generate_featurespace();
}

template<class T> 
void FSIterationTest2D<T>::TearDown()
{
    delete m_mean;
    
    delete m_center;
    
    delete m_deviation;
    
    FSTestBase<T>::TearDown();
}

#pragma mark -
#pragma mark Test parameterization

template<class T> 
FSIterationTest2D<T>::FSIterationTest2D() : m_center(NULL), m_mean(NULL), m_deviation(NULL), m_cloudSize(5000)
{
    this->m_settings = new FSTestSettings( 2, 1, NUMBER_OF_GRIDPOINTS, FSTestBase<T>::filename_from_current_testcase() );
    
    this->m_mean = new vector<T>( 3, 0.0 );
    
    this->m_deviation = new vector<T>( 3, 0.25 );
    
    this->m_center = new vector<T>( 3, 0.0 );

    // Run the iteration tests with a series of bandwidths
    
    size_t bw_size = this->m_settings->fs_dim();
    
    for ( int i = 1; i <= 5; i++ )
    {
        T radius = 0.25 + static_cast<float>(i) * 0.1;
        
        // increasing bandwidth
        
        vector<T> h( bw_size, radius );
        
        this->m_bandwidths.push_back( h );
        
        // random starting points in relation to the bandwidth 
        // parameter. The further outward the less dense the points 
        // in the cloud, so I associate larger bandwidths with 
        // points lying out more
        
        vector<T> origin( bw_size );

        vector<T> coord( 2 );

        // spiral outwards
        float alpha = static_cast<float>(i) * 2.0 * M_PI / 5.0;
        
        coord[0] = origin[0] = radius * cos( alpha );
        
        coord[1] = origin[1] = radius * sin( alpha );
        
        origin[2] = FS_VALUE_MAX;
        
        cfa::meanshift::Point<T> point( coord, origin );
        
        this->m_origins.push_back( point );
    }
}

template<class T> 
FSIterationTest3D<T>::FSIterationTest3D() 
{
    // The superconstructor already put some stuff
    // up, so get rid of it.
    // TODO: perhaps move this shit into SetUp() ?
    
    this->m_bandwidths.clear();
    
    this->m_origins.clear();
    
    // Settings
    
    this->m_settings = new FSTestSettings( 3, 1, NUMBER_OF_GRIDPOINTS, FSTestBase<T>::filename_from_current_testcase() );

    this->m_cloudSize = 10000;

    this->m_mean = new vector<T>( 4, 0.0 );
    
    this->m_deviation = new vector<T>( 4, 0.25 );
    
    this->m_center = new vector<T>( 4, 0.0 );

    // Run the iteration tests with a series of bandwidths
    
    size_t bw_size = this->m_settings->fs_dim();
    
    for ( int i = 1; i <= 5; i++ )
    {
        float radius = static_cast<float>(i) * 0.2;
        
        // increasing bandwidth
        
        vector<T> h( bw_size, radius );
        
        this->m_bandwidths.push_back( h );
        
        // random starting points in relation to the bandwidth 
        // parameter. The further outward the less dense the points 
        // in the cloud, so I associate larger bandwidths with 
        // points lying out more
        
        vector<T> origin( bw_size, 0.0 );
        
        vector<T> coord( 2 );

        float theta = static_cast<float>(i) * M_2_PI / 5.0; 

        float phi = static_cast<float>(i) * M_PI_2 / 5.0; 

        coord[0] = origin[0] = radius * cos( phi ) * sin( theta );
        
        coord[1] = origin[1] = radius * sin( phi ) * sin( theta );
        
        coord[2] = origin[2] = radius * cos( phi );
        
        origin[3] = FS_VALUE_MAX;
        
        cfa::meanshift::Point<T> point( coord, origin );
        
        this->m_origins.push_back( point );
    }
}
// 2D
#if RUN_2D

TYPED_TEST_CASE( FSIterationTest2D, DataTypes );

TYPED_TEST( FSIterationTest2D, FS_Iteration_2D_Test ) 
{
    using cfa::meanshift::GaussianNormalKernel;
    using cfa::meanshift::IterationOperation;
    using cfa::meanshift::RangeSearchParams;
    using cfa::meanshift::vector_norm;
    using cfa::utils::VisitUtils;
    
#if WRITE_FEATURESPACE
    const ::testing::TestInfo* const test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    string fs_filename = string(test_info->test_case_name()) + "_featurespace.vtk";
    boost::replace_all( fs_filename, "/", "_" );
    cout << "Writing Featurespace to " + fs_filename << " ... ";
    VisitUtils<TypeParam>::write_featurespace_vtk( fs_filename, this->m_featureSpace );
    cout << "done." << endl;
#endif
    
    // G'is a kernel
    
    GaussianNormalKernel<TypeParam> kernel( 1.0f );
    
    cout << setiosflags(ios::fixed) << setprecision(TEST_PRINT_PRECISION);
    
    for ( size_t i = 0; i < this->m_origins.size(); i++ )
    {
        //cout << "Entering iteration at " << this->m_origins[i].values << " with bandwidths " << this->m_bandwidths[i] << " ... " << endl;
        
        typename FeatureSpace<TypeParam>::Trajectory *trajectory = NULL;
        
        IterationOperation<TypeParam> op( this->m_featureSpace, this->m_featureSpaceIndex );
        
        RangeSearchParams<TypeParam> *params = new RangeSearchParams<TypeParam>( this->m_bandwidths[i] );
        
        trajectory = op.iterate( &this->m_origins[i],
                                params,
                                &kernel, 
                                NULL,
                                TERMCRIT_EPSILON,
                                TERMCRIT_ITER);
        
        vector<TypeParam> last = trajectory->back();
        
        //cout << "done. (" << trajectory->size() << " points, ending at " << last << endl;
        
        // Results should lie around 0 with tolerance termcrit_epsilon
        EXPECT_NEAR( vector_norm(last), FS_VALUE_MAX, this->coordinate_system()->resolution_norm() );
        
        // There should be no more than termcrit_iter points in the 
        // trajectory
        EXPECT_LT( trajectory->size(), TERMCRIT_ITER );
        
#if WRITE_TRAJECTORIES
        
        // Write trajectory out
        const ::testing::TestInfo* const test_info = ::testing::UnitTest::GetInstance()->current_test_info();
        string fn = std::string(test_info->test_case_name()) + "_trajectory_" + boost::lexical_cast<string>( i ) + ".vtk";
        boost::replace_all( fn, "/", "_" );
        VisitUtils<TypeParam>::write_pointlist_vtk( fn, trajectory, "trajectory" );
        
#endif
        
#if WRITE_BANDWIDTH
        vector<TypeParam> start = trajectory->front();
        size_t dim = trajectory->front().size();
        string bw_fn = "trajectory_bw_" + boost::lexical_cast<string>(i) + (dim == 2 ? ".curve" : ".3D" );
        if ( dim == 2 )
        {
            VisitUtils<TypeParam>::write_ellipsis_2d( bw_fn.c_str(), this->m_bandwidths[i], 1000, &start );
        }
        else if ( dim == 3 )
        {
            VisitUtils<TypeParam>::write_ellipsis_3d( bw_fn.c_str(), this->m_bandwidths[i], 1000, &start );
        }
#endif   
        // clean up
        
        delete trajectory;
    }
}
#endif

// 3D
#if RUN_3D

TYPED_TEST_CASE( FSIterationTest3D, DataTypes );

TYPED_TEST( FSIterationTest3D, FS_Iteration_3D_Test ) 
{
    using cfa::meanshift::IterationOperation;
    using cfa::meanshift::GaussianNormalKernel;
    using cfa::meanshift::RangeSearchParams;
    using cfa::meanshift::vector_norm;
    using cfa::utils::VisitUtils;
    
#if WRITE_FEATURESPACE
    const ::testing::TestInfo* const test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    string fs_filename = string(test_info->test_case_name()) + "_featurespace.vtk";
    boost::replace_all( fs_filename, "/", "_" );
    cout << "Writing Featurespace to " + fs_filename << " ... ";
    VisitUtils<TypeParam>::write_featurespace_vtk( fs_filename, this->m_featureSpace );
    cout << "done." << endl;
#endif

    // G'is a kernel
    
    GaussianNormalKernel<TypeParam> kernel( 1.0f );
    
    cout << setiosflags(ios::fixed) << setprecision(TEST_PRINT_PRECISION);
    
#if WRITE_ITERATION_ORIGINS
    static size_t test_count = 0;
    string fn = "origins_" + boost::lexical_cast<string>(test_count++) + ".3D";
    ofstream f(fn.c_str());
    f << "x\ty\tz\tv" << endl;
    f << fixed << setprecision(4);
    
    cout << "Starting points:" << endl;
    for ( size_t i = 0; i < this->m_origins.size(); i++ )
    {
        cout << this->m_origins[i].values << endl;
        for ( size_t j=0; j < this->m_origins[i].values.size(); j++ )
        {
            f << this->m_origins[i].values[j] << "\t";
        }
        f << endl;
        cout << endl;
    }
    f.close();
#endif
    
    for ( size_t i = 0; i < this->m_origins.size(); i++ )
    {
        typename FeatureSpace<TypeParam>::Trajectory *trajectory = NULL;
        
        //cout << "Entering iteration at " << this->m_origins[i].values << " with bandwidths " << this->m_bandwidths[i] << " ... " << endl;
        
        IterationOperation<TypeParam> op( this->m_featureSpace, this->m_featureSpaceIndex );
        
        RangeSearchParams<TypeParam> *params = new RangeSearchParams<TypeParam>( this->m_bandwidths[i] );
        
        trajectory = op.iterate( &this->m_origins[i],
                                params,
                                &kernel, 
                                NULL,
                                TERMCRIT_EPSILON,
                                TERMCRIT_ITER);
        
        vector<TypeParam> last = trajectory->back();

        //cout << "done. (" << trajectory->size() << " points, ending at " << last << endl;
         
        // Results should lie around 0 with tolerance termcrit_epsilon
        EXPECT_NEAR( vector_norm(last), FS_VALUE_MAX, this->coordinate_system()->resolution_norm() );
        
        // There should be no more than termcrit_iter points in the 
        // trajectory
        EXPECT_LT( trajectory->size(), TERMCRIT_ITER );
        
        
#if WRITE_TRAJECTORIES
        
        // Write trajectory out
        
        string fn = "trajectory_" + boost::lexical_cast<string>( i ) + ".vtk";
        
        VisitUtils<TypeParam>::write_pointlist_vtk( fn, trajectory, "trajectory" );
#endif
        
#if WRITE_BANDWIDTH
        
        typename FeatureSpace<TypeParam>::TrajectoryIterator ti;
        
        size_t iter_point = 0;
        
        for ( ti = trajectory->begin(); ti != trajectory->end(); ti++ )
        {
            vector<TypeParam> p = *ti;
            
            string bw_fn = "bandwidth_" + boost::lexical_cast<string>( iter_point++ ) + ".3D";
            
            VisitUtils<TypeParam>::write_ellipsis_3d( bw_fn.c_str(), this->m_bandwidths[i], 250, &p );
        }
#endif   
        
        // clean up
        
        delete trajectory;
    }
}
#endif

#endif

