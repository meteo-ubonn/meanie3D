#ifndef M3D_TEST_FS_ITERATION_IMPL_H
#define M3D_TEST_FS_ITERATION_IMPL_H

template<class T>
void FSIterationTest2D<T>::write_cloud(const NcVar &variable) {
    using namespace m3D;
    using m3D::utils::box_muller;

    // allocate a cursor
    vector <size_t> cursor(this->coordinate_system()->rank(), 0);

    // start generating random points (with random values between 0 and 1)
    size_t numPoints = 0;
    INFO << "Writing a cloud:" << endl;
    do {
        // genrate a random coordinate
        for (size_t d = 0; d < this->coordinate_system()->rank(); d++) {
            bool valid = false;
            while (!valid) {
                NcDim dim = this->coordinate_system()->dimensions()[d];
                NcVar var = this->m_file->getVar(dim.getName());
                float rand = box_muller(m_mean->at(d), m_deviation->at(d));
                // re-transform to grid coordinates
                float min, max;
                var.getAtt("valid_min").getValues(&min);
                var.getAtt("valid_max").getValues(&max);
                long n = (long) round((dim.getSize() - 1) * (rand - min) / (max - min));
                if (n >= 0 && n < dim.getSize()) {
                    cursor[d] = n;
                    valid = true;
                }
            }
        }

        // generate a random value
        typename CoordinateSystem<T>::Coordinate coordinate;
        coordinate = this->coordinate_system()->newCoordinate();
        vector<int> gridpoint(coordinate.size(), 0);
        for (size_t i = 0; i < coordinate.size(); i++)
            gridpoint[i] = (int) cursor[i];
        this->coordinate_system()->lookup(gridpoint, coordinate);
        static float sigma_square = 0.5;
        T norm = vectors::vector_norm(coordinate);
        T value = FS_VALUE_MAX * (1.0 / (sqrt(2 * M_PI * sigma_square))) * exp(-0.5 * norm * norm / sigma_square);
        variable.putVar(cursor, value);
        numPoints++;

    } while (numPoints < m_cloudSize);
}

template<class T>
void FSIterationTest2D<T>::SetUp() {
    FSTestBase<T>::SetUp();

    // Generate dimensions and dimension variables according to
    // the current settings
    this->generate_dimensions();

    // Create a variable
    // test case 1 : ellipsis for unweighed sample mean
    NcVar var = this->add_variable("iteration_test", 0.0, FS_VALUE_MAX);
    write_cloud(var);
    FSTestBase<T>::generate_featurespace();
}

template<class T>
void FSIterationTest2D<T>::TearDown() {
    delete m_mean;
    delete m_center;
    delete m_deviation;
    FSTestBase<T>::TearDown();
}

#pragma mark -
#pragma mark Test parameterization

template<class T>
FSIterationTest2D<T>::FSIterationTest2D() : m_center(NULL), m_mean(NULL), m_deviation(NULL), m_cloudSize(5000) {
    this->m_settings = new FSTestSettings(2, 1, NUMBER_OF_GRIDPOINTS, FSTestBase<T>::filename_from_current_testcase());
    this->m_mean = new vector<T>(3, 0.0);
    this->m_deviation = new vector<T>(3, 0.25);
    this->m_center = new vector<T>(3, 0.0);

    // Run the iteration tests with a series of bandwidths
    size_t bw_size = this->m_settings->fs_dim();

    for (int i = 1; i <= 5; i++) {
        T radius = 0.25 + static_cast<float> (i) * 0.1;
        // increasing bandwidth
        vector <T> h(bw_size, radius);
        h[h.size() - 1] = 2 * FS_VALUE_MAX;
        this->m_bandwidths.push_back(h);

        // random starting points in relation to the bandwidth 
        // parameter. The further outward the less dense the points 
        // in the cloud, so I associate larger bandwidths with 
        // points lying out more
        vector <T> origin(bw_size);
        vector <T> coord(2);

        // spiral outwards
        float alpha = static_cast<float> (i) * 2.0 * M_PI / 5.0;
        coord[0] = origin[0] = radius * cos(alpha);
        coord[1] = origin[1] = radius * sin(alpha);
        origin[2] = FS_VALUE_MAX * box_muller(0, 1);
        Point <T> point(coord, origin);
        this->m_origins.push_back(point);
    }
}

template<class T>
FSIterationTest3D<T>::FSIterationTest3D() {
    // The superconstructor already put some stuff
    // up, so get rid of it.
    // TODO: perhaps move this shit into SetUp() ?
    this->m_bandwidths.clear();
    this->m_origins.clear();

    // Settings
    this->m_settings = new FSTestSettings(3, 1, NUMBER_OF_GRIDPOINTS, FSTestBase<T>::filename_from_current_testcase());
    this->m_cloudSize = 10000;
    this->m_mean = new vector<T>(4, 0.0);
    this->m_deviation = new vector<T>(4, 0.25);
    this->m_center = new vector<T>(4, 0.0);

    // Run the iteration tests with a series of bandwidths
    size_t bw_size = this->m_settings->fs_dim();

    for (int i = 1; i <= 5; i++) {
        float radius = static_cast<float> (i) * 0.2;
        // increasing bandwidth
        vector <T> h(bw_size, radius);
        this->m_bandwidths.push_back(h);
        // random starting points in relation to the bandwidth 
        // parameter. The further outward the less dense the points 
        // in the cloud, so I associate larger bandwidths with 
        // points lying out more
        vector <T> origin(bw_size, 0.0);
        vector <T> coord(2);
        float theta = static_cast<float> (i) * M_2_PI / 5.0;
        float phi = static_cast<float> (i) * M_PI_2 / 5.0;
        coord[0] = origin[0] = radius * cos(phi) * sin(theta);
        coord[1] = origin[1] = radius * sin(phi) * sin(theta);
        coord[2] = origin[2] = radius * cos(phi);
        origin[3] = FS_VALUE_MAX * box_muller(0, 1);
        Point <T> point(coord, origin);
        this->m_origins.push_back(point);
    }
}
// 2D
#if RUN_2D

TYPED_TEST_CASE(FSIterationTest2D, DataTypes);

TYPED_TEST(FSIterationTest2D, FS_Iteration_2D_Test)
{
    using namespace m3D::utils::vectors;

#if WRITE_FEATURESPACE
    const ::testing::TestInfo * const test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    string fs_filename = string(test_info->test_case_name()) + "_featurespace.vtk";
    boost::replace_all(fs_filename, "/", "_");
    INFO << "Writing Featurespace to " + fs_filename << " ... ";
    VisitUtils<TypeParam>::write_featurespace_vtk(fs_filename, this->m_featureSpace);
    if (INFO_ENABLED)
        cout << "done." << endl;
#endif

    // G'is a kernel
    GaussianNormalKernel<TypeParam> kernel(1.0f);
    cout << setiosflags(ios::fixed) << setprecision(TEST_PRINT_PRECISION);
    for (size_t i = 0; i < this->m_origins.size(); i++) {
        
        typename FeatureSpace<TypeParam>::Trajectory *trajectory = NULL;
        
        IterationOperation<TypeParam> op(this->m_featureSpace, this->m_featureSpaceIndex);
        
        RangeSearchParams<TypeParam> *searchParams 
                = new RangeSearchParams<TypeParam>(this->m_bandwidths[i]);

        detection_params_t<TypeParam> params = Detection<TypeParam>::defaultParams();
        params.variables = this->m_variables;
        params.dimensions = this->m_dimensions;
        params.dimension_variables = this->m_dimension_variables;
        params.min_cluster_size = 20;
               
        detection_context_t<TypeParam> ctx;
        Detection<TypeParam>::initialiseContext(ctx);
        
        Detection<TypeParam>::initialiseContext(ctx);
        ctx.fs = this->m_featureSpace;
        ctx.coord_system = this->m_coordinate_system;
        ctx.data_store = this->m_data_store;

        WeightFunction<TypeParam> *weight 
                = new DefaultWeightFunction<TypeParam>(params,ctx);

        trajectory = op.get_trajectory(&this->m_origins[i],
                searchParams,
                &kernel,
                weight,
                TERMCRIT_EPSILON,
                TERMCRIT_ITER);

        // Results should lie around 0 with tolerance termcrit_epsilon

        vector<TypeParam> x2 = trajectory->at(trajectory->size() - 1);
        vector<TypeParam> x1 = trajectory->at(trajectory->size() - 2);
        EXPECT_NEAR(vector_norm(x2 - x1), 0, TERMCRIT_EPSILON);

        // There should be no more than termcrit_iter points in the 
        // trajectory
        EXPECT_LT(trajectory->size(), TERMCRIT_ITER);

        // clean up
        delete weight;
        delete trajectory;
    }
}
#endif

// 3D
#if RUN_3D

TYPED_TEST_CASE(FSIterationTest3D, DataTypes);

TYPED_TEST(FSIterationTest3D, FS_Iteration_3D_Test)
{
    using namespace m3D::utils::vectors;

#if WRITE_FEATURESPACE
    const ::testing::TestInfo * const test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    string fs_filename = string(test_info->test_case_name()) + "_featurespace.vtk";
    boost::replace_all(fs_filename, "/", "_");
    INFO << "Writing Featurespace to " + fs_filename << " ... ";
    VisitUtils<TypeParam>::write_featurespace_vtk(fs_filename, this->m_featureSpace);
    if (INFO_ENABLED) cout << "done." << endl;
#endif

    // G'is a kernel
    GaussianNormalKernel<TypeParam> kernel(1.0f);
    cout << setiosflags(ios::fixed) << setprecision(TEST_PRINT_PRECISION);

    for (size_t i = 0; i < this->m_origins.size(); i++) {
        typename FeatureSpace<TypeParam>::Trajectory *trajectory = NULL;

        INFO << "Entering iteration at " << this->m_origins[i].values << " with bandwidths " << this->m_bandwidths[i] << " ... " << endl;

        IterationOperation<TypeParam> op(this->m_featureSpace, this->m_featureSpaceIndex);

        detection_params_t<TypeParam> params;
        detection_context_t<TypeParam> ctx;
        
        Detection<TypeParam>::initialiseContext(ctx);
        ctx.fs = this->m_featureSpace;
        ctx.coord_system = this->m_coordinate_system;
        ctx.data_store = this->m_data_store;

        WeightFunction<TypeParam> *weight 
                = new DefaultWeightFunction<TypeParam>(params,ctx);

        RangeSearchParams<TypeParam> *searchParams = new RangeSearchParams<TypeParam>(this->m_bandwidths[i]);

        trajectory = op.get_trajectory(&this->m_origins[i],
                searchParams,
                &kernel,
                weight,
                TERMCRIT_EPSILON,
                TERMCRIT_ITER);

        vector<TypeParam> last = trajectory->back();

        if (INFO_ENABLED)
            cout << "done. (" << trajectory->size() << " points, ending at " << last << endl;

        // Results should lie around 0 with tolerance termcrit_epsilon

        vector<TypeParam> x2 = trajectory->at(trajectory->size() - 1);
        vector<TypeParam> x1 = trajectory->at(trajectory->size() - 2);
        EXPECT_NEAR(vector_norm(x2 - x1), 0, TERMCRIT_EPSILON);

        // There should be no more than termcrit_iter points in the 
        // trajectory
        EXPECT_LT(trajectory->size(), TERMCRIT_ITER);

        delete weight;

#if WRITE_TRAJECTORIES

        // Write trajectory out

        string fn = "trajectory_" + boost::lexical_cast<string>(i) + ".vtk";

        VisitUtils<TypeParam>::write_pointlist_vtk(fn, trajectory, "trajectory");
#endif

#if WRITE_BANDWIDTH

        typename FeatureSpace<TypeParam>::TrajectoryIterator ti;

        size_t iter_point = 0;

        for (ti = trajectory->begin(); ti != trajectory->end(); ti++) {
            vector<TypeParam> p = *ti;

            string bw_fn = "bandwidth_" + boost::lexical_cast<string>(iter_point++) + ".3D";

            VisitUtils<TypeParam>::write_ellipsis_3d(bw_fn.c_str(), this->m_bandwidths[i], 250, &p);
        }
#endif   

        // clean up

        delete trajectory;
    }
}
#endif

#endif

