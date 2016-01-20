#ifndef M3D_TEST_FS_CLUSTERING_IMPL_H
#define M3D_TEST_FS_CLUSTERING_IMPL_H

#include <meanie3D/meanie3D.h>

#include "variable_weighed_impl.h"

template<class T>
void FSClusteringTest2D<T>::write_cloud(const NcVar &var, vector<T> mean, vector<T> deviation)
{
    using namespace netCDF;
    using namespace m3D;
    using namespace m3D::utils::vectors;

    // start generating random points (with random values between 0 and 1)
    size_t numPoints = 0;

    // Use a simple gaussian normal function to simulate a value distribution
    GaussianNormal<T> gauss;
    T gauss_zero = gauss(vector<T>(this->coordinate_system()->rank(), 0));

    typename CoordinateSystem<T>::GridPoint gridpoint = this->coordinate_system()->newGridPoint();
    do {
        // generate a random coordinate
        for (size_t d = 0; d < this->coordinate_system()->rank(); d++) {
            bool valid = false;
            while (!valid) {

                NcDim dim = this->coordinate_system()->dimensions()[d];
                NcVar dimVar = this->coordinate_system()->dimension_variable(dim);

                T min, max;
                dimVar.getAtt("valid_min").getValues(&min);
                dimVar.getAtt("valid_max").getValues(&max);
                float rand = box_muller(mean[d], deviation[d]);

                // re-transform to grid coordinates
                size_t n = (size_t) round((dim.getSize() - 1) * (rand - min) / (max - min));
                if (n < dim.getSize()) {
                    gridpoint[d] = n;
                    valid = true;
                }
            }
        }

        // value goes to
        vector<T> x = this->coordinate_system()->newCoordinate();
        this->coordinate_system()->lookup(gridpoint, x);
        T value = (T) (FS_VALUE_MAX * gauss(x - mean)) / gauss_zero;

        // cout << "x=" << x << " mean=" << mean << " g(x-m)=" << value << endl;
        vector<size_t> gp(gridpoint.begin(), gridpoint.end());
        var.putVar(gp, value);
        numPoints++;

        this->m_pointCount++;

    } while (numPoints < m_cloudSize);
}

template<class T>
void FSClusteringTest2D<T>::create_clouds_recursive(const NcVar &var, size_t dimensionIndex, typename CoordinateSystem<T>::GridPoint &gridpoint)
{
    using namespace netCDF;

    NcDim dim = var.getDim(dimensionIndex);

    size_t increment = this->m_division_increments[dim];

    if (dimensionIndex < (var.getDimCount() - 1)) {
        for (int index = 1; index < m_divisions; index++) {
            gridpoint[dimensionIndex] = index * increment;

            create_clouds_recursive(var, dimensionIndex + 1, gridpoint);
        }
    } else {
        for (int index = 1; index < m_divisions; index++) {
            gridpoint[dimensionIndex] = index * increment;

            // get the variables together and construct the cartesian coordinate
            // of the current point. If it's on the ellipse, put it in the variable

            vector<T> coordinate(var.getDimCount());

            this->coordinate_system()->lookup(gridpoint, coordinate);

            cout << "\tWriting cloud at gridpoint=" << gridpoint << " coordinate=" << coordinate << " deviation=" << m_deviation << endl;

            write_cloud(var, coordinate, m_deviation);
        }
    }
}

template<class T>
void FSClusteringTest2D<T>::create_clouds(const NcVar &var)
{
    // calculate the divisions in terms of grid points
    // calculate the deviations
    vector<NcDim *>::iterator dim_iter;
    for (size_t index = 0; index < this->coordinate_system()->rank(); index++) {
        NcDim dim = this->coordinate_system()->dimensions()[index];
        NcVar dim_var = this->coordinate_system()->dimension_variable(dim);
        size_t number_gridpoints = utils::netcdf::num_vals(dim_var) / m_divisions;
        m_division_increments[ dim ] = number_gridpoints;
        m_deviation.push_back(0.4 / m_divisions);
    }
    typename CoordinateSystem<T>::GridPoint gridpoint = this->coordinate_system()->newGridPoint();
    this->m_pointCount = 0;
    cout << "Creating clusters at the intersection of " << m_divisions << " lines per axis ..." << endl;
    create_clouds_recursive(var, 0, gridpoint);
    cout << "done. (" << this->m_pointCount << " points)" << endl;
    this->m_totalPointCount += this->m_pointCount;
}

template<class T>
void FSClusteringTest2D<T>::SetUp()
{
    FSTestBase<T>::SetUp();
    m_smoothing_scale = 0.01;
    // Set the bandwidths
    size_t bw_size = this->m_settings->fs_dim();
    for (size_t i = 0; i < bw_size; i++) {
        m_fuzziness.push_back(1.0 / m_divisions);
    }
    this->m_bandwidths.push_back(vector<T>(bw_size, 1.0 / m_divisions));

    // Generate dimensions and dimension variables according to
    // the current settings
    this->generate_dimensions();

    // Create a variable
    // test case 1 : ellipsis for unweighed sample mean
    NcVar var = this->add_variable("cluster_test", 0.0, FS_VALUE_MAX);
    create_clouds(var);
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
    this->m_settings = new FSTestSettings(2, 1, NUMBER_OF_GRIDPOINTS, FSTestBase<T>::filename_from_current_testcase());
    this->m_divisions = 3;
    this->m_cloudSize = 50 * NUMBER_OF_GRIDPOINTS / m_divisions;
    this->m_smoothing_scale = 0.01;
}

template<class T>
FSClusteringTest3D<T>::FSClusteringTest3D() : FSClusteringTest2D<T>()
{
    this->m_settings = new FSTestSettings(3, 1, NUMBER_OF_GRIDPOINTS, FSTestBase<T>::filename_from_current_testcase());
    this->m_divisions = 4;
    this->m_cloudSize = 50 * NUMBER_OF_GRIDPOINTS / this->m_divisions;
    this->m_smoothing_scale = 0.01;
}

// 2D
#if RUN_2D

TYPED_TEST_CASE(FSClusteringTest2D, DataTypes);

TYPED_TEST(FSClusteringTest2D, FS_Clustering_2D_Range_Test)
{
    using namespace m3D;
    using namespace m3D::utils;
    using namespace m3D::utils::vectors;

    GaussianNormal<TypeParam> gauss;
    // TypeParam gauss_zero = gauss( vector<TypeParam>(this->coordinate_system()->size(),0) );

    cout << setiosflags(ios::fixed) << setprecision(TEST_PRINT_PRECISION);

    for (size_t i = 0; i < this->m_bandwidths.size(); i++) {
        
        detection_params_t<TypeParam> params = Detection<TypeParam>::defaultParams();
        params.ranges = this->m_bandwidths.at(i);
        params.kernel_name = "gauss";
        params.variables = this->m_variables;
        params.dimensions = this->m_dimensions;
        params.dimension_variables = this->m_dimension_variables;
        params.min_cluster_size = 20;
        
        start_timer();
        detection_context_t<TypeParam> ctx;
        Detection<TypeParam>::run(params,ctx);
        double time = stop_timer();
        cout << "done. (" << time << " seconds, " 
             << ctx.clusters->size() << " modes found:" << endl;

        EXPECT_EQ(ctx.clusters->clusters.size(), 4);

        ctx.clusters->print();
        Detection<TypeParam>::cleanup(params,ctx);
        ClusterList<TypeParam>::reset_clustering(this->m_featureSpace);
    }
}
#endif

// 3D
#if RUN_3D

TYPED_TEST_CASE(FSClusteringTest3D, DataTypes);

TYPED_TEST(FSClusteringTest3D, FS_Clustering_3D_Test)
{
    using namespace m3D;
    using namespace m3D::utils;
    using namespace m3D::utils::vectors;

    const ::testing::TestInfo * const test_info = ::testing::UnitTest::GetInstance()->current_test_info();

    cout << setiosflags(ios::fixed) << setprecision(TEST_PRINT_PRECISION);

    for (size_t i = 0; i < this->m_bandwidths.size(); i++) {
        
        detection_params_t<TypeParam> params = Detection<TypeParam>::defaultParams();
        
        params.ranges = this->m_bandwidths.at(i);
        params.ranges *= ((TypeParam) 0.25);
        
        params.kernel_name = "gauss";
        params.variables = this->m_variables;
        params.dimensions = this->m_dimensions;
        params.dimension_variables = this->m_dimension_variables;
        params.min_cluster_size = 20;
        
        start_timer();
        detection_context_t<TypeParam> ctx;
        Detection<TypeParam>::run(params,ctx);
        double time = stop_timer();
        cout << "done. (" << time << " seconds, " 
             << ctx.clusters->size() << " modes found:" << endl;

        EXPECT_EQ(ctx.clusters->clusters.size(), 27);

        ctx.clusters->print();
        Detection<TypeParam>::cleanup(params,ctx);
        ClusterList<TypeParam>::reset_clustering(this->m_featureSpace);
    }
}

#endif


#endif

