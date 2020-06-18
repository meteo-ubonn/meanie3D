#ifndef M3D_TEST_FSCIRCULARPATTERN_IMPL_H
#define M3D_TEST_FSCIRCULARPATTERN_IMPL_H

using namespace m3D;
using namespace m3D::utils::vectors;

#include "../testcase_base.h"

#pragma mark -
#pragma mark Test case 1 - ellipsoid sample means - data generation

/** This function checks, if a coordinate is on the surface of the ellipsoid
 * described by the given set of axis
 * @param coordinate
 * @param axis
 * @return yes or no
 */
template<class T>
bool FSCircularPatternTest2D<T>::isPointOnEllipse(vector<T> coordinate, vector<T> axis) {
    assert(axis.size() >= coordinate.size());
    float value = 0;
    for (size_t index = 0; index < coordinate.size(); index++) {
        value += ((coordinate[index] * coordinate[index]) / (axis[index] * axis[index]));
    }
    return (fabs(value - 1.0) <= ELLIPSE_FUZZINESS);
}

template<class T>
void FSCircularPatternTest2D<T>::create_ellipsoid_recursive(NcVar &var,
                                                            vector<T> &h,
                                                            size_t dimensionIndex,
                                                            vector<int> &gridpoint) {
    NcDim dim = var.getDim(dimensionIndex);
    if (dimensionIndex < (this->coordinate_system()->rank() - 1)) {
        for (int index = 0; index < dim.getSize(); index++) {
            gridpoint[dimensionIndex] = index;
            create_ellipsoid_recursive(var, h, dimensionIndex + 1, gridpoint);
        }
    } else {
        for (int index = 0; index < dim.getSize(); index++) {
            gridpoint[dimensionIndex] = index;

            // get the variables together and construct the cartesian coordinate
            // of the current point. If it's on the ellipse, put it in the variable
            typename CoordinateSystem<T>::Coordinate coordinate = this->coordinate_system()->newCoordinate();
            this->coordinate_system()->lookup(gridpoint, coordinate);
            vector<size_t> gp(gridpoint.begin(), gridpoint.end());
            if (isPointOnEllipse(coordinate, h)) {
                T value = (T) FS_VALUE_MAX;
                var.putVar(gp, value);
                this->m_pointCount = this->m_pointCount + 1;
            }
        }
    }
}

template<class T>
void FSCircularPatternTest2D<T>::create_ellipsoid(NcVar &var, vector<T> h) {
    vector<int> gridpoint(this->coordinate_system()->rank(), 0);
    this->m_pointCount = 0;
    create_ellipsoid_recursive(var, h, 0, gridpoint);
    this->m_totalPointCount += this->m_pointCount;
}

template<class T>
void FSCircularPatternTest2D<T>::SetUp() {
    FSTestBase<T>::SetUp();

    // Set the bandwidths
    size_t bw_size = this->m_settings->num_dimensions();
    this->m_bandwidths.push_back(vector<T>(bw_size, 1.0));
    this->m_bandwidths.push_back(vector<T>(bw_size, 2.0));
    this->m_bandwidths.push_back(vector<T>(bw_size, 3.0));
    this->m_bandwidths.push_back(vector<T>(bw_size, 4.0));

    // Figure out the bandwidth maximum
    vector<float> max_h(this->m_settings->fs_dim(), 0.0);
    for (size_t i = 0; i < m_bandwidths.size(); i++) {
        vector<T> h = m_bandwidths.at(i);
        assert(h.size() == bw_size);
        for (size_t j = 0; j < bw_size; j++) {
            assert(h[j] > 0);
            if (h[j] > max_h[j]) {
                max_h[j] = h[j];
            }
        }
    }

    // Make the actual values 1.1 times bigger, so that in graphical
    // representations there is some boundary
    max_h *= 1.1f;
    this->m_settings->set_axis_bound_values(max_h);

    // Generate dimensions and dimension variables according to
    // the current settings
    this->generate_dimensions();

    // Create a variable
    // test case 1 : ellipsis for unweighed sample mean
    NcVar var = this->add_variable("circular_pattern_test", 0.0, FS_VALUE_MAX);

    // Create an ellipse with each bandwidth set in the new variable
    typename vector<vector<T> >::iterator it;
    for (it = m_bandwidths.begin(); it != m_bandwidths.end(); it++) {
        INFO << "Creating ellipsoid at origin with ranges " << (*it) << " ... ";
        create_ellipsoid(var, *it);
        if (INFO_ENABLED)
            cout << " (" << this->m_pointCount << " points)." << endl;
    }

    // Calculate KNN
    for (size_t i = 0; i < m_bandwidths.size(); i++) {
        this->m_nearest_neighbours.push_back((i + 1) * (this->m_totalPointCount / m_bandwidths.size()));
    }

    this->generate_featurespace();
}

template<class T>
void FSCircularPatternTest2D<T>::TearDown() {
    FSTestBase<T>::TearDown();
}

#pragma mark -
#pragma mark Test parameterization

template<class T>
FSCircularPatternTest2D<T>::FSCircularPatternTest2D() : FSTestBase<T>() {
    INFO << "Setting up 2D test with typeid " << typeid(T).name() << endl;
    std::string filename = FSTestBase<T>::filename_from_current_testcase();
    INFO << "Test filename = " << filename << endl;
    this->m_settings = new FSTestSettings(2, 1, NUMBER_OF_GRIDPOINTS, FSTestBase<T>::filename_from_current_testcase());
}

template<class T>
FSCircularPatternTest3D<T>::FSCircularPatternTest3D() : FSCircularPatternTest2D<T>() {
    INFO << "Setting up 3D test with typeid " << typeid(T).name() << endl;
    std::string filename = FSTestBase<T>::filename_from_current_testcase();
    INFO << "Test filename = " << filename << endl;
    this->m_settings = new FSTestSettings(3, 1, NUMBER_OF_GRIDPOINTS, FSTestBase<T>::filename_from_current_testcase());
}

// 2D
#if RUN_2D

TYPED_TEST_CASE(FSCircularPatternTest2D, DataTypes);

TYPED_TEST(FSCircularPatternTest2D, FS_CicrularPattern_2D_Test)
{
    // origin
    vector<TypeParam> x(this->m_settings->fs_dim(), 0);
    x[ x.size() - 1 ] = FS_VALUE_MAX;

    // G'is a kernel
    GaussianNormalKernel<TypeParam> kernel(1.0f);
    cout << setiosflags(ios::fixed) << setprecision(TEST_PRINT_PRECISION);

    // iterate over the bandwidths
    for (size_t i = 0; i < this->m_bandwidths.size(); i++) {
        vector<TypeParam> h = this->m_bandwidths.at(i);
        h.push_back(FS_VALUE_MAX);
        MeanshiftOperation<TypeParam> op(this->m_featureSpace, this->m_featureSpaceIndex);
        RangeSearchParams<TypeParam> *params = new RangeSearchParams<TypeParam>(((TypeParam) 1.1) * h);
        vector<TypeParam> m = op.meanshift(x, params, &kernel);
        EXPECT_NEAR(vector_norm(m), 0.0, this->coordinate_system()->resolution_norm());
    }
}
#endif

// 3D
#if RUN_3D

TYPED_TEST_CASE(FSCircularPatternTest3D, DataTypes);

TYPED_TEST(FSCircularPatternTest3D, FS_CicrularPattern_3D_Test)
{
    // origin
    vector<TypeParam> x(this->m_settings->fs_dim(), 0);
    x[ x.size() - 1 ] = FS_VALUE_MAX;

    // G'is a kernel
    Kernel<TypeParam> *kernel = new GaussianNormalKernel<TypeParam>(1.0f);
    cout << setiosflags(ios::fixed) << setprecision(TEST_PRINT_PRECISION);

    // iterate over the bandwidths
    for (size_t i = 0; i < this->m_bandwidths.size(); i++) {
        vector<TypeParam> h = this->m_bandwidths.at(i);
        h.push_back(FS_VALUE_MAX);
        MeanshiftOperation<TypeParam> op(this->m_featureSpace, this->m_featureSpaceIndex);
        RangeSearchParams<TypeParam> *params = new RangeSearchParams<TypeParam>(h);
        vector<TypeParam> m = op.meanshift(x, params, kernel);
        EXPECT_NEAR(vector_norm(m), 0.0, this->coordinate_system()->resolution_norm());
    }
    delete kernel;
}

#endif

#endif