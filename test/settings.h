#ifndef M3D_TEST__FS_SETTINGS_H
#define M3D_TEST__FS_SETTINGS_H

#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

#pragma mark -
#pragma mark Constants

/** The value of the generated variables. 
 */
static const float FS_VALUE_MAX = 3.0;

/** Termination criterium for meanshift iterations
 */
static const float TERMCRIT_EPSILON = 0.01;

/** Termination criterium for meanshift iterations
 */
static const size_t TERMCRIT_ITER = 1000;

/** Number of gridpoints in each dimension. Cranking this up will increase
 * calculation time for the tests significantly!
 */
static const size_t NUMBER_OF_GRIDPOINTS = 100;

#pragma mark -
#pragma mark Declaration

/** Class for bundling settings for the featurespace / clustering tests 
 */
class FSTestSettings
{
protected:

    size_t m_num_dimensions;

    size_t m_num_vars;

    size_t m_num_gridpoints;

    vector<float> m_axis_bound_values;

    std::string m_test_filename;

public:

    /** @constructor
     * @param Number of dimensions
     * @param Number of variables
     * @param Number of gridpoints on each axis
     * @param Values on each axis go from -axis_bound_value to +axis_bound_value in each dimension
     * @param Filename of the generated test data
     */
    FSTestSettings(size_t num_dims,
                   size_t num_vars,
                   vector<float> axis_bound_value,
                   size_t num_gridpoints = NUMBER_OF_GRIDPOINTS,
                   std::string test_filename = "fs_test.nc");

    /** @constructor
     * @param Number of dimensions
     * @param Number of variables
     * @param Number of gridpoints on each axis
     * @param Filename of the generated test data
     */
    FSTestSettings(size_t num_dims,
                   size_t num_vars,
                   size_t num_gridpoints = NUMBER_OF_GRIDPOINTS,
                   std::string test_filename = "fs_test.nc");

    /** Number of dimensions and dimension variables in the testfile (max 4) 
     */
    size_t num_dimensions();

    /** Number of non-dimension variables in the testfile 
     */
    size_t num_vars();

    /** Number of grid points in each direction 
     */
    size_t num_gridpoints();

    /** Convenience accessor. 
     * @return num_dimensions() + num_vars()
     */
    size_t fs_dim();

    /** Numerical bound on the axis of each dimension 
     */
    vector<float> &axis_bound_values();

    /** Setter 
     */
    void set_axis_bound_values(vector<float> &newValue);

    /** Name of the file used for the test data 
     */
    std::string test_filename();

    void set_test_filename(std::string value);
};

#endif