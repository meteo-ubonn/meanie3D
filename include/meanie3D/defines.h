#ifndef M3D_DEFINES_H
#define M3D_DEFINES_H

// If enabled, this switches on vector size assertion
// in vector utils
#define ASSERT_VECTOR_SIZES 0

// ---------------------------------------------------- //
// Configuration
// ---------------------------------------------------- //

// Enable the replacement of zero-shift vectors
// with an average of neighboring points
#define REPLACE_ZEROSHIFT_VECTORS 0

// Figure out the strongest neighbor to the
// zero-shift cluster and add it to the cluster
#define ADD_STRONGEST_NEIGHBOUR 0

// If switched on, the scale-space filter
// excludes all points that were marked as
// 'off limits' in feature-space construction
#define SCALE_SPACE_SKIPS_NON_ORIGINAL_POINTS 0

// Method for rounding vectors to grid resolution

#define GRID_ROUNDING_METHOD_FLOOR 0
#define GRID_ROUNDING_METHOD_ROUND 1
#define GRID_ROUNDING_METHOD_CEIL 0
#define GRID_ROUNDING_METHOD_RINT 0
#define GRID_ROUNDING_METHOD_NONE 0

// ---------------------------------------------------- //
// Debugging Flags(stdout)
// ---------------------------------------------------- //

#define DEBUG_DIMENSION_DATA 0
#define DEBUG_INDEX 0
#define DEBUG_INDEX_SEARCHES 0
#define DEBUG_FEATURESPACE 0
#define DEBUG_ITERATION 0
#define DEBUG_MEANSHIFT_GRAPH 0
#define DEBUG_MEANSHIFT_CALCULATION 0
#define DEBUG_MEANSHIFT_RESULT_CALCULATION 0
#define DEBUG_MEANSHIFT_SAMPLING 0
#define DEBUG_CLUSTER_MERGING 0
#define DEBUG_CLUSTER_MERGING_DECISION 0
#define DEBUG_TRACKING 0
#define DEBUG_HISTOGRAM_CORRELATION 0
#define DEBUG_GRAPH_AGGREGATION 0

// ---------------------------------------------------- //
// Debugging Flags (files)
// ---------------------------------------------------- //

// Write out the constructed index data. Only works for FLANN 
// at this point
#define WRITE_INDEX 0

// Write the feature-space out in .vtk file format
#define WRITE_FEATURESPACE 1

// Write out each mean-shift sample as vtk file
// WARNING: this produces a LOT of data
#define WRITE_MEANSHIFT_SAMPLES 0

// Write out the weights at a set of sample points. If this
// feature is enabled, set the sample points using the property
// FeatureSpace::weight_sample_points
#define WRITE_MEANSHIFT_WEIGHTS 0

// When using the iterative mean-shift, this feature allows
// writing out the individual trajectories as .vtk files.
// WARNING: this produces a LOT of data
#define WRITE_TRAJECTORIES 0

// Write out modes found by the iterations or graph approach
#define WRITE_MODES 0

// Writes out a vtk file that visualizes the selected bandwidth
// parameters
#define WRITE_BANDWIDTH 0

// In the iterative approach, writes out a file that shows iteration
// origins only
#define WRITE_ITERATION_ORIGINS 0

// Write out a VTK file that contains 1 and 0 for
// all points. 1 means, that the point was marked
// as 'off limits' in the original feature-space
// (one of the variables outside of valid_range).
#define WRITE_OFF_LIMITS_MASK 0

// Write out the center of the cluster in a single
// file per cluster
#define WRITE_CLUSTER_CENTERS 0

// Write out the result of the initial clustering
// step into separate files.
#define WRITE_ZEROSHIFT_CLUSTERS 0

// Write out files containing the mean-shift
// vectors for individual clusters alone
#define WRITE_CLUSTER_MEANSHIFT 0

#endif
