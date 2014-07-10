#ifndef _M3D_DEFINES_H_
#define _M3D_DEFINES_H_

// Threading

#define WITH_TBB 1
#define WITH_BOOST_THREADS 0

#define PROVIDE_THREADSAFETY 0

#if WITH_BOOST_THREADS
	#define PROVIDE_MUTEX 1
#endif

#if WITH_TBB
	#define PROVIDE_MUTEX 1
#endif

// Featurespace debugging

#define DEBUG_MEANSHIFT_GRAPH 0
#define DEBUG_CLUSTER_MERGING 0
#define DEBUG_CLUSTER_MERGING_DECISION 0
#define DEBUG_TRACKING 0
#define DEBUG_HISTOGRAM_CORRELATION 0

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

// Enable the replacement of zero-shift vectors
// with an average of neighboring points
#define REPLACE_ZEROSHIFT_VECTORS 0

// Figure out the strongest neighbor to the
// zero-shift cluster and add it to the cluster
#define ADD_STRONGEST_NEIGHBOUR 0

#endif
