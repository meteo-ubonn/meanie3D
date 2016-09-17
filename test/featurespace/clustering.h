#ifndef M3D_TEST_FS_CLUSTERING_H
#define M3D_TEST_FS_CLUSTERING_H

//
//  clustering.h
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#include "../testcase_base.h"

#pragma mark -
#pragma mark Test Fixture

template<class T>
class FSClusteringTest2D : public FSTestBase<T>
{
protected:

    //
    // Protected member variables
    //

    /** Run the clustering with a set of bandwidths
     */
    vector<vector<T> > m_bandwidths;

    /** Run the clustering with a set of bandwidths
     */
    vector<T> m_fuzziness;

    /** Size of the cloud at each point 
     */
    size_t m_cloudSize;

    vector<T> m_deviation;

    /** Number of divisions (same for all axis's )
     */
    size_t m_divisions;

    /** Divider increments in terms of gridpoints 
     */
    map<NcDim, size_t> m_division_increments;

    //
    // Protected methods
    //

    /** Writes a gaussian cloud at the position 'mean' with spread 'deviation'.
     * @param variable
     * @param mean / position
     * @param deviation / spread
     */
    void write_cloud(const NcVar &var, vector<T> mean, vector<T> deviation);

    /** Creates gaussian clouds at the intersection of the lines
     * that result by dividing each axis into m_divisions parts. 
     * (this excludes the intersections at the boundaries)
     * Example1: 2D, 3 divisions => 4 nodes.
     * Example2: 3D, 3 divisions => 9 nodes.
     * Example3: 3D, 4 divisions => 27 nodes
     */
    void create_clouds(const NcVar &var);

    void create_clouds_recursive(const NcVar &, size_t, typename CoordinateSystem<T>::GridPoint &);

public:

    FSClusteringTest2D();

    virtual void SetUp();

    virtual void TearDown();

};

template<class T>
class FSClusteringTest3D : public FSClusteringTest2D<T>
{
public:
    FSClusteringTest3D();
};

#include "clustering_impl.h"

#endif
