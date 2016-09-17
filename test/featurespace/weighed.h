#ifndef M3D_TEST_FS_WEIGHED_H
#define M3D_TEST_FS_WEIGHED_H

//
//  weighed.h
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#include "../testcase_base.h"

#pragma mark -
#pragma mark Test Fixture

template<class T>
class FSWeighedTest2D : public FSTestBase<T>
{
protected:

    //
    // Protected member variables
    //

    // The half-axis numbers for the ellipsoids 
    vector<vector<T> > m_bandwidths;

    GaussianNormal<T> m_distribution;

    //
    // Protected methods
    //

    void create_uniform_distribution_recursive(const NcVar &var,
                                               size_t modulo,
                                               size_t dimensionIndex,
                                               typename CoordinateSystem<T>::GridPoint &gridpoint);

    void create_uniform_distribution(const NcVar &var,
                                     size_t modulo);

public:

    FSWeighedTest2D();

    virtual void SetUp();

    virtual void TearDown();

};

template<class T>
class FSWeighedTest3D : public FSWeighedTest2D<T>
{
public:
    FSWeighedTest3D();
};

#include "weighed_impl.h"

#endif
