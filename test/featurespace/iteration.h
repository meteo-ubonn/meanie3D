#ifndef M3D_TEST_FS_ITERATION_H
#define M3D_TEST_FS_ITERATION_H

//
//  iteration.h
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#include "../testcase_base.h"

#pragma mark -
#pragma mark Test Fixture

template <class T>
class FSIterationTest2D : public FSTestBase<T>
{
protected:

    //
    // Protected member variables
    //

    typename CoordinateSystem<T>::Coordinate *m_center;

    vector<T> *m_mean;

    vector<T> *m_deviation;

    size_t m_cloudSize;

    vector< vector<T> > m_bandwidths;

    vector< Point<T> > m_origins;

    //
    // Protected methods
    //

    void write_cloud(const NcVar &variable);

public:

    FSIterationTest2D();

    virtual void SetUp();

    virtual void TearDown();
};

template <class T>
class FSIterationTest3D : public FSIterationTest2D<T>
{
public:
    FSIterationTest3D();
};

#include "iteration_impl.h"

#endif
