#ifndef M3D_TEST_FS_CIRCULARPATTERN_H
#define M3D_TEST_FS_CIRCULARPATTERN_H

//
//  circular_pattern.h
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

// #include <cf-algorithms/cf-algorithms.h>

#include "testcase_base.h"

#pragma mark -
#pragma mark Constants

/** 'Sharpness' of the test data ellipsoids 
 */
static const float ELLIPSE_FUZZINESS = 0.05f;

#pragma mark -
#pragma mark Test Fixture

template <class T> 
class FSCircularPatternTest2D : public FSTestBase<T>
{
    
protected:
    
    //
    // Protected member variables
    //

    // The half-axis numbers for the ellipsoids 
    vector< vector<T> > m_bandwidths;
    
    // Calculates from the number of points created.
    vector<size_t> m_nearest_neighbours;
    
    //
    // Protected methods
    //
    
    bool isPointOnEllipse( vector<T> coordinate, vector<T> axis );
    
    void create_ellipsoid_recursive( NcVar &var, vector<T> &h, size_t dimensionIndex, typename CoordinateSystem<T>::GridPoint &gridpoint );
    
    void create_ellipsoid( NcVar &var, vector<T> h );

public:
    
    FSCircularPatternTest2D();
    
    virtual void SetUp();
    
    virtual void TearDown();
    
};

template <class T> 
class FSCircularPatternTest3D : public FSCircularPatternTest2D<T>
{
public:
    FSCircularPatternTest3D();
};

#include "circular_pattern_impl.h"

#endif
