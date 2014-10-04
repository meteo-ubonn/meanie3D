#ifndef cf_algorithms_kernel_testcases_h
#define cf_algorithms_kernel_testcases_h

//
//  kernel_testcases_h
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#include <meanie3D/meanie3D.h>

#include <gtest/gtest.h>
#include <limits>

using namespace std;
using namespace m3D;
using testing::Types;
using testing::Test;

/** Fixture for testing circular patterns. Circular patterns in feature space 
 * around the origin of the coordinate system must result in a near null sample 
 * meanshift.
 */
template <class T> 
class DoubleKernelTest : public Test
{
public:
    
    T *kernel;
    
    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp()
    {
        kernel = new T(1.0);
    }
    
    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown()
    {
        delete kernel;
    }
};

template <class T> 
class FloatKernelTest : public Test
{
public:
    
    T *kernel;
    
    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp()
    {
        kernel = new T(1.0);
    }
    
    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown()
    {
        delete kernel;
    }
};


typedef Types< GaussianNormalKernel<double>,
               EpanechnikovKernel<double>,
               UniformKernel<double> > DoubleKernels;

TYPED_TEST_CASE( DoubleKernelTest, DoubleKernels );

// Then use TYPED_TEST(TestCaseName, TestName) to define a typed test,
// similar to TEST_F.
TYPED_TEST( DoubleKernelTest, TestDouble ) 
{
    // Perform within limits
    
    // 1-D (profile)
    
    EXPECT_GE( this->kernel->apply( numeric_limits<double>::max()), 0.0 );
    
    EXPECT_GE( this->kernel->apply( numeric_limits<double>::min()), 0.0 );

    // N-D
    
    for ( size_t dim = 2; dim <= 5; dim++ )
    {
        vector<double> maxVector(dim,numeric_limits<double>::max());
        
        EXPECT_GE( this->kernel->apply( maxVector ), 0.0 );
        
        
        vector<double> minVector(dim,numeric_limits<double>::min());
        
        EXPECT_GE( this->kernel->apply( minVector ), 0.0 );
    }
    
    // Positive definit and monotonically decreasing
    
    double previousValue = this->kernel->apply(0.0);
    
    for ( int i=1; i < 1000; i++ )
    {
        double value = this->kernel->apply((double) i);
        
        EXPECT_LE( value, previousValue );
        
        previousValue = value;
    }
}

typedef Types< GaussianNormalKernel<float>, 
               EpanechnikovKernel<float>, 
               UniformKernel<float> > FloatKernels;

TYPED_TEST_CASE( FloatKernelTest, FloatKernels );

// Then use TYPED_TEST(TestCaseName, TestName) to define a typed test,
// similar to TEST_F.
TYPED_TEST( FloatKernelTest, TestFloat ) 
{
    // Perform within limits
    
    // 1-D (profile)
    
    EXPECT_GE( this->kernel->apply( numeric_limits<float>::max()), 0.0 );
    
    EXPECT_GE( this->kernel->apply( numeric_limits<float>::min()), 0.0 );
    
    // N-D
    
    for ( size_t dim = 2; dim <= 5; dim++ )
    {
        vector<float> maxVector(dim,numeric_limits<float>::max());
        
        EXPECT_GE( this->kernel->apply( maxVector ), 0.0 );
        
        
        vector<float> minVector(dim,numeric_limits<float>::min());
        
        EXPECT_GE( this->kernel->apply( minVector ), 0.0 );
    }
    
    // Positive definit and monotonically decreasing
    
    float previousValue = this->kernel->apply(0.0);
    
    for ( int i=1; i < 1000; i++ )
    {
        float value = this->kernel->apply((float) i);
        
        EXPECT_LE( value, previousValue );
        
        previousValue = value;
    }
}





#endif
