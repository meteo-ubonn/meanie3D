/* 
 * File:   test_linear_mapping.h
 * Author: simon
 *
 * Created on September 17, 2014, 12:55 PM
 */

#ifndef M3D_TEST_LINEAR_MAPPING_H
#define    M3D_TEST_LINEAR_MAPPING_H

#include <meanie3D/array/linear_index_mapping.h>

#include <gtest/gtest.h>
#include <vector>
#include <set>

using namespace testing;
using namespace m3D;

class LinearMappingTest : public testing::Test
{
public:
};

TEST(LinearMappingTest, LinearMappingTest) {
    // 1D .. 5D

    for (size_t N = 1; N < 5; N++) {
        // 4x4x....
        std::vector<size_t> dim_sizes(N);

        for (size_t i = 0; i < N; i++) {
            dim_sizes[i] = 5 * (i + 1);
        }

        LinearIndexMapping mapping(dim_sizes);

        size_t size = mapping.size();

        std::set<std::vector<int> > points;

        for (size_t i = 0; i < size; i++) {
            vector<int> p = mapping.linear_to_grid(i);

            for (size_t j = 0; j < N; j++) {
                ASSERT_LT(p[j], dim_sizes[j]);
            }

            points.insert(p);
        }

        // there must be 'size' numbers of points
        // in the set. 

        EXPECT_EQ(size, points.size());

        // Run the process in reverse and make
        // sure every point is found in the set

        for (size_t linear_index = 0; linear_index < mapping.size(); linear_index++) {
            vector<int> point = mapping.linear_to_grid(linear_index);

            std::set<vector < int> > ::iterator
            fi;

            fi = points.find(point);

            EXPECT_NE(fi, points.end());
        }
    }
}

#endif	/* TEST_LINEAR_MAPPING_H */

