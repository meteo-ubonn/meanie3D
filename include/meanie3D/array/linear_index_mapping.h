/* The MIT License (MIT)
 * 
 * (c) Jürgen Simon 2014 (juergen.simon@uni-bonn.de)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#ifndef M3D_LINEAR_INDEX_MAPPING_H
#define	M3D_LINEAR_INDEX_MAPPING_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

namespace m3D { 

    /** Encapsulates the mapping behavior between linear and 
     * grid point indexes. Note that this is not storage related
     * in any way, so there are no strides. This simply provides
     * a bijective projection from grid points to a number in order
     * to facilitate parallelization.
     */
    class LinearIndexMapping
    {
    private:
        std::vector<size_t> m_dimension_sizes;
        std::vector<size_t> m_slice_sizes;
        size_t              m_size;

    public:

        /** Constructs a mapping object 
         */
        LinearIndexMapping(const std::vector<size_t> &dimension_sizes) 
                : m_dimension_sizes(dimension_sizes)
                , m_slice_sizes(std::vector<size_t>(dimension_sizes.size(),1))
                , m_size(1)
        {
            size_t N = m_dimension_sizes.size();

            for (size_t i=0; i<N; i++)
                m_size *= m_dimension_sizes[i];

            // pre-calculate slice sizes. 

            m_slice_sizes[N-1] = 1;

            for (int n=(N-2); n >= 0; n--)
            {
                m_slice_sizes[n] = m_slice_sizes[n+1] * m_dimension_sizes[n+1];
            }

//            cout << "Constructing mapping for " << m_dimension_sizes << endl;
//            cout << "Slice-sizes              " << m_slice_sizes << endl;

        };

        /** @return number of possible points in this mapping. 
         */
        size_t size()
        {
            return m_size;
        }

        /** Maps linear index to grid index components 
         * @param linear index
         * @return grid index
         */
        vector<int> 
        linear_to_grid(size_t linear_index)
        {
            size_t N = m_dimension_sizes.size();

            vector<int> g(N,0);

            size_t position_sum = 0;

            for (int i=0; i < N; i++)
            {
                if (i>0)
                {
                    position_sum += g[i-1] * m_slice_sizes[i-1];
                }

                int idx = (linear_index - position_sum) / m_slice_sizes[i];

                g[i] = idx;
            }

            return g;
        }

        size_t 
        grid_to_index(const vector<int> &g)
        {
            size_t linear_index = 0;

            size_t N = m_dimension_sizes.size();

            for (int i = (N-1); i >=0; i++)
            {
                linear_index += m_slice_sizes[i] * g[i];
            }

            return linear_index;
        }
    };
}

#endif 
