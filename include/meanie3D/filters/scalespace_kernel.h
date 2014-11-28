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

#ifndef M3D_SCALESPACEKERNEL_H
#define M3D_SCALESPACEKERNEL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

namespace m3D {

        using namespace std;

    /** Kernel profile g for scale-space filtering. Allows for 
     * pre-rendering at a number of points for speed-up.
     */
    template <class T>
    class ScaleSpaceKernel
    {

    private:

        T          m_t;
        T          m_gauging_factor;
        vector<T>  m_values;

    public:

#pragma mark -
#pragma mark Constructor/Destructor

        /** Constructor
         * @param scale parameter t
         */
        ScaleSpaceKernel(const T &t);

        /** Constructor. Pre-calculates the values at the given
         * distances. 
         *
         * @param scale parameter t
         * @param distances
         */
        ScaleSpaceKernel(const T &t, const vector<T> &distances);

        /** Copy constructor
         * @param other kernel
         */
        ScaleSpaceKernel(const ScaleSpaceKernel<T> &o);

        /** Destructor
         */
        virtual ~ScaleSpaceKernel();

        /** Copy operator 
         */
        ScaleSpaceKernel<T>
        operator = ( const ScaleSpaceKernel<T>& other );


#pragma mark -
#pragma mark Accessors

        /** @return pre-calculated values.
         */
        const vector<T> &values();

        /** @return true if the kernel was constructed with pre-sampled values, false else.
         */
        bool isPreSampled();

        const T scale() {return m_t;};

#pragma mark -
#pragma mark Calculation

        /** Obtains the pre-calculated value at index i
         * @param index
         * @return value
         * @throws
         */
        const T value(size_t index);

        /** Calculates the value at distance r
         * @param distance r
         * @return calculated value;
         */
        const T value(T r);
    };
}

#endif
