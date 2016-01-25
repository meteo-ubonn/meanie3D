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

#ifndef M3D_WEIGHTFUNCTION_H
#define M3D_WEIGHTFUNCTION_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <meanie3D/featurespace/point.h>

namespace m3D {

    /** Weight function interface. The weight function plays an important
     * role in the mean-shift algorithm. It replaces density as criterion
     * in places, where the actual spatial density is homogenous, as is
     * often the case in scientific, gridded data sets.
     */
    template <typename T>
    class WeightFunction
    {
    public:

        /** Weight at given point in feature-space 
         * @param point object
         * @return weight
         */
        virtual T operator()(const typename Point<T>::ptr p) const = 0;

        virtual ~WeightFunction()
        {
        }
    };
}

#endif
