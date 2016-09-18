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

#ifndef M3D_GAUSSIAN_NORMAL_H
#define M3D_GAUSSIAN_NORMAL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/featurespace.h>

namespace m3D {
    namespace utils {

        using namespace std;

        /** Specialization of weight function as gaussian normal function
         * around a given point. Used for test purposes.
         */
        template<class T>
        class GaussianNormal : public WeightFunction<T>
        {
            static T FS_VALUE_MAX;

        public:

            virtual T operator()(const vector<T> &values) const {
                static float sigma_square = 0.5;

                T norm = vectors::vector_norm(values);

                T value =
                        FS_VALUE_MAX * (1.0 / (sqrt(2 * M_PI * sigma_square))) * exp(-0.5 * norm * norm / sigma_square);

                return value;
            }

            virtual T operator()(const typename Point<T>::ptr p) const {
                return this->operator()(p->values);
            }

            virtual T operator()(const vector<int> &gridpoint) const {
                throw "not implemented";
            }

        };

        template<typename T>
        T GaussianNormal<T>::FS_VALUE_MAX = 3.0;
    }
}

#endif
