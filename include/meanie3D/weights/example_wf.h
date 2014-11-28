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

/* 
 * File:   example_wf.h
 * Author: simon
 *
 * Created on November 6, 2014, 11:29 AM
 */

#ifndef M3D_EXAMPLE_WF_H
#define M3D_EXAMPLE_WF_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils.h>
#include <meanie3D/weights/weight_function.h>
#include <vector>

namespace m3D {
    
    template <class T>
    class ExampleWF : public WeightFunction<T>
    {
    private:

        MultiArray<T>       *m_weight;

    public:

        ExampleWF(FeatureSpace<T> *fs, const vector<T> &center)
        {
            using namespace utils::vectors;
            
            m_weight = new MultiArrayBlitz<T>(fs->coordinate_system->get_dimension_sizes(),0.0);
            
            for (size_t i=0; i < fs->points.size(); i++)
            {
                Point<T> *p = fs->points[i];
                T distance = vector_norm(center - p->coordinate);
                m_weight->set(p->gridpoint, distance);
            }
        }

        ~ExampleWF()
        {
            if (m_weight!=NULL)
            {
                delete m_weight;
                m_weight = NULL;
            }
        }

        T operator()(const typename Point<T>::ptr p) const
        {
            return m_weight->get(p->gridpoint);
        }
    };
}


#endif	/* EXAMPLE_WF_H */

