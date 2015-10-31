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


#ifndef M3D_ID
#define M3D_ID

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <map>
#include <set>
#include <vector>

extern "C" {
    namespace m3D {
        
        /**
         * Unique cluster id (in one tracking run)
         */
        typedef unsigned long uuid_t;

        /** Data type for object identifiers
         */
        typedef unsigned long id_t;

        typedef std::set<id_t> id_set_t;

        typedef std::vector<id_t> id_vec_t;

        typedef std::map<id_t, id_set_t> id_map_t;

        // Constants

        static const id_t NO_ID = std::numeric_limits<id_t>::max();

        static const id_t MIN_ID = 0;

        static const id_t MAX_ID = std::numeric_limits<id_t>::max() - 1;

        // Methods

        /** Increment the ID, rotating around if necessary.
         * @param current id, incremented after the call
         * @return next id
         */
        m3D::id_t next_id(m3D::id_t &current)
        {
            if (current >= m3D::MAX_ID) {
                current = m3D::MIN_ID;
            } else {
                current++;
            }

            return current;
        }
    }
}

#endif
