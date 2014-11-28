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

#ifndef M3D_SET_UTILS_H
#define M3D_SET_UTILS_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>
#include <math.h>
#include <cmath>
#include <assert.h>
#include <set>

namespace std {

    // Convenience operator << for printing vectors out to streams

    template < class T >
    std::ostream& operator << (std::ostream& os, const set<T>& v)
    {
        size_t count = 0;
        os << "{";
        for (typename set<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
        {
            os << *ii;

            if ( count < (v.size()-1))
            {
                os << ",";
            }
            count++;
        }
        os << "}";
        return os;
    }
}

namespace m3D { namespace utils { namespace sets {

    template <typename T>
    std::string
    to_string(const std::set<T> &v)
    {
        stringstream str( stringstream::in | stringstream::out );

        str << "{";

        typename std::set<T>::const_iterator ci;

        size_t i=0;

        for (ci = v.begin(); ci != v.end(); ci++ )
        {
            str << (*ci);

            if ( i < (v.size()-1) )
            {
                str << ",";
            }

            i++;
        };

        str << "}";

        return str.str();
    }

    template <typename T>
    std::set<T>
    from_string(const std::string &const_str)
    {
        std::set<T> result;

        string str = const_str;

        // Decode mode

        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

        boost::char_separator<char> sep(",");

        char chars[] = "{}";

        for (unsigned int i = 0; i < strlen(chars); ++i)
        {
            // you need include <algorithm> to use general algorithms like std::remove()
            str.erase (std::remove(str.begin(), str.end(), chars[i]), str.end());
        }

        tokenizer tokens(str, sep);

        vector<T> mode;

        for ( tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter )
        {
            string token = *tok_iter;

            result.insert( boost::lexical_cast<T>(token) );
        }

        return result;
    }
}}}

#endif
