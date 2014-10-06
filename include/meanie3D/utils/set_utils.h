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
