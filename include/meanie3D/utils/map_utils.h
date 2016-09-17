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

#ifndef M3D_MAPUTILS_H
#define M3D_MAPUTILS_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/clustering/id.h>
#include <meanie3D/utils/set_utils.h>

#include <boost/cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>
#include <math.h>
#include <cmath>
#include <assert.h>
#include <vector>
#include <map>

namespace std {

    // Convenience operator << for printing maps out to streams

    template<typename K, typename T>
    std::ostream &operator<<(std::ostream &os, const map<K, T> &v) {
        //        size_t count = 0;
        //        os << "(";
        //        for (typename vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
        //        {
        //            os << *ii;
        //            
        //            if ( count < (v.size()-1))
        //            {
        //                os << ",";
        //            }
        //            count++;
        //        }
        //        os << ")";

        return os;
    }
}

namespace m3D {
    namespace utils {
        namespace maps {

            std::string
            id_map_to_string(const id_map_t &m) {
                using namespace ::m3D::utils::sets;
                stringstream str(stringstream::in | stringstream::out);

                str << "[";

                id_map_t::const_iterator mi;

                for (mi = m.begin(); mi != m.end(); mi++) {
                    str << mi->first;
                    str << ":";
                    str << mi->second;
                    str << ";";
                }

                str << "]";

                return str.str();
            }

            id_map_t
            id_map_from_string(const std::string &const_str) {
                id_map_t result;

                string str = const_str;

                // Decode mode

                typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;

                typedef boost::tokenizer<boost::char_separator<char> > kv_tokenizer_t;

                boost::char_separator<char> sep(";");

                boost::char_separator<char> kv_sep(":");

                char chars[] = "[]";

                for (unsigned int i = 0; i < strlen(chars); ++i) {
                    str.erase(std::remove(str.begin(), str.end(), chars[i]), str.end());
                }

                tokenizer_t tokens(str, sep);

                for (tokenizer_t::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
                    string token = *tok_iter;

                    kv_tokenizer_t kvt(token, kv_sep);

                    for (kv_tokenizer_t::iterator kv_tok_iter = kvt.begin(); kv_tok_iter != kvt.end(); ++kv_tok_iter) {
                        std::string key_string = *kv_tok_iter;

                        id_t key = boost::lexical_cast<id_t>(key_string);

                        kv_tok_iter++;

                        std::string value_string = *kv_tok_iter;

                        id_set_t value = utils::sets::from_string<m3D::id_t>(value_string);

                        result.insert(std::pair<id_t, id_set_t>(key, value));
                    }
                }

                return result;
            }
        }
    }
}

#endif
