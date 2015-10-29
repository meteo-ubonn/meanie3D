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

#ifndef M3D_SEARCHPARAMS_H
#define M3D_SEARCHPARAMS_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

namespace m3D {

    typedef enum
    {
        SearchTypeKNN,
        SearchTypeRange
    } SearchType;

    class SearchParameters
    {
    private:

        SearchType m_searchType;

        /** Private default constructor 
         */
        SearchParameters()
        {
        };

    protected:

        /** Constructor is protected to prevent direct construction 
         */
        SearchParameters(SearchType type) : m_searchType(type)
        {
        };

    public:

        /** Destructor */
        virtual ~SearchParameters()
        {
        };

        /** @returns type
         */
        const SearchType search_type() const
        {
            return m_searchType;
        };

    };

    /** K nearest neighbours search.
     */
    template <typename T>
    class KNNSearchParams : public SearchParameters
    {
    public:
        size_t k;
        vector<T> resolution;

        KNNSearchParams(const size_t _k)
        : SearchParameters(SearchTypeKNN), k(_k)
        {
        };

        KNNSearchParams(const size_t _k, const vector<T> &_resolution)
        : SearchParameters(SearchTypeKNN), k(_k), resolution(_resolution)
        {
        };

        KNNSearchParams(const KNNSearchParams& other)
        : SearchParameters(other.search_type())
        {
            k = other.k;

            resolution = other.resolution;
        }

        KNNSearchParams operator=(const KNNSearchParams& other)
        {
            KNNSearchParams copy(other);

            return copy;
        }

        ~KNNSearchParams()
        {
        };
    };

    /** Range search
     */
    template <typename T>
    class RangeSearchParams : public SearchParameters
    {
    public:
        vector<T> bandwidth;

        RangeSearchParams(const vector<T>& bandwidth) : SearchParameters(SearchTypeRange)
        {
            this->bandwidth = bandwidth;
        };

        RangeSearchParams(const RangeSearchParams& other) : SearchParameters(other.search_type())
        {
            bandwidth = other.bandwidth;
        }

        RangeSearchParams operator=(const RangeSearchParams& other)
        {
            RangeSearchParams copy(other);

            return copy;
        }

        ~RangeSearchParams()
        {
        };
    };
}

#endif
