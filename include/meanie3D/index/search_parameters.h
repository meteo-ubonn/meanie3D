#ifndef M3D_SEARCHPARAMS_H
#define M3D_SEARCHPARAMS_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

namespace m3D {

    typedef enum {
        SearchTypeKNN,
        SearchTypeRange
    } SearchType;

    class SearchParameters
    {
    private:

        SearchType  m_searchType;

        /** Private default constructor 
         */
        SearchParameters() {};

    protected:

        /** Constructor is protected to prevent direct construction 
         */
        SearchParameters( SearchType type ) : m_searchType(type) {};

    public:

        /** Destructor */
        virtual ~SearchParameters() {};

        /** @returns type
         */
        const SearchType search_type() const { return m_searchType; };

    };

    /** K nearest neighbours search.
     */
    template <typename T>
    class KNNSearchParams : public SearchParameters
    {
    public:
        size_t      k;
        vector<T>   resolution;

        KNNSearchParams(const size_t _k)
        : SearchParameters(SearchTypeKNN), k(_k)
        {
        };

        KNNSearchParams(const size_t _k, const vector<T> &_resolution)
        : SearchParameters(SearchTypeKNN), k(_k), resolution(_resolution)
        {
        };

        KNNSearchParams( const KNNSearchParams& other )
        : SearchParameters( other.search_type() )
        {
            k = other.k;

            resolution = other.resolution;
        }

        KNNSearchParams operator=(const KNNSearchParams& other) 
        {
            KNNSearchParams copy(other);

            return copy;
        }

        ~KNNSearchParams() {};
    };

    /** Range search
     */
    template <typename T>
    class RangeSearchParams : public SearchParameters
    {
    public:
        vector<T>   bandwidth;

        RangeSearchParams(const vector<T>& bandwidth) : SearchParameters(SearchTypeRange)
        {
            this->bandwidth = bandwidth;
        };

        RangeSearchParams( const RangeSearchParams& other ) : SearchParameters( other.search_type() )
        {
            bandwidth = other.bandwidth;
        }

        RangeSearchParams operator=(const RangeSearchParams& other) 
        {
            RangeSearchParams copy(other);

            return copy;
        }

        ~RangeSearchParams() {};
    };
}
    
#endif
