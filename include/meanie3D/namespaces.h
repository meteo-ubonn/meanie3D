#ifndef _M3D_NAMESPACE_DEFINITIONS_H_
#define _M3D_NAMESPACE_DEFINITIONS_H_

#include <iostream> // namespace std
#include <netcdf>   // namespace netCDF
#include <boost/locale.hpp>  // boost::locale
#include <cf-algorithms/cf-algorithms.h>

// This file manages the imports of namespaces into the 
// various namespaces defined in this software.
// That way all dependencies can be changed in one place
// should the need arise. Also, this is a great way of 
// actually seeing all namespaces and dependencies.

// dependencies of cfa namespaces to other namespaces

namespace m3D {

	using namespace ::std;
	using ::cfa::meanshift::FeatureSpace;
	using ::cfa::meanshift::FeatureSpaceIndex;
	using ::cfa::meanshift::SearchParameters;
	using ::cfa::meanshift::Kernel;

	namespace detection {
	};

	namespace utils { namespace visit {
		using namespace ::std;
		using namespace ::boost::locale;
	}};
};

#endif
