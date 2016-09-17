
#include <meanie3D/exceptions.h>

namespace m3D {

    /** Construct a new exception of this type
     */
    CFFileConversionException::CFFileConversionException(const char *message) {
        m_message = std::string(message);
    }

    CFFileConversionException::~CFFileConversionException() throw() {}

    /** @override 
     */
    const char *CFFileConversionException::what() const throw() {
        return m_message.c_str();
    }
}
