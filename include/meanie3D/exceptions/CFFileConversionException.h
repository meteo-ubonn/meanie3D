#ifndef M3D_FILECONVERSIONEXCEPTION_H
#define M3D_FILECONVERSIONEXCEPTION_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <exception>
#include <string>

namespace m3D { 

    /** This exception is to be used when converting files from
     * different formats into the CF Metadata format or vica versa.
     */
    class CFFileConversionException : public std::exception
    {
    private:

        std::string   m_message;

    public:

        CFFileConversionException( const char* message );

        ~CFFileConversionException() throw();

        virtual const char* what() const throw();
    };
}

#endif
