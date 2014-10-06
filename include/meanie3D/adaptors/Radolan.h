#ifndef M3D_ADAPTORS_RADOLAN_H
#define M3D_ADAPTORS_RADOLAN_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/exceptions.h>

#include <string>
#include <netcdf>
#include <radolan/radolan.h>

using namespace Radolan;

namespace m3D { 

    /** Converts the radolan file at path into a CF-Metadata compliant NetCDF-File.
     * @param radolanPath full path to the radolan file
     * @param netcdfPath full path to the netcdf file to be created
     * @param threshold minimum treshold for values to make it in the NetCDF file
     * @param mode NcFile::Mode for opening the netcdf file with
     * @param omitOutside @see RDReadScan
     * @return NCFile* NetCDF-Filehandler
     * @throw CFFileConversionException
     */
    netCDF::NcFile * CFConvertRadolanFile( const char* radolanPath,
                                  const char* netcdfPath,
                                  const RDDataType *threshold = NULL,
                                  netCDF::NcFile::FileMode mode = netCDF::NcFile::replace, 
                                  bool omitOutside = true)
        throw (CFFileConversionException);

    /** Converts the radolan file at path into a CF-Metadata compliant NetCDF-File. 
     * @param radolanPath full path to the radolan file
     * @param netcdfPath full path to the netcdf file to be created
     * @param threshold minimum treshold for values to make it in the NetCDF file
     * @param mode NcFile::Mode for opening the netcdf file with
     * @return NCFile* NetCDF-Filehandler
     * @throw CFFileConversionException
     */
    netCDF::NcFile * CFConvertRadolanScan( RDScan *scan,
                                  const char* netcdfPath, 
                                  const RDDataType *threshold = NULL,
                                  netCDF::NcFile::FileMode mode = netCDF::NcFile::write)
        throw (CFFileConversionException);

    /** Simple function to get a visual rep of the file with ascii characters 
     * on terminal.
     * @param NcFile* netcdf Radolan file in CF-Metadata format 
     * @param print values every latVertices points in y
     * @param print values every lonVertices points in x
     */
    void CFPrintConvertedRadolanScan( netCDF::NcFile *file, int latVertices=20, int lonVertices=20 );

    /** CF-Metadata 'standard_name' for the given scan type 
     * @param scanType
     * @return standard_name
     */
    const char * CFRadolanDataStandardName( RDScanType scanType );

}

#endif
