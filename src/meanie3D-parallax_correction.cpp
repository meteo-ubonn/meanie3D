//
//  meanie3D-detect
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/date_time.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <sstream>

#include <meanie3D/meanie3D.h>
#include <radolan/radolan.h>

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <exception>
#include <locale>
#include <limits>
#include <stdlib.h>
#include <netcdf>
#include <time.h>

using namespace std;
using namespace boost;
using namespace netCDF;
using namespace m3D;

namespace fs = boost::filesystem;

#pragma mark -
#pragma mark Definitions

/** Feature-space data type
 */
typedef double T;

#pragma mark -
#pragma mark Command line parsing

void parse_commmandline ( program_options::variables_map vm,
                          string &filename )
{
	if ( vm.count ( "file" ) == 0 ) {
		cerr << "Missing 'file' argument" << endl;

		exit ( 1 );
	}

	filename = vm["file"].as<string>();
}

#pragma mark -
#pragma mark Worker Methods

/** This is a c++ translation of the method sent to me by Marianne Koenig (in F90):
 *
 * <cite>
 * Subroutine does a parallax correction for something seen at some
 * height in a position lat/lon by the satellite given by satheight,
 * satlat, satlon. The new coordinates are returned in latcorr and loncorr
 * </cite>
 *
 * @param satheight (REAL): height of the satellite in km
 * @param satlat (REAL): subsatellite latitude (deg, N is positive)
 * @param satlon (REAL): subsatellite longitude (deg, E is positive)
 * @param height (REAL): height of the cloud (km)
 * @param lat (REAL): latitude of the satellite pixel (N is positive)
 * @param lon (REAL): longitude of the satellite pixel (E is positive)
 */
template <typename T>
void parallax ( T satheight, T satlat, T satlon, T height, T lat, T lon, T& latcorr, T& loncorr )
{
	T dpi;
	T radius_eq;
	T radius_pole;
	T radius_ratio;
	T mean_radius;
	T dheight;
	T alat,alon;
	T asatlat,asatlon;
	T satlat_geod,satlon_geod;
	T xsat,ysat,zsat;
	T xsurf,ysurf,zsurf;
	T alat_geod;
	T radius_surf;
	T radius_ratio_local;
	T xdiff,ydiff,zdiff;
	T xfact,zen;
	T e1,e2,e3;
	T corr;
	T xcorr,ycorr,zcorr;

	dpi = 3.14159265;

	// varius earth radius information

	radius_eq = 6378.077;
	dheight = satheight;
	radius_pole = 6356.577;
	radius_ratio = radius_eq/radius_pole;
	mean_radius = 0.5* ( radius_eq+radius_pole );
	zdiff = 0.0;

//     angle conversion to radians

	asatlat = satlat * dpi/180.0;
	asatlon = satlon * dpi/180.0;
	alat = lat*dpi/180.0;
	alon = lon*dpi/180.0;

//     cartesian coordinates for the satellite
//     satlat_geod is the geodetic satellite latitude

	satlat_geod = atan ( tan ( asatlat ) * ( radius_ratio*radius_ratio ) );
	xsat = dheight * cos ( satlat_geod ) * sin ( asatlon );
	ysat = dheight * cos ( satlat_geod );
	zsat = dheight * cos ( satlat_geod ) * cos ( asatlon );

//     cartesian coordinates of the surface point

	alat_geod = atan ( tan ( alat ) *radius_ratio*radius_ratio );
	radius_surf = radius_eq/sqrt ( cos ( alat_geod ) *cos ( alat_geod ) + radius_ratio*radius_ratio * sin ( alat_geod ) *sin ( alat_geod ) );
	xsurf = radius_surf * cos ( alat_geod ) * sin ( alon );
	ysurf = radius_surf * sin ( alat_geod );
	zsurf = radius_surf * cos ( alat_geod ) * cos ( alon );

//     compute new radius ratio depending on height

	radius_ratio_local = ( ( radius_eq+height ) / ( radius_pole+height ) );
	radius_ratio_local *= radius_ratio_local;

//     Satellite minus surface location

	xdiff = xsat - xsurf;
	ydiff = ysat - ysurf;
	zdiff = zsat - zdiff;

//     compute local zenith angle

	xfact = sqrt ( xdiff*xdiff + ydiff*ydiff + zdiff*zdiff );
	zen = ( xdiff*xsurf+ydiff*ysurf+zdiff*zsurf ) / ( mean_radius*xfact );
	zen = acos ( zen );
	zen = zen*180.0/dpi;

//     equation to solve for the line of sight at height Z

	e1 = xdiff*xdiff + radius_ratio_local*ydiff*ydiff + zdiff*zdiff;
	e2 = 2.0 * ( xsurf*xdiff + radius_ratio_local*ysurf*ydiff + zsurf*zdiff );
	e3 = xsurf*xsurf + zsurf*zsurf + radius_ratio_local*ysurf*ysurf - ( radius_eq+height ) * ( radius_eq+height );

	corr = ( sqrt ( e2*e2 - 4.0*e1*e3 ) - e2 ) /2.0/e1;

//     corrected surface coordinates

	xcorr = xsurf + corr*xdiff;
	ycorr = ysurf + corr*ydiff;
	zcorr = zsurf + corr*zdiff;

//     convert back to latitude and longitude

	latcorr = atan ( ycorr/sqrt ( xcorr*xcorr + zcorr*zcorr ) );
	latcorr = atan ( tan ( latcorr ) /radius_ratio*radius_ratio ) * 180.0/dpi;

	loncorr = atan2 ( xcorr,zcorr ) * 180.0/dpi;
}

template <typename T>
T** allocate_array ( size_t dim_y, size_t dim_x )
{
	T **array = new T*[dim_y];
	for ( int i = 0; i < dim_y; ++i )
		array[i] = new T[dim_x];
	return array;
}

#define deallocate_array(array,dim) for (int i=0; i<dim; i++) delete[] array[i]; delete[] array;

/** Corrects the parallax on all seviri satellite variables
 * in the national 2D OASE composite
 * @param in_path path to the netcdf file to be corrected. The data
 * in the file is overwritten.
 */
void correct_parallax ( boost::filesystem::path in_path )
{
	// Some constants
	const double SAT_LON = 9.5;			// longitude of METEOSAT-9
	const double SAT_LAT = 0.0;			// latitute of METEOSAT-9
	const double SAT_HEIGHT = 35785.83;	// height of METEOSAT-9 [km]
	
	#define dim_x 900
	#define dim_y 900

	try 
	{
		NcFile file ( in_path.generic_string(), NcFile::write );

		typedef std::multimap<std::string,NcVar> vmap_t;

		vmap_t variables = file.getVars();

		// Cloud-Top-Height is needed as input

		static int cloud_top_height[dim_y][dim_x];

		vmap_t::iterator fi = variables.find ( "msevi_l2_nwcsaf_cth" );
		if ( fi == variables.end() ) {
			cerr << "ERROR: could not find cloud top height (msevi_l2_nwcsaf_cth) variable" << endl;
			return;
		}

		fi->second.getVar ( &cloud_top_height[0][0] );

		/*
		cout << endl;
		for ( size_t iy = 0; iy < dim_y; iy++ )
		{
			if (iy % 10 == 0)
			{
				for ( size_t ix = 0; ix < dim_x; ix++ )
				{
					if (ix % 10 == 0)
					{
						char c = (cloud_top_height[iy][ix] > 0) ? '*' : ' ';
						cout << c;
					}
				}

				cout << endl;
			}
		}
		cout << endl;
		*/

		// create a variable for input and one for output

		static int corrected_iy[dim_y][dim_x];
		static int corrected_ix[dim_y][dim_x];

		// initialize output data with flag to find pixels
		// later that have not been set

		for ( size_t iy = 0; iy < dim_y; iy++ ) {
			for ( size_t ix = 0; ix < dim_x; ix++ ) {
				corrected_ix[iy][ix] = 0;
				corrected_iy[iy][ix] = 0;
			}
		}

		// Coordinate system for lat/lon transformation
		RDCoordinateSystem rcs ( RD_RX );

		// correct the parallax

		for ( size_t iy = 0; iy < dim_y; iy++ ) {
			for ( size_t ix = 0; ix < dim_x; ix++ ) {
				RDGridPoint gp = rdGridPoint ( ix,iy );

				// get lat/lon for this pixel
				RDGeographicalPoint coord = rcs.geographicalCoordinate ( gp );

				// Get Marianne Koenig's correction values

				T cth = boost::numeric_cast<float> ( cloud_top_height[iy][ix] ) / 1000.0f;

				T lat_corrected = 0;
				T lon_corrected = 0;

				parallax<double> ( SAT_HEIGHT, SAT_LAT, SAT_LON, cth, coord.latitude, coord.longitude, lat_corrected, lon_corrected );

				// Figure out the grid point again and set
				// data at corrected position

				RDGeographicalPoint coord_coorected;
				coord_coorected.latitude = -lat_corrected;
				coord_coorected.longitude = -lon_corrected;

				bool is_inside = false;
				RDGridPoint gp_corrected = rcs.gridPoint ( coord_coorected,is_inside );

				// TODO: what if two pixels are moved to the same place?
				// The way things are now is 'last write wins'

				if ( is_inside ) {
					corrected_ix[iy][ix] = gp_corrected.ix;
					corrected_iy[iy][ix] = gp_corrected.iy;
				}
			}
		}

		static int input_data[dim_y][dim_x];
		static int output_data[dim_y][dim_x];

		for ( vmap_t::iterator vi = variables.begin(); vi != variables.end(); vi++ ) 
		{
			NcVar variable = vi->second;
			
			int fill_value = 0;
			try
			{
				// Get the official _FillValue value if the
				// variable has one
				NcVarAtt fillValue = variable.getAtt("_FillValue");
				fillValue.getValues(&fill_value);
			}
			catch (netCDF::exceptions::NcException e)
			{
				// if not, put the value just outside the valid range
				int valid_min = std::numeric_limits<int>::min();
				fill_value = valid_min - 1;
			}

			if ( boost::starts_with ( variable.getName(), "msevi_" ) ) {
				// Initialize arrays

				for ( size_t iy = 0; iy < dim_y; iy++ ) {
					for ( size_t ix = 0; ix < dim_x; ix++ ) {
						input_data[iy][ix] = 0;
						output_data[iy][ix] = fill_value;
					}
				}

				// Read the satellite variable

				variable.getVar ( &input_data[0][0] );

				// apply the correction derived from cloud top height

				for ( size_t iy = 0; iy < dim_y; iy++ ) {
					for ( size_t ix = 0; ix < dim_x; ix++ ) {
						// get corrected indicees

						size_t iy_corr = corrected_iy[iy][ix];
						size_t ix_corr = corrected_ix[iy][ix];

						// copy data over

						output_data[iy_corr][ix_corr] = input_data[iy][ix];
					}
				}

				// TODO: post-processing of points that got no values

				// Write data back

				variable.putVar ( &output_data[0][0] );
			}
		}

	} catch ( netCDF::exceptions::NcException &e ) {
		cerr << e.what() << endl;
		return;
	}

}

#pragma mark -
#pragma mark MAIN

/* MAIN
 */
int main ( int argc, char** argv )
{
	using namespace m3D;

	// Declare the supported options.

	program_options::options_description desc ( "Applies parallax correction to mseviri satellite data in OASE composite files." );
	desc.add_options()
	( "help", "Produces this help." )
	( "file,f", program_options::value<string>(), "A single file or a directory to be processed. Only files ending in .nc will be processed." );

	program_options::variables_map vm;

	try {
		program_options::store ( program_options::parse_command_line ( argc, argv, desc ), vm );
		program_options::notify ( vm );
	} catch ( std::exception &e ) {
		cerr << "ERROR:parsing command line caused exception: " << e.what() << endl;
		cerr << "Check meanie3D-trackplot --help for command line options" << endl;
		exit ( -1 );
	}

	if ( vm.count ( "help" ) ==1 || argc < 2 ) {
		cout << desc << "\n";
		return 1;
	}

	// Evaluate user input

	string source_path;

	try {
		parse_commmandline ( vm, source_path );
	} catch ( const std::exception &e ) {
		cerr << e.what() << endl;
		exit ( -1 );
	}

	typedef set<fs::path> fset_t;

	fset_t files;

	if ( fs::is_directory ( source_path ) ) {
		fs::directory_iterator dir_iter ( source_path );
		fs::directory_iterator end;

		while ( dir_iter != end ) {
			fs::path f = dir_iter->path();

			if ( fs::is_regular_file ( f ) && fs::extension ( f ) == ".nc" ) {
				//cout << "Adding " << f.generic_string() << endl;
				files.insert ( f );
			} else {
				cout << "Skipping " << f.generic_string() << endl;
			}

			dir_iter++;
		}
	} else {
		fs::path f = fs::path ( source_path );

		std::string extension = fs::extension ( f );

		if ( fs::is_regular_file ( f ) && extension == ".nc" ) {
			files.insert ( f );
		}
	}

	fset_t::iterator it;

	boost::progress_display *progress = NULL;

	if ( files.size() > 1 ) {
		progress = new progress_display ( files.size() );
	}

	for ( it = files.begin(); it != files.end(); ++it ) {
		if ( progress != NULL ) {
			progress->operator++();
		}

		boost::filesystem::path path = *it;

		// Correct

		try {
			correct_parallax ( path );
		} catch ( std::exception &e ) {
			cerr << "Exception processing " << path.filename().generic_string() << endl;
			cerr << "Cause: " << e.what() << endl;
		}
	}

	if ( progress != NULL ) {
		delete progress;
	}

	return 0;
};
