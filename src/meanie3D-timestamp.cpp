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

#pragma mark -
#pragma mark Definitions

/** Feature-space data type 
 */
typedef double T;

/** Filename formats 
 */
typedef enum {
    TimestampFormatRadolan,
    TimestampFormatOASE2D,
    TimestampFormatOASE3D,
    TimestampFormatHErZTB4,
    TimestampFormatAutomatic,
    
} TimestampFormat;


#pragma mark -
#pragma mark Command line parsing

void parse_commmandline(program_options::variables_map vm,
                        string &source,
                        TimestampFormat &ts_format)
{
    if ( vm.count("source") == 0 )
    {
        cerr << "Missing 'source' argument" << endl;
        
        exit( 1 );
    }
    
    source = vm["source"].as<string>();
    
    std::string fmt = vm["format"].as<string>();
    
    if (fmt == "oase-2d")
    {
        ts_format = TimestampFormatOASE2D;
    }
    else if (fmt == "oase-3d")
    {
        ts_format = TimestampFormatOASE3D;
    }
    else if (fmt == "radolan")
    {
        ts_format = TimestampFormatRadolan;
    }
    else if (fmt == "tb4")
    {
        ts_format = TimestampFormatHErZTB4;
    }
    
    else if (fmt == "auto")
    {
        ts_format = TimestampFormatAutomatic;
    }
    else
    {
        cerr << "Unknown value '" << fmt << "' for --format. Allowed values are oase-2d oase-3d radolan" << endl;
        exit(-1);
    }
}

#pragma mark -
#pragma mark Helper Methods

typedef struct tm tm_t;

tm_t * initialize_tm_struct()
{
    tm_t *ts = (tm_t *) malloc(sizeof(tm_t));
    
    ts->tm_sec = 0;
    ts->tm_min = 0;		/* minutes after the hour [0-59] */
    ts->tm_hour = 0;     /* hours since midnight [0-23] */
    ts->tm_mday = 0;     /* day of the month [1-31] */
    ts->tm_mon = 0;      /* months since January [0-11] */
    ts->tm_year = 0;     /* years since 1900 */
    ts->tm_wday = 0;     /* days since Sunday [0-6] */
    ts->tm_yday = 0;     /* days since January 1 [0-365] */
    ts->tm_isdst = 0;    /* Daylight Savings Time flag */
    ts->tm_gmtoff = 0;   /* offset from CUT in seconds */
    ts->tm_zone = NULL;  /* timezone abbreviation */
    
    return ts;
}

// some constants

const std::string RADOLAN_PREFIX = "raa01-rx_10000-";
const std::string RADOLAN_FORMAT = "%y%m%d%H%M";
const std::string RADOLAN_EXAMPLE = "1307010740";

const std::string OASE_2D_PREFIX = "oase-";
const std::string OASE_2D_FORMAT = "%Y%m%d%tH%Mz";
const std::string OASE_2D_EXAMPLE = "20110622t1500z";

const std::string OASE_3D_PREFIX = "herz-oase-";
const std::string OASE_3D_FORMAT = "%Y%m%dt%H%M";
const std::string OASE_3D_EXAMPLE = "20110605t1555";

const std::string TB4_PREFIX = "tb_";
const std::string TB4_FORMAT = "%Y%m%d%H";
const std::string TB4_EXAMPLE = "2011060618";

/** 
 */
timestamp_t parse_timestamp(std::string filename, TimestampFormat format)
{
    timestamp_t result = 0;
    
    std::string prefix;
    std::string dateformat;
    std::string example;
    
    switch (format)
    {
        case TimestampFormatAutomatic:
        {
            if (boost::starts_with(filename,RADOLAN_PREFIX))
            {
                // raa01-rx_10000-1307010740-dwd---bin.nc
                
                prefix = RADOLAN_PREFIX;
                dateformat = RADOLAN_FORMAT;
                example = RADOLAN_EXAMPLE;
            }
            else if (boost::starts_with(filename,OASE_2D_PREFIX))
            {
                prefix = OASE_2D_PREFIX;
                dateformat = OASE_2D_FORMAT;
                example = OASE_2D_EXAMPLE;
            }
            else if (boost::starts_with(filename,OASE_3D_PREFIX))
            {
                prefix = OASE_3D_PREFIX;
                dateformat = OASE_3D_FORMAT;
                example = OASE_3D_EXAMPLE;
            }
            else if (boost::starts_with(filename,TB4_PREFIX))
            {
                prefix = TB4_PREFIX;
                dateformat = TB4_FORMAT;
                example = TB4_EXAMPLE;
            }
            else
            {
                cerr << "ERROR:could not detect format for filename " << filename << ". Please advise format with --format switch." << endl;
                exit(-1);
            }
        } break;
            
        case TimestampFormatRadolan:
        {
            // raa01-rx_10000-1307010740-dwd---bin.nc
            
            prefix = RADOLAN_PREFIX;
            dateformat = RADOLAN_FORMAT;
            example = RADOLAN_EXAMPLE;
            
        } break;

        case TimestampFormatOASE2D:
        {
            // oase-20110622t1500z-1km-germany-2d-v01a.nc
            
            prefix = OASE_2D_PREFIX;
            dateformat = OASE_2D_FORMAT;
            example = OASE_2D_EXAMPLE;
            
        } break;

        case TimestampFormatOASE3D:
        {
            // herz-oase-20110605t1555utc-0500m-bonnjue-3d-v01a.nc
            
            prefix = OASE_3D_PREFIX;
            dateformat = OASE_3D_FORMAT;
            example = OASE_3D_EXAMPLE;
            
        } break;

        case TimestampFormatHErZTB4:
        {
            // tb_2011060618_150ghz.nc
            
            prefix = TB4_PREFIX;
            dateformat = TB4_FORMAT;
            example = TB4_EXAMPLE;
            
        } break;
    }
    
    std::string str = filename.substr( prefix.size(),example.size());
    
    tm_t *ts = initialize_tm_struct();
    
    if (strptime( str.c_str(), dateformat.c_str(), ts) == NULL)
    {
        cerr << "Error parsing datetime string " << str << endl;
        exit(-1);
    }
    
    // HErZ-TB4 has no minute/second
    
    if (format == TimestampFormatHErZTB4)
    {
        ts->tm_min = 0;
        ts->tm_sec = 0;
    }

    result = timegm(ts);
    
    delete ts;

    
    return result;
}

#pragma mark -
#pragma mark MAIN

/* MAIN
 */
int main(int argc, char** argv)
{
    using namespace m3D;
    
    // Declare the supported options.
    
    program_options::options_description desc("Extracts a timestamp from filename(s) and adds it as a NetCDF-variable 'time' to the file(s).");
    desc.add_options()
    ("help", "Produces this help.")
    ("source", program_options::value<string>(), "A single file or a directory to be processed. Only files ending in .nc will be processed.")
    ("format", program_options::value<string>()->default_value("auto"), "Timestamp format specifier <oase-2d,oase-3d,radolan,auto(default)>");
    
    program_options::variables_map vm;
    
    try
    {
        program_options::store( program_options::parse_command_line(argc, argv, desc), vm);
        program_options::notify(vm);
    }
    catch (std::exception &e)
    {
        cerr << "ERROR:parsing command line caused exception: " << e.what() << endl;
        cerr << "Check meanie3D-trackplot --help for command line options" << endl;
        exit(-1);
    }
    
    if ( vm.count("help")==1 || argc < 2 )
    {
        cout << desc << "\n";
        return 1;
    }
    
    // Evaluate user input
    
    string source_path;
    TimestampFormat format;
    
    namespace fs = boost::filesystem;
    
    try
    {
        parse_commmandline(vm, source_path, format);
    }
    catch (const std::exception &e)
    {
        cerr << e.what() << endl;
        exit(-1);
    }
    
    typedef set<fs::path> fset_t;
    
    fset_t files;
    
    if (fs::is_directory(source_path))
    {
        fs::directory_iterator dir_iter(source_path);
        fs::directory_iterator end;
        
        while (dir_iter != end)
        {
            fs::path f = dir_iter->path();
            
            if (fs::is_regular_file(f) && fs::extension(f) == ".nc")
            {
                //cout << "Adding " << f.generic_string() << endl;
                files.insert(f);
            }
            else
            {
                cout << "Skipping " << f.generic_string() << endl;
            }
            
            dir_iter++;
        }
    }
    else
    {
        fs::path f = fs::path(source_path);
        
        std::string extension = fs::extension(f);
        
        if (fs::is_regular_file(f) && extension == ".nc")
        {
            files.insert(f);
        }
    }
    
    fset_t::iterator it;
    
    boost::progress_display *progress = NULL;
    
    if ( files.size() > 1 )
    {
        progress = new progress_display(files.size());
    }
    
    for (it = files.begin(); it != files.end(); ++it)
    {
        if (progress != NULL)
        {
            progress->operator++();
        }
        
        fs::path path = *it;
        
        std::string fn = path.generic_string();
        
        // Add a dimension 'time'
        
        timestamp_t ts = parse_timestamp(boost::filesystem::basename(fn), format);
        ::cfa::utils::netcdf::add_time(fn, ts);
        
    }
    
    if (progress != NULL)
    {
        delete progress;
    }

    return 0;
};
