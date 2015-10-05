__author__ = "juergen.simon@uni-bonn.de"

# -------------------------------------------------------------------
# Define some executables to be called from python
# -------------------------------------------------------------------

import glob
import os
import os.path
import shutil

## Parses a JSON configuration file
# \param filename
# \returns configuration dictionary
def load_configuration(filename):
    json_data=open(filename)
    data = json.load(json_data)
    json_data.close()
    return data;


## Counts the number of netcdf files in the given
# directory
# \param directory
# \returns number of netcdf files
def number_of_netcdf_files(source_dir):
    netcdf_pattern = source_dir + "/*.nc"
    netcdf_list=sorted(glob.glob(netcdf_pattern))
    return len(netcdf_list)


## Creates an output filename based on given filename by
# appending -<slicenum>.nc at the end.
# \param basic filename
# \param slice num
# \returns filename-1.nc
def numbered_filename(filename,index):
    basename = os.path.basename(filename)
    return os.path.splitext(basename)[0]+"-"+str(index)+".nc"


## Deletes the directories 'log' and 'netcdf' underneath
# base path. Removes previous ones if they do exist
#
# \param base path
# -------------------------------------------------------------------
def create_ouput_directories(base_path):
    # base path
    if not os.path.exists(base_path):
        os.makedirs(base_path)

    # logs
    log_dir = base_path+"/log"
    if os.path.exists(log_dir):
        shutil.rmtree(log_dir)
    os.makedirs(log_dir)

    # netcdf results
    netcdf_dir = base_path+"/netcdf"
    if os.path.exists(netcdf_dir):
        shutil.rmtree(netcdf_dir)
    os.makedirs(netcdf_dir)

    return
