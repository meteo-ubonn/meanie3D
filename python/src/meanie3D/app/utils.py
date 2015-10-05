__author__ = "juergen.simon@uni-bonn.de"

# -------------------------------------------------------------------
# Define some executables to be called from python
# -------------------------------------------------------------------

import glob
import json
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

##
# Tests if the dictionary has the given key. If so, it retrieves the
# value and returns it. If not, returns None.
# \param:dictionary
# \param:key
# \param value or None
#
def safeGet(dict,key):
    if key in dict and dict[key]:
        return dict[key]
    else:
        return None

##
# Sets the given object values along a keypath separated
# by periods. Example: axes2D.yAxis.title.units.
# \param:object
# \param:keypath
# \param:value
#
def setValueForKeyPath(object,keypath,value):
    keys = keypath.split(".")
    if len(keys) > 1:
        if type(object) is dict:
            nextObject = object[keys[0]]
        else:
            nextObject = getattr(object,keys[0])
        keys.pop(0)
        remainingPath = ".".join(keys)
        setValueForKeyPath(nextObject,remainingPath,value)
    else:
        # The specific enumerated values in Visit are configured
        # as strings in configuration dictionary, but are int values
        # in visit objects. Try to set an enumerated value when hitting
        # this combination
        value = None
        if type(value) is str and type(getattr(object,keypath)) is int:
            value = getattr(object,value)

        if type(object) is dict:
            object[keypath] = value
        else:
            setattr(object,keypath,value)

    return

##
# Gets the given object values along a keypath separated
# by periods. Example: axes2D.yAxis.title.units
# \param:object
# \param:keypath
# \return:value
#
def getValueForKeyPath(object,keypath):
    keys = keypath.split(".")
    if len(keys) > 1:
        if type(object) is dict:
            nextObject = safeGet(object,keys[0])
        else:
            nextObject = getattr(object,keys[0])
        keys.pop(0)
        remainingPath = ".".join(keys)
        return getValueForKeyPath(nextObject,remainingPath)
    else:
        # The specific enumerated values in Visit are configured
        # as strings in configuration dictionary, but are int values
        # in visit objects. Try to set an enumerated value when hitting
        # this combination
        value = None
        if type(object) is dict:
            value = safeGet(object,keypath)
        else:
            value = getattr(object,keypath)
        return value

##
# Sets the values from the dictionary on the given object
# \param:object
# \param:dictionary
# TODO: resolve key paths
#
def setValuesFromDictionary(object,dictionary):
    for key,value in dictionary.items():
        if type(value) is list:
            setValueForKeyPath(object,key,tuple(value))
        else:
            setValueForKeyPath(object,key,value)
    return
