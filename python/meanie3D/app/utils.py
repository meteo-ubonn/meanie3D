__author__ = "juergen.simon@uni-bonn.de"

# -------------------------------------------------------------------
# Define some executables to be called from python
# -------------------------------------------------------------------

import glob
import json
import os
import os.path
import shutil
import sys

import external

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
def getSafe(dict,key):
    if key in dict and dict[key]:
        return dict[key]
    else:
        return None

##
# Courtesy of http://stackoverflow.com/questions/1549641/how-to-capitalize-the-first-letter-of-each-word-in-a-string-python
#
def capitalize(line):
    return ' '.join(s[0].upper() + s[1:] for s in line.split(' '))

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
        # in visualisation objects. Try to set an enumerated value when hitting
        # this combination
        if type(value) is str and type(getattr(object,keypath)) is int:
            print "Attempting to resolve constant"
            value = getattr(object,value)

        if type(object) is dict:
            object[keypath] = value
        else:
            try:
                setattr(object,keypath,value)
            except:
                try:
                    # try a setter
                    setter = getattr(object,'Set'+capitalize(keypath))
                    if setter:
                        setter(value)
                except:
                    sys.stderr.write("Can't set value %s for key '%s' on object" % (str(value),keypath))
                    raise

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
            nextObject = getSafe(object,keys[0])
        else:
            nextObject = getattr(object,keys[0])
        keys.pop(0)
        remainingPath = ".".join(keys)
        return getValueForKeyPath(nextObject,remainingPath)
    else:
        # The specific enumerated values in Visit are configured
        # as strings in configuration dictionary, but are int values
        # in visualisation objects. Try to set an enumerated value when hitting
        # this combination
        value = None
        if type(object) is dict:
            value = getSafe(object,keypath)
        else:
            value = getattr(object,keypath)
        return value

##
# Sets the values from the dictionary on the given object
# \param:object
# \param:dictionary
#
def setValuesFromDictionary(object,dictionary):
    for key,value in dictionary.items():
        setValueForKeyPath(object,key,value)
    return

##
# Recursive directory search
# \param:path
# \param:filename
# \param:must have one component in the path to the result that matches this
# \return: fully qualified path to result or None
#
def find(path,filename,requiredComponent=None):
    if os.path.exists(path) and os.path.isdir(path):
        files = os.listdir(path)
        for file_ in files:
            if file_ == filename:
                if requiredComponent:
                    components = path.split(os.path.sep)
                    if not requiredComponent in components:
                        continue
                    return os.path.abspath(path + os.sep + filename)
                else:
                    return os.path.abspath(path + os.sep + filename)

    for f in files:
        full_path = os.path.abspath(path + os.sep + f)
        if os.path.isdir(full_path):
            result = find(full_path, filename, requiredComponent)
            if result:
                return result
    return None


##
# Finds the paths to the visit executable and python package
# \return (visit python package path or None, visit binary path or None)
#
def findVisitPaths():
    visitPath = None
    visitImportPath = None
    cmdMap = external.locateCommands(["visit"],True)
    if cmdMap['visit']:
        visitPath = os.path.abspath(os.path.join(cmdMap['visit'], os.pardir + os.sep + os.pardir))
        if visitPath:
            visitImportPath = find(visitPath,"site-packages")
    return visitPath,visitImportPath


##
# Converts a string into a bool using the following sets:
# True: ['true', '1', 't', 'y', 'yes']:
# False: ['false', '0', 'f', 'n', 'no']:
# \param:s
# \returns:True if s is in true value set, False if s is None or in false value set.
# \throws:ValueError s is neither None nor in either true or false value set.
def strToBool(s):
    if s:
        test = s.lower()
        if test in ['true', '1', 't', 'y', 'yes']:
            return True
        elif test in ['false', '0', 'f', 'n', 'no']:
            return False
        else:
            raise ValueError
    else:
        return False

