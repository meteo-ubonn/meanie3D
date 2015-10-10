##
# This modules contains a number of utility methods for visualisation
# purposes, such as commonly used plots, setup of color tables, adding
# annotations to plots etc.
#

import glob
import os
import os.path
import string
from os.path import basename
from os.path import dirname
from subprocess import call
import sys
import visit

# Own modules
import meanie3D.app.utils
import meanie3D.app.external

# make sure external commands are available
meanie3D.app.external.find_ext_cmds(['convert','composite'])

# ---------------------------------------------------------
#
# Global configuration handling
#
# ---------------------------------------------------------

# Global variable to help making sure the global visit configuration
# is only ran twice to avoid adding multiple color tables etc.
did_global_config_execute=False

##
# Creates any configuration items regarded as global, such
# as creating color tables and named views.
# \param:configuration
#
def runGlobalVisitConf(configuration):
    # Only run this once
    global did_global_config_execute
    if did_global_config_execute:
        return

    # Color tables
    if getValueForKeyPath(configuration,'postprocessing.visit.colorTables'):
        createColorTables(configuration,'postprocessing.visit.colorTables')

    saveWindowAttributes = getValueForKeyPath(configuration,'postprocessing.visit.saveWindowAttributes');
    if saveWindowAttributes:
        s = visit.GetSaveWindowAttributes()
        updateVisitObjectFromDictionary(s,saveWindowAttributes)
        visit.SetSaveWindowAttributes(s)

    # Level or verbosity
    visit.SuppressMessages(2)
    visit.SuppressQueryOutputOn()


    did_global_config_execute = True
    return

# ---------------------------------------------------------
#
# Plotting functions, helpers
#
# ---------------------------------------------------------

##
# Add the data of one variable from the given file to the current window.
# \param:databaseFile
# \param:configuration
#   {
#       "variable":"X",
#       "PseudocolorAttributes": {
#               see visit.PseudocolorAttributes()
#       }
#   }
#
def addPseudolorPlot(databaseFile,configuration):
    variable = getValueForKeyPath(configuration,'variable')
    attributes = getValueForKeyPath(configuration,'PseudocolorAttributes')
    if variable and attributes:
        # Open data file
        visit.OpenDatabase(databaseFile)
        # Add the plot
        visit.AddPlot("Pseudocolor", variable)
        p = visit.PseudocolorAttributes()
        updateVisitObjectFromDictionary(p,attributes)
        visit.SetPlotOptions(p)
        # Threshold?
        threshold = getValueForKeyPath(configuration,"ThresholdAttributes")
        if threshold:
            visit.AddOperator("Threshold")
            t = visit.ThresholdAttributes();
            updateVisitObjectFromDictionary(t,threshold)
            visit.SetOperatorOptions(t)
    return

##
# Adds a vector plot from the given file.
# \param:databaseFile
# \param:configuration
#   {
#       "variable":"X",
#       "VectorAttributes": {
#               see visit.VectorAttributes()
#       }
#   }
#
#
def addVectorPlot(databaseFile,configuration):
    variable = getValueForKeyPath(configuration,'variable')
    attributes = getValueForKeyPath(configuration,'VectorAttributes')
    if variable and attributes:
        visit.OpenDatabase(databaseFile)
        visit.AddPlot("Vector",variable)
        p=visit.VectorAttributes();
        updateVisitObjectFromDictionary(p,attributes)
        visit.SetPlotOptions(p)
    return

# -------------------------------------------------------------------
# Adds a label plot
# \param:file
# \param:configuration
#   {
#       "variable":"X",
#       "LabelAttributes": {
#               see visit.LabelAttributes()
#       }
#   }
#
#
# -------------------------------------------------------------------
def addLabelPlot(file,configuration):
    variable = getValueForKeyPath(configuration,'variable')
    attributes = getValueForKeyPath(configuration,'VectorAttributes')
    if variable and attributes:
        visit.OpenDatabase(file)
        visit.AddPlot("Label",variable)
        a = visit.LabelAttributes()
        updateVisitObjectFromDictionary(a,configuration)
        visit.SetActivePlots(visit.GetNumPlots()-1)
        visit.SetPlotOptions(a)
    return

##
# Add a 2D text annotation of height 2%
# \param:x position
# \param:y position
# \param:text message
#
def addTextAnnotation(x,y,message):
    try:
        textAnnotation = visit.GetAnnotationObject("Text2D1")
    except visit.VisItException:
        textAnnotation = visit.CreateAnnotationObject("Text2D")

    if textAnnotation:
        textAnnotation.text = message;
        textAnnotation.position = (x,y)
        textAnnotation.height = 0.02
    return

##
# Sets the view up with the given perspective object from the configuration.
# \param configuration Depending on what you wish to use, you can use all keys
# found in GetView2D() or GetView3D() respectively.
#
def setView(configuration,path):
    viewConfig = meanie3D.app.utils.getValueForKeyPath(configuration,path)
    if viewConfig:
        if 'windowCoords' in viewConfig or 'viewportCoords' in viewConfig:
            view = visit.GetView2D()
            updateVisitObjectFromDictionary(view,viewConfig)
            visit.SetView2D(view)
        else:
            view = visit.GetView3D()
            updateVisitObjectFromDictionary(view,viewConfig)
            visit.SetView3D(view)
    return

##
# Plots mapdata according to the configuration given. Note that $ variables
# will be replaced with the value found in the environment.
# \param:configuration Configuration dictionary
# \param:path Path to the map configuration
#
def addPseudocolorPlots(databaseFile,configuration,path):
    plots = getValueForKeyPath(configuration,path())
    if plots:
        visit.OpenDatabase(databaseFile)
        for plot in plots:
            addPseudolorPlot(databaseFile,plot)
    return

##
# Set annotation attributes from defaults and configuration.
# \param:configuration
#
def setAnnotations(configuration,path):
    ca = meanie3D.app.utils.getValueForKeyPath(configuration,path)
    if ca:
        a = visit.GetAnnotationAttributes()
        updateVisitObjectFromDictionary(a,ca)
        visit.SetAnnotationAttributes(a)
    return

##
# Creates a list of colortables from the configuration.
# \param:configuration
# \param:paht to list of colortables
#
def createColorTables(configuration,path):
    colorTables = getValueForKeyPath(configuration,path())
    if colorTables:
        for colorTable in colorTables:
            colors = colorTable['colors']
            positions = colorTable['positions']
            ccpl = visit.ColorControlPointList()
            for i in range(0,len(positions)):
                controlPoint = visit.ColorControlPoint()
                controlPoint.colors = colors[i]
                controlPoint.position = positions[i]
                ccpl.AddControlPoints(controlPoint)
            name = colorTable['name']
            visit.AddColorTable(name, ccpl)


##
# Closes the databases of all files matching the given pattern.
# \param:pattern
#
def close_pattern(pattern):
    list = glob.glob(pattern)
    for file in list:
        visit.CloseDatabase(file)
    return

##
# Plots mapdata according to the configuration given. Note that $ variables
# will be replaced with the value found in the environment.
# \param:configuration Configuration dictionary
# \param:path Path to the map configuration
#
def plotMapdata(configuration,path):
    mapConfig = getValueForKeyPath(configuration,path)
    if mapConfig:
        # Find the map data file
        mapFile = os.path.expandvars(mapConfig['mapDataFile'])

        # Check if data file exists
        if os.path.exists(mapFile):
            addPseudocolorPlots(mapFile, configuration, path+".plots")
        else:
            print "ERROR:could not find map file at " + mapFile
    return

# ---------------------------------------------------------
#
# File utilities
#
# ---------------------------------------------------------

##
# Extracts the path part from the given filename
# \param:filename
# \returns path
#
def path(filename):
    path = dirname(filename)
    return path

##
# extracts the filename WITHOUT extension from the
# the given filename
# \param:filename
# \returns stripped filename
# -------------------------------------------------------------------
def naked_name(filename):
    stripped = os.path.splitext(filename)[0]
    return stripped

# ---------------------------------------------------------
#
# Handling .png exports / images
#
# ---------------------------------------------------------

##
# Iterates over the given view objects, sets each in turn and
# saves an image. If there is only one view, the resulting image
# filenames are 'basename_<number>.png'. If there are more than
# one view, the resulting image filenames look like
# 'p<viewNum>_basename_<number>.png'
# \param:views
# \param:basename
#
def saveImagesForViews(views,basename):
    if views:
        for i in range(0,len(views)):
            view = views[i]
            setView(view)
            if len(views) > 1:
                filename = "p%d_%s" % (i,basename)
            else:
                filename = basename
            saveImage(filename,1)
    else:
        saveImage(basename,1)
    return

##
# Saves the current window to an image file.
# \param:basename Filename of the image. An underscore is appended.
# \param:progressive Is number appended?
#
def saveImage(basename,progressive):
    s = visit.GetSaveWindowAttributes()
    s.progressive = progressive
    s.filename = basename + "_"
    visit.SetSaveWindowAttributes(s)
    visit.SaveWindow()
    return


##
# Creates movies for the formats given from the images
# based on the basename given. Depending on the views,
# the filenames used start with 'p<viewNum>_basename_' or
# 'basename_'. See also saveImagesForViews(..).
# \param: views
# \param: basename
# \param extension, determining the movie format. Must be understood
# by the 'convert' program. Examples: "gif","mpeg", ...
#
def createMovieForViews(views,basename,format):
    if views:
        for i in range(0,len(views)):
            if len(views) > 1:
                movie_fn = "p%d_%s" % (i,basename)
            else:
                movie_fn = basename
            image_fn = movie_fn + "_"
            create_movie(image_fn,movie_fn + os.path.extsep + format)
    else:
        image_fn = basename+ "_"
        create_movie(image_fn,basename + os.path.extsep + format)


##
# \param: views
# \param: basename
# \param: movie formats ("gif","m4v" etc.)
#
def createMoviesForViews(views, basename, formats):
    for k in range(0,len(formats)):
        createMovieForViews(views, basename, formats[k])


##
# Checks if the file with the given number exists. Takes
# perspectives into account. If any perspective exists,
# all images from this number are deleted
# \param:configuration
# \param:basename ('source_','tracking_' etc.)
# \param:number number of file
#
def delete_images(conf,basename,image_count):
    number_postfix = str(image_count).rjust(4,'0') + ".png";
    result = False
    if 'PERSPECTIVES' in conf.keys():
        perspective_nr = 1
        for perspective in conf['PERSPECTIVES']:
            fn = "p"+str(perspective_nr)+"_"+basename+"_"+number_postfix
            if (os.path.exists(fn)):
                os.remove(fn)
            perspective_nr = perspective_nr + 1
    else:
        fn = basename+"_"+number_postfix
        if (os.path.exists(fn)):
            os.remove(fn)

##
# Checks if the image with the given number exists. Takes
# perspectives into account
# \param:views - visit view configurations
# \param:basename ('source_','tracking_' etc.)
# \param:number number of file
# \returns "all","none" or "partial"
#
def images_exist(views, basename, image_count):
    number_postfix = str(image_count).rjust(4,'0') + ".png";
    num_found = 0
    if views:
        for i in range(0,len(views)):
            view = views[i]
            fn = None
            if len(view) > 1:
                fn = "p"+str(i)+"_"+basename+"_"+number_postfix
            else:
                fn = basename+"_"+number_postfix
            if (os.path.exists(fn)):
                num_found = num_found + 1

        if num_found == 0:
            return "none"
        elif num_found == image_count * len(views):
            return "all"
        else:
            return "partial"

##
# Creates a movie from a .png image series.
# \param:basename of the image series
# \param:name of the resulting movie (with extension)
#
def create_movie(basename,moviename):
    print "Creating movie '" +moviename+"' from files '" + basename+"*.png ..."
    args = "-limit memory 4GB -delay 50 -quality 100 -dispose Background %s *.png %s" % (basename,moviename)
    meanie3D.app.external.execute_command('convert',args,False)
    print "done."
    return

# ---------------------------------------------------------
#
# Key-value coding, dictionary utils etc.
#
# ---------------------------------------------------------

##
# Visit's objects are a bit of an arse. Standard getattr does
# not seem to work in a lot of times.
# \param:object
# \param:key
# \returns value retrieved from visit object
#
def getVisitValue(object,key):
    value = getattr(object,key)
    if not value:
        # This may require a Getter
        #print "Attempting to resolve value for %s via getter" % key
        getter = getattr(object,'Get'+meanie3D.app.utils.capitalize(key))
        if getter:
            value = getter()
            #print "Success: object value for %s is %s" % (key,str(value))
    return value


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
        objectValue = getVisitValue(object,keypath)

        # The specific enumerated values in Visit are configured
        # as strings in configuration dictionary, but are int values
        # in visualisation objects. Try to set an enumerated value when hitting
        # this combination
        if type(value) != type(objectValue):
            value = objectValue

        if type(object) is dict:
            object[keypath] = value
        else:
            try:
                setattr(object,keypath,value)
            except:
                try:
                    # try a setter
                    #print "Attempting to set value via setter"
                    setter = getattr(object,'Set'+meanie3D.app.utils.capitalize(keypath))
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
    return meanie3D.app.utils.getValueForKeyPath(object,keypath)


##
# Sets the values from the dictionary on the given object
# \param:object
# \param:dictionary
#
def updateVisitObjectFromDictionary(object,dictionary):
    for key,value in dictionary.items():
        if type(value) is list:
            setValueForKeyPath(object,key,tuple(value))
        else:
            setValueForKeyPath(object,key,value)
    return

# ---------------------------------------------------------
#
# Deprecated code to be refactored or phased out
#
# ---------------------------------------------------------

# -------------------------------------------------------------------
# Extracts the date/time from a filename according
# to the oase 3D composite format and adds it to
# the currenty image
# -------------------------------------------------------------------
@PendingDeprecationWarning
def add_datetime(conf,filename):
    # Only use the leaf node
    fn = basename(filename);

    # OASE 3D
    # herz-oase-20110605t2355utc-0500m-bonnjue-3d-v01a.nc
    baseIndex=string.find(fn,"herz-oase")
    if baseIndex >= 0:
        year = fn[baseIndex+10:baseIndex+14]
        month = fn[baseIndex+14:baseIndex+16]
        day = fn[baseIndex+16:baseIndex+18]
        hour = fn[baseIndex+19:baseIndex+21]
        minute = fn[baseIndex+21:baseIndex+23]
        text = day+"."+month+"."+year+" "+hour+":"+minute+" UTC"
        addTextAnnotation(0.725,0.95,text);
        return

    # OASE 2D/3D
    #oase-20110622t2055z-1km-germany-2d-v01a.nc
    #oase-20130620t2100z-1km-germany-3d-test2
    baseIndex=string.find(fn,"oase-")

    if baseIndex >= 0:
        year = fn[baseIndex+5:baseIndex+9]
        month = fn[baseIndex+9:baseIndex+11]
        day = fn[baseIndex+11:baseIndex+13]
        hour = fn[baseIndex+14:baseIndex+16]
        minute = fn[baseIndex+16:baseIndex+18]
        text = day+"."+month+"."+year+" "+hour+":"+minute+" UTC"
        addTextAnnotation(0.725,0.95,text);
        return

    # RADOLAN
    #raa01-rx_10000-YYMMDDhhmm-dwd---bin.nc

    baseIndex=string.find(fn,"raa01-rx_10000")
    if baseIndex >= 0:
        year = fn[baseIndex+15:baseIndex+17]
        month = fn[baseIndex+17:baseIndex+19]
        day = fn[baseIndex+19:baseIndex+21]
        hour = fn[baseIndex+21:baseIndex+23]
        minute = fn[baseIndex+23:baseIndex+25]
        text = day+"."+month+".'"+year+" "+hour+":"+minute+" UTC"
        addTextAnnotation(0.725,0.95,text);
        return

# -------------------------------------------------------------------
# Creates a series of images combined from two series of images
# as one. It removes the date stamp on the left side and cuts
# the right side such, that no legend is visible
#
# \param:basename of the left image series
# \param:basename of the right image series
# \param:basename of the combined image series
# -------------------------------------------------------------------
@PendingDeprecationWarning
def create_dual_panel(basename_left,basename_right,basename_combined):
    left_files=sorted(glob.glob(basename_left+"*.png"))
    right_files=sorted(glob.glob(basename_right+"*.png"))
    if len(left_files) != len(right_files):
        print "ERROR:the two image series "+basename_left+"*.png and "+basename_right+"*.png have different lengths!"
        return

    # create a small image to blank the datestamp with
    call("/usr/local/bin/convert -size 280x50 xc:white dateblind.png",shell=True)

    # Iterate over the series
    for i in range(len(left_files)):

        # create backdrop
        combined=basename_combined+str(i)+".png"
        meanie3D.app.external.execute_command("convert","-size 1852x1024 xc:white "+combined,shell=True)

        # copy images and blank date
        meanie3D.app.external.execute_command("composite","-geometry +826+0 "+right_files[i]+" "+combined+" "+combined)
        meanie3D.app.external.execute_command("composite","-geometry +0+0 "+left_files[i]+" "+combined+" "+combined)
        meanie3D.app.external.execute_command("composite","-geometry +735+20 dateblind.png "+combined+" "+combined)

    # delete the dateblind
    call("rm dateblind.png",shell=True)
    return

