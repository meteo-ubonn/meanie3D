#!/usr/bin/python

# -------------------------------------------------------------------
# Filename: utils.py
# @author Juergen Simon (juergen_simon@mac.com)
#
# This modules contains a number of utility methods used throughout
# the project.
# -------------------------------------------------------------------

import glob
import os
import os.path
import string
from os.path import basename
from subprocess import call
# Own modules
from meanie3D import external
# Visit modules

# make sure external commands are available
external.find_ext_cmds(['convert','composite'])

# -------------------------------------------------------------------
# This module bundles python routines for handling Visit3D
# plotting more comfortably from the various visualization
# routines
# module variable holding a reference
# -------------------------------------------------------------------
datetime_annotation="not_set"

# ------------------------------------------------------------------------------
# Add the data of one variable from the given file to the current window.
# \param filename
# \param variable name
# \param configuration A configuration dictionary containing a complete
# configuration for a Pseudocolor plot (see visit python API)
#
# ------------------------------------------------------------------------------
def add_pseudocolor(vtk_file,configuration):
    OpenDatabase(vtk_file)
    AddPlot("Pseudocolor", configuration['variable'])
    p = PseudocolorAttributes()
    p.__dict__.update(configuration.__dict__)
    SetPlotOptions(p)
    return

# -------------------------------------------------------------------
# Saves PNG file
# @param basename for the file
# @param is numbering appended?
# TODO: make image dimensions configurable
# -------------------------------------------------------------------
def save_window(basename,progressive):
    s = GetSaveWindowAttributes()
    s.outputToCurrentDirectory = 1
    s.outputDirectory = "."
    s.fileName=basename
    s.width = 1024
    s.height = 1024
    s.quality = 100
    s.progressive=progressive
    SetSaveWindowAttributes(s)
    SaveWindow()
    return

# -------------------------------------------------------------------
# Closes the databases associated with the cluster files
# parameters follow the same rules as add_clusters_2D
# -------------------------------------------------------------------
def close_pattern(cluster_pattern):
    cluster_list=glob.glob(cluster_pattern)
    for cluster_file in cluster_list:
        CloseDatabase(cluster_file)
    return

# -------------------------------------------------------------------
# Extracts the path part from the given filename
# @param filename
# @return path
# -------------------------------------------------------------------
def path(filename):
    path = dirname(filename)
    return path

# -------------------------------------------------------------------
# extracts the filename WITHOUT extension from the
# the given filename
# @param filename
# @return stripped filename
# -------------------------------------------------------------------
def naked_name(filename):
    base = s.basename(filename)
    stripped = os.path.splitext(filename)[0]
    return stripped

# -------------------------------------------------------------------
# Create a color table
# @param name
# @param control points
# -------------------------------------------------------------------
def make_rgb_colortable(name, ct):
    ccpl = ColorControlPointList()
    index=0;
    for pt in ct:
        p = ColorControlPoint()
        p.colors = (pt[0], pt[1], pt[2], pt[3])
        p.position = 0.2 + (0.2 * index)
        ccpl.AddControlPoints(p)
        index=index+1
    AddColorTable(name, ccpl)
    return

# -------------------------------------------------------------------
# Add a 2D text annotation of height 2%
# @param x position
# @param y position
# @param text message
# -------------------------------------------------------------------
def add_text_annotation(x,y,message):
    text = "";
    try:
        text = GetAnnotationObject("Text2D1")

    except visit.VisItException:
        text = CreateAnnotationObject("Text2D")

    text.text = message;
    text.position = (x,y)
    text.height = 0.02
    return


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
        nextObject = getattr(object,keys[0])
        keys.pop(0)
        remainingPath = ".".join(keys)
        setValueForKeyPath(nextObject,remainingPath,value)
    else:
        # The specific enumerated values in Visit are configured
        # as strings in configuration dictionary, but are int values
        # in visit objects. Try to set an enumerated value when hitting
        # this combination
        if type(value) is str and type(getattr(object,keypath)) is int:
            setattr(object,keypath,getattr(object,value))
        else:
            setattr(object,keypath,value)
    return

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

# -------------------------------------------------------------------
# Add background gradient from dark (middle)
# to light gray (outside)
# -------------------------------------------------------------------
def add_background_gradient(configuration):
    a = AnnotationAttributes()
    if configuration['backgroundMode']:
        backgroundMode = a.__dict__.get(configuration['backgroundMode'])
        a.backgroundMode = backgroundMode

    if a.backgroundMode == a.Gradient:
        if configuration['gradientBackgroundStyle']:
            gradientBackgroundStyle = configuration['gradientBackgroundStyle']
            a.gradientBackgroundStyle = gradientBackgroundStyle

    a.__dict__.update(configuration.__dict__)

    SetAnnotationAttributes(a)
    return

# -------------------------------------------------------------------
# Creates a movie from a .png image series.
# @param basename of the image series
# -------------------------------------------------------------------
def create_movie(basename,moviename):
    print "Creating movie '" +moviename+"' from files '" + basename+"*.png ..."
    args = "-limit memory 4GB -delay 50 -quality 100 -dispose Background %s *.png %s" % (basename,moviename)
    external.execute_command('convert',args,False)
    print "done."
    return

# -------------------------------------------------------------------
# Creates the color table for clusters
# -------------------------------------------------------------------
@PendingDeprecationWarning("This code is specific to project OASE and will be deprecated")
def create_cluster_colortable(name):
    opacity=1*255
    rgb_colors = ((255, 0,      0,      opacity),   #plain red
                  (255, 255,    0,      opacity),   #yellow
                  (255, 0,      255,    opacity),   #purple
                  (0,   255,    0,      opacity),   #green
                  (0,   153,    255,    opacity))   #blue

    # Each control point is: (r,g,b,t)
    make_rgb_colortable(name,rgb_colors)
    return len(rgb_colors)

# -------------------------------------------------------------------
# Creates a colortable for visualizing the topography data from
# the mapstuff file
# -------------------------------------------------------------------
@PendingDeprecationWarning("This code is specific to project OASE and will be deprecated")
def create_topography_colortable():
    opacity=255
    rgb = ((0,   0,    100,    opacity),   #sea
           (120, 120,   83,    opacity),   #lowlands
           (49,  22,     0,    opacity),   #midlands
           (100, 100,  100,    opacity),   #foothills
           (255, 255,  255,    opacity))   #mountains

    pos=(0.0, 0.0375, 0.30, 0.50, 0.85)

    ccpl = ColorControlPointList()
    index=0;
    for pt in rgb:
        p = ColorControlPoint()
        p.colors = (pt[0], pt[1], pt[2], pt[3])
        p.position = pos[index]
        ccpl.AddControlPoints(p)
        index=index+1
    AddColorTable("topography", ccpl)
    return

# -------------------------------------------------------------------
# Extracts the date/time from a filename according
# to the oase 3D composite format and adds it to
# the currenty image
# -------------------------------------------------------------------
@PendingDeprecationWarning("This code is specific to project OASE and will be deprecated")
def add_datetime(filename):
    
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
        add_text_annotation(0.725,0.95,text);
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
        add_text_annotation(0.725,0.95,text);
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
        add_text_annotation(0.725,0.95,text);
        return

# -------------------------------------------------------------------
# Creates a series of images combined from two series of images
# as one. It removes the date stamp on the left side and cuts
# the right side such, that no legend is visible
#
# @param basename of the left image series
# @param basename of the right image series
# @param basename of the combined image series
# -------------------------------------------------------------------
@PendingDeprecationWarning("This code is specific to project OASE and will be deprecated")
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
        call("/usr/local/bin/convert -size 1852x1024 xc:white "+combined,shell=True)

        # copy images and blank date
        call("/usr/local/bin/composite -geometry +826+0 "+right_files[i]+" "+combined+" "+combined,shell=True)
        call("/usr/local/bin/composite -geometry +0+0 "+left_files[i]+" "+combined+" "+combined,shell=True)
        call("/usr/local/bin/composite -geometry +735+20 dateblind.png "+combined+" "+combined,shell=True)

    # delete the dateblind
    call("rm dateblind.png",shell=True)
    return

# -------------------------------------------------------------------
# Adds a label file
# @param file
# @param label variable name
# -------------------------------------------------------------------
@PendingDeprecationWarning("This code is specific to project OASE and will be deprecated")
def add_labels(file,variable):
    OpenDatabase(file)
    AddPlot("Label",variable)
    a = LabelAttributes()
    a.legendFlag = 0
    #a.showNodes = 0
    #a.showCells = 1
    a.restrictNumberOfLabels = 0
    #a.numberOfLabels = 200
    a.specifyTextColor1 = 0
    a.textColor1 = (255, 0, 0, 0)
    #a.textHeight1 = 0.035
    a.textHeight1 = 0.02
    #a.specifyTextColor2 = 0
    #a.textColor2 = (0, 0, 255, 0)
    #a.textHeight2 = 0.03
    a.formatTemplate = "%g"
    SetActivePlots(GetNumPlots()-1)
    SetPlotOptions(a)

    return


















# End of utils.py
