#!/usr/bin/python
# Filename: visitUtils.py

version = 'v1.3'

from visit import *
import glob
import os
import os.path
import string

# This module bundles python routines for handling Visit3D
# plotting more comfortably from the various visualization
# routines

# module variable holding a reference
datetime_annotation="not_set"

# Saves PNG file
# @param basename for the file
# @param is numbering appended?
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


# Closes the databases associated with the cluster files
# parameters follow the same rules as add_clusters_2D
def close_pattern(cluster_pattern):
    cluster_list=glob.glob(cluster_pattern)
    for cluster_file in cluster_list:
        CloseDatabase(cluster_file)
    return

# Adds a label file
# @param file
# @param label variable name
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
    a.textHeight1 = 0.035
    #a.specifyTextColor2 = 0
    #a.textColor2 = (0, 0, 255, 0)
    #a.textHeight2 = 0.03
    a.formatTemplate = "%g"
    SetPlotOptions(a)
    return

# Extracts the path part from the given filename
# @param filename
# @return path
def path(filename):
    path = dirname(filename)
    return path

# extracts the filename WITHOUT extension from the
# the given filename
# @param filename
# @return stripped filename
def naked_name(filename):
    base = s.basename(filename)
    stripped = os.path.splitext(filename)[0]
    return stripped

# Create a color table
# @param name
# @param control points
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

## Add a 2D text annotation of height 2%
# @param x position
# @param y position
# @param text message
def add_text_annotation(x,y,message):

    text="";
    try:
        text = GetAnnotationObject("Text2D1")

    except visit.VisItException:
        text = CreateAnnotationObject("Text2D")

    text.text = message;
    text.position = (x,y)
    text.height = 0.02
    return

## Add background gradient from dark (middle)
# to light gray (outside)
def add_background_gradient():
    a = AnnotationAttributes()
    a.backgroundMode = a.Gradient
    a.gradientBackgroundStyle = a.Radial
    a.gradientColor1 = (0,0,0,255)
    a.gradientColor2 = (192,192,192,255)
    SetAnnotationAttributes(a)
    return

# Extracts the date/time from a filename according
# to the oase 3D composite format and adds it to
# the currenty image
def add_datetime(filename):
    
    # OASE 3D
    # herz-oase-20110605t2355utc-0500m-bonnjue-3d-v01a.nc
    baseIndex=string.find(filename,"herz-oase")
    if baseIndex >= 0:
        year = filename[baseIndex+10:baseIndex+14]
        month = filename[baseIndex+14:baseIndex+16]
        day = filename[baseIndex+16:baseIndex+18]
        hour = filename[baseIndex+19:baseIndex+21]
        minute = filename[baseIndex+21:baseIndex+23]
        text = day+"."+month+"."+year+" "+hour+":"+minute+" UTC"
        add_text_annotation(0.725,0.95,text);
        return
    
    # OASE 2D
    #oase-20110622t2055z-1km-germany-2d-v01a.nc
    baseIndex=string.find(filename,"oase-")
    if baseIndex >= 0:
        year = filename[baseIndex+5:baseIndex+9]
        month = filename[baseIndex+9:baseIndex+11]
        day = filename[baseIndex+11:baseIndex+13]
        hour = filename[baseIndex+14:baseIndex+16]
        minute = filename[baseIndex+16:baseIndex+18]
        text = day+"."+month+"."+year+" "+hour+":"+minute+" UTC"
        add_text_annotation(0.725,0.95,text);
        return

    # RADOLAN
    #raa01-rx_10000-YYMMDDhhmm-dwd---bin.nc
    
    baseIndex=string.find(filename,"raa01-rx_10000")
    if baseIndex >= 0:
        year = filename[baseIndex+15:baseIndex+17]
        month = filename[baseIndex+17:baseIndex+19]
        day = filename[baseIndex+19:baseIndex+21]
        hour = filename[baseIndex+21:baseIndex+23]
        minute = filename[baseIndex+23:baseIndex+25]
        text = day+"."+month+".'"+year+" "+hour+":"+minute+" UTC"
        add_text_annotation(0.725,0.95,text);
        return

# End of visitUtils.py