#!/usr/bin/python

# Parameters
SOURCE_DIR = "SOURCE_DIR_P"
M3D_HOME = "MEANIE3D_HOME_P"

# Import modules
import sys
sys.path.append(M3D_HOME+"/scripts/python-modules")
import glob
import os
import time
import visit2D
import meanie3D_visit_utils
from subprocess import call

# Silent
SuppressMessages(True)
SuppressQueryOutputOn()

# Set view and annotation attributes

WITH_BACKGROUND_GRADIENT=False
WITH_TOPOGRAPHY=False
WITH_RIVERS_AND_BOUNDARIES=True
WITH_DATETIME=True
OPACITY=1.0

print "Setting annotation attributes:"
visit2D.set_annotations()

#print "Cleaning up *.vtk *.vtr"
return_code=call("rm -f *.vtk *.vtr", shell=True)

#print "Cleaning up "+VAR_NAME+"*.png"
return_code=call("rm -f *.png", shell=True)

# Setting 3D view parameters
print "Setting 2D view parameters"
visit2D.set_view_to_radolan();

# Add gray/black background gradient
if WITH_BACKGROUND_GRADIENT:
    print "Setting background gradient"
    meanie3D_visit_utils.add_background_gradient();

if WITH_TOPOGRAPHY:
    print "Creating colortables"
    meanie3D_visit_utils.create_topography_colortable()

# General parameters
a = GetAnnotationAttributes()
a.databaseInfoFlag = 0
a.timeInfoFlag = 0
a.userInfoFlag = 0
a.legendInfoFlag = 1
a.axes2D.visible=1
a.axes2D.xAxis.title.visible=0
a.axes2D.yAxis.title.visible=0
a.axes2D.autoSetScaling=0
SetAnnotationAttributes(a)

# Glob the netcdf directory
print "Processing files in directory " + SOURCE_DIR
netcdf_files = glob.glob(SOURCE_DIR+"/*.nc");

# Counter for re-setting the engine periodically
image_count = 0

VARIABLES=("cband_radolan_rx",
           "linet_oase_tl",
           "msevi_l2_cmsaf_cph",
           "msevi_l2_cmsaf_cwp",
           "msevi_l2_cmsaf_cot",
           "msevi_l2_cmsaf_reff",
           "msevi_l2_nwcsaf_cth",
           "msevi_l2_nwcsaf_ct",
           "msevi_l2_nwcsaf_cma",
           "msevi_l15_hrv",
           "msevi_l15_vis006",
           "msevi_l15_vis008",
           "msevi_l15_wv_062",
           "msevi_l15_wv_073",
           "msevi_l15_ir_016",
           "msevi_l15_ir_039",
           "msevi_l15_ir_087",
           "msevi_l15_ir_097",
           "msevi_l15_ir_108",
           "msevi_l15_ir_120",
           "msevi_l15_ir_134")

#VARIABLES=("msevi_l15_wv_062",
#           "msevi_l15_wv_073")

for VAR_NAME in VARIABLES:
    
    print "------------------------------------------------------------------------"
    print VAR_NAME
    print "------------------------------------------------------------------------"
    
    debug_count = 0;

    for netcdf_file in netcdf_files:
        
        print "Proessing " + netcdf_file
        
        if WITH_TOPOGRAPHY:
            print "-- Adding topography data --"
            visit2D.add_topography("national_topo_2D")
        
        if WITH_RIVERS_AND_BOUNDARIES:
            print "-- Adding map data --"
            visit2D.add_map_rivers("national")
            visit2D.add_map_borders("national")
        
        if WITH_DATETIME:
            print "-- Adding timestamp --"
            meanie3D_visit_utils.add_datetime(netcdf_file)

        # now plot the data
        OpenDatabase(netcdf_file)
        
        if VAR_NAME == "cband_radolan_rx":
            
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"hot_desaturated",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="hot_desaturated"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,65
            SetPlotOptions(p)
            AddOperator("Threshold")
            t = ThresholdAttributes();
            t.lowerBounds=(15)
            SetOperatorOptions(t)

        elif VAR_NAME == "linet_oase_tl":
        
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"YlOrRd",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="YlOrRd"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=1,20
            SetPlotOptions(p)
            AddOperator("Threshold")
            t = ThresholdAttributes();
            t.lowerBounds=(1)
            SetOperatorOptions(t)
        
        elif VAR_NAME == "msevi_l2_cmsaf_cph":

            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"Purples",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="Purples"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,2
            SetPlotOptions(p)
            AddOperator("Threshold")
            t = ThresholdAttributes();
            t.lowerBounds=(0)
            SetOperatorOptions(t)

        elif VAR_NAME == "msevi_l2_cmsaf_cwp":
            
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"hot_desaturated",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="hot_desaturated"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,5000
            SetPlotOptions(p)
            AddOperator("Threshold")
            t = ThresholdAttributes();
            t.lowerBounds=(250)
            SetOperatorOptions(t)

        elif VAR_NAME == "msevi_l2_cmsaf_cot":
    
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"hot_desaturated",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="hot_desaturated"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,250
            SetPlotOptions(p)
            AddOperator("Threshold")
            t = ThresholdAttributes();
            t.lowerBounds=(100)
            SetOperatorOptions(t)

        elif VAR_NAME == "msevi_l2_cmsaf_reff":
            
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"hot_desaturated",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="hot_desaturated"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,50
            SetPlotOptions(p)
            AddOperator("Threshold")
            t = ThresholdAttributes();
            t.lowerBounds=(5)
            SetOperatorOptions(t)

        elif VAR_NAME == "msevi_l2_nwcsaf_cth":
            
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"bluehot",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="bluehot"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=2500,15000
            SetPlotOptions(p)
            AddOperator("Threshold")
            t = ThresholdAttributes();
            t.lowerBounds=(2500)
            SetOperatorOptions(t)

        elif VAR_NAME == "msevi_l2_nwcsaf_ct":
            
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"levels",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="levels"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,20
            SetPlotOptions(p)
            AddOperator("Threshold")
            t = ThresholdAttributes();
            t.lowerBounds=(7)
            t.upperBounds=(14)
            SetOperatorOptions(t)

        elif VAR_NAME == "msevi_l2_nwcsaf_cma":
            
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"levels",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="levels"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,5
            SetPlotOptions(p)
            AddOperator("Threshold")
            t = ThresholdAttributes();
            t.lowerBounds=(3)
            t.upperBounds=(4)
            SetOperatorOptions(t)

        elif VAR_NAME == "msevi_l15_hrv":
            
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"gray",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="gray"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,13
            SetPlotOptions(p)

        elif VAR_NAME == "msevi_l15_vis006" or VAR_NAME == "msevi_l15_vis008":
            
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"gray",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="gray"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,12
            SetPlotOptions(p)

        elif VAR_NAME == "msevi_l15_wv_062":
            
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"xray",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="xray"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,3.5
            SetPlotOptions(p)

        elif VAR_NAME == "msevi_l15_wv_073":
    
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"xray",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="xray"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,12
            SetPlotOptions(p)

        elif VAR_NAME == "msevi_l15_ir_016":
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"gray",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="gray"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,7
            SetPlotOptions(p)

        elif VAR_NAME == "msevi_l15_ir_039":
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"xray",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="xray"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,1
            SetPlotOptions(p)

        elif VAR_NAME == "msevi_l15_ir_087":
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"xray",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="xray"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=0,70
            SetPlotOptions(p)

        elif VAR_NAME == "msevi_l15_ir_097":
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"xray",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="xray"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=10,50
            SetPlotOptions(p)

        elif VAR_NAME == "msevi_l15_ir_108":
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"xray",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="xray"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=10,120
            SetPlotOptions(p)

        elif VAR_NAME == "msevi_l15_ir_120":
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"xray",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="xray"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=20,125
            SetPlotOptions(p)

        elif VAR_NAME == "msevi_l15_ir_134":
            visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"xray",OPACITY,1)
            p = PseudocolorAttributes()
            p.colorTableName="xray"
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=25,100
            SetPlotOptions(p)

        DrawPlots()
        CloseDatabase(netcdf_file)
        visit2D.close_topography();
        meanie3D_visit_utils.save_window(VAR_NAME+"_",1)
        DeleteAllPlots()
        ClearWindow()
        
        # Visit has a bug that causes it to stop
        # producing images after a certain amount
        # of runs. This will periodically restart
        # the rendering engine, thereby preventing
        # that from happening
        
        image_count = image_count + 1
    
        if image_count > 100:
            CloseComputeEngine()
            image_count = 0
    
    #debug_count = debug_count + 1;
    #    if debug_count == 10:
    #        break

    print "Creating movie ..."
    return_code=call("mkdir "+VAR_NAME, shell=True)
    return_code=call("mkdir "+VAR_NAME+"/images", shell=True)

    convert_cmd="/usr/local/bin/convert -limit memory 4GB -delay 50 -quality 100 "+VAR_NAME+"_*.png "+VAR_NAME+".gif"
    return_code=call(convert_cmd, shell=True)

    convert_cmd="/usr/local/bin/convert -limit memory 4GB -delay 50 -quality 100 "+VAR_NAME+"_*.png "+VAR_NAME+".m4v"
    return_code=call(convert_cmd, shell=True)

    return_code=call("mv "+VAR_NAME+"*.png "+VAR_NAME+"/images", shell=True)
    return_code=call("mv "+VAR_NAME+"*.gif "+VAR_NAME, shell=True)
    return_code=call("mv "+VAR_NAME+"*.m4v "+VAR_NAME, shell=True)

return_code=call("rm visitlog.py", shell=True)

print "Done."
quit()