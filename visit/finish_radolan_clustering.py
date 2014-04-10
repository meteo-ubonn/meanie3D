#!/usr/bin/python

# Import modules
import sys
import os
from subprocess import call

# create loops
print "Creating animated gifs ..."
convert_cmd="/usr/local/bin/convert -limit memory 4GB -delay 50 -quality 100 source_*.png source.gif"
return_code=call(convert_cmd, shell=True)

convert_cmd="/usr/local/bin/convert -limit memory 4GB -delay 50 -quality 100 tracking_*.png tracking.gif"
return_code=call(convert_cmd, shell=True)

print "Creating mpegs ..."
convert_cmd="/usr/local/bin/convert -limit memory 4GB -delay 50 -quality 100 source_*.png source.mpeg"
return_code=call(convert_cmd, shell=True)

convert_cmd="/usr/local/bin/convert -limit memory 4GB -delay 50 -quality 100 tracking_*.png tracking.mpeg"
return_code=call(convert_cmd, shell=True)

# clean up
print "Cleaning up ..."
return_code=call("mkdir images", shell=True)
return_code=call("mv *.png images", shell=True)
