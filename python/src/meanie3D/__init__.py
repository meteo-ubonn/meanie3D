__author__ = 'juergen.simon@uni-bonn.de'
__version__ = '1.5.4'

# Check if all required external commands are available
from meanie3D.app import external

try:
    external.find_ext_cmds(["meanie3D-detect","meanie3D-track","meanie3D-cfm2vtk","meanie3D-trackstats"])
except IOError:
    print "ERROR:could not locate one or more of meanie3D binaries. Are the C++ binaries installed?"
    exit(-1)

##
# \return meanie3D package version
#
def getVersion():
    return __version__

##
# \return meanie3D package location
#
def getHome():
    import meanie3D
    return meanie3D.__file__
