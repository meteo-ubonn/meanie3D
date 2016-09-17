__author__ = 'Juergen Simon'
__email__ = 'juergen.simon@uni-bonn.de'
__version__ = '1.6.0'
__url__ = 'http://git.meteo.uni-bonn.de/projects/meanie3d'
__all__ = ['app', 'visualisation', 'resources']

import os.path
import sys

def getVersion():
    '''
    :return:meanie3D package version
    '''
    from . import __version__
    return __version__

def getHome():
    '''
    meanie3D package location
    :return:
    '''
    return os.path.abspath(os.path.dirname(__file__))

def appendSystemPythonPath():
    '''
    Returns the system's python path to import outside modules.
    '''
    system_python_path = ":/Library/Frameworks/SQLite3.framework/Versions/C/Python/2.7:/Library/Frameworks/GEOS.framework/Versions/3/Python/2.7:/Library/Frameworks/GDAL.framework/Versions/1.11/Python/2.7/site-packages:/usr/local/lib/python2.7/site-packages/scipy-0.14.0-py2.7-macosx-10.9-x86_64.egg:/usr/local/lib/python2.7/site-packages/meanie3D-1.6.0-py2.7.egg:/usr/local/Cellar/python/2.7.11/Frameworks/Python.framework/Versions/2.7/lib/python27.zip:/usr/local/Cellar/python/2.7.11/Frameworks/Python.framework/Versions/2.7/lib/python2.7:/usr/local/Cellar/python/2.7.11/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin:/usr/local/Cellar/python/2.7.11/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac:/usr/local/Cellar/python/2.7.11/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages:/usr/local/Cellar/python/2.7.11/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk:/usr/local/Cellar/python/2.7.11/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old:/usr/local/Cellar/python/2.7.11/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload:/usr/local/lib/python2.7/site-packages:/usr/local/Cellar/matplotlib/1.5.1/libexec/lib/python2.7/site-packages:/usr/local/Cellar/numpy/1.9.2/libexec/nose/lib/python2.7/site-packages:/usr/local/lib/python2.7/site-packages/wx-3.0-osx_cocoa:/Library/Python/2.7/site-packages/numpy-override:/Library/Python/2.7/site-packages/matplotlib-override:/Library/Python/2.7/site-packages:/Library/Python/2.7/site-packages/pip-1.5.2-py2.7.egg".strip()
    paths = system_python_path.split(':')
    for path in paths:
        if path:
            sys.path.append(path)
