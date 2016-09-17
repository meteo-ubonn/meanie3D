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
    system_python_path = ":/Users/simon/anaconda/lib/python35.zip:/Users/simon/anaconda/lib/python3.5:/Users/simon/anaconda/lib/python3.5/plat-darwin:/Users/simon/anaconda/lib/python3.5/lib-dynload:/Users/simon/anaconda/lib/python3.5/site-packages:/Users/simon/anaconda/lib/python3.5/site-packages/Sphinx-1.4.1-py3.5.egg:/Users/simon/anaconda/lib/python3.5/site-packages/aeosa:/Users/simon/anaconda/lib/python3.5/site-packages/setuptools-23.0.0-py3.5.egg".strip()
    paths = system_python_path.split(':')
    for path in paths:
        if path:
            sys.path.append(path)
