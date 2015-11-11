__author__ = 'Juergen Simon'
__email__ = 'juergen.simon@uni-bonn.de'
__version__ = '1.5.4'
__url__ = 'http://git.meteo.uni-bonn.de/projects/meanie3d'
__all__ = ['app', 'visualisation']

import os.path

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
