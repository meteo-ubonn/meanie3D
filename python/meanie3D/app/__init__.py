__author__ = 'juergen.simon@uni-bonn.de'
__all__ = ["external", "postprocessing","utils","tracking"]
__have_visit__ = None
__visitPath__ = None
__visitImportPath__ = None

import sys
from meanie3D.app.utils import findVisitPaths

# Import visit package
__visitPath__, __visitImportPath__ = findVisitPaths()
if __visitImportPath__ and __visitPath__:
    __have_visit__ = True
    sys.path.append(__visitImportPath__)
else:
    __have_visit__ = False
    sys.stderr.write("ERROR:could not find visit import path. Visualisation is switched off.")
