__all__ = ["clusters", "tracks","utils"]
__have_visit = None
__visitPath__ = None
__visitImportPath__ = None

import sys
import meanie3D.app.utils

# Import visit package
__visitPath__, __visitImportPath__ = meanie3D.app.utils.findVisitPaths()
if __visitImportPath__ and __visitPath__:
    __have_visit = True
    sys.path.append(__visitImportPath__)
else:
    __have_visit = False
    sys.stderr.write("ERROR:could not find visit import path. Visualisation is switched off.")
