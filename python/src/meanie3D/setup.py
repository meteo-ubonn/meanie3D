from distutils.core import setup

setup(
    name='meanie3D',
    version='1.5.4',
    packages=['app','visualisation'],
    include_package_data=True,
    url='http://git.meteo.uni-bonn.de/projects/meanie3d',
    license='MIT License',
    author='Juergen Simon',
    author_email='juergen.simon@uni-bonn.de',
    description='A python script for running the meanie3D clustering and tracking algorithms as well as processing and visualising the results.'
)
