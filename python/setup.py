from setuptools import find_packages, setup
setup(
    name='meanie3D',
    version='1.5.4',
    packages=["meanie3D","meanie3D.app", "meanie3D.visualisation"],
    url='http://git.meteo.uni-bonn.de/projects/meanie3d',
    license='MIT License',
    author='Juergen Simon',
    author_email='juergen.simon@uni-bonn.de',
    zip_safe = True,
    description='A python script for running the meanie3D clustering and tracking algorithms as well as processing and visualising the results.',
    entry_points = {
        'console_scripts' : [
            'meanie3D = meanie3D.meanie3D_run:main'
        ]
    }
)
