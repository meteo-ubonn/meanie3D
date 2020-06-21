FROM debian:stable

ENV PACKAGES=
RUN apt-get -y update --fix-missing
RUN apt-get -y upgrade 
RUN apt-get -y dist-upgrade
RUN apt-get -y install software-properties-common
RUN apt-get -y update 
RUN apt-get -y install \
wget git cmake \
gcc g++ python python-pip \
libboost-all-dev libflann1.9 libflann-dev blitz++ \
shapelib libhdf5-dev netcdf-bin libnetcdf-dev libnetcdf-c++4 libnetcdf-c++4-dev zlib1g zlib1g-dev 
RUN pip install setuptools netcdf4 external utils

# Visualisation
# RUN pip install Cython h5py netcdf4
# RUN apt-get -y --fix-missing  install gnuplot vtk6 libvtk6-dev
# RUN wget --quiet http://portal.nersc.gov/project/visit/releases/2.10.0/visit2_10_0.linux-x86_64-rhel6-wmesa.tar.gz
# RUN wget --quiet http://portal.nersc.gov/project/visit/releases/2.10.0/visit-install2_10_0
# RUN chmod a+x visit-install2_10_0
# RUN echo "1" | ./visit-install2_10_0 2.10.0 linux-x86_64-rhel6-wmesa /usr/local/visit
# ENV VISIT_EXECUTABLE=/usr/local/visit/bin/visit
# RUN rm -rf visit*

# libradolan
RUN git clone https://github.com/JuergenSimon/radolan.git
RUN cd radolan && cmake . && make install && cd .. && rm -rf radolan

# Meanie3D
RUN git clone --recurse-submodules --depth=1 https://github.com/JuergenSimon/meanie3D
WORKDIR /meanie3D
RUN git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*" && git fetch --all
RUN git checkout --track remotes/origin/dockerize && git pull
RUN git pull && cmake -DFOR_DOCKER=YES . && make install 

# Cleanup
WORKDIR /
RUN rm -rf meanie3D
RUN apt-get remove -y wget git cmake
RUN apt autoremove -y

# Prepare for runtime
RUN mkdir /data
ENV LD_LIBRARY_PATH=/usr/local/lib
ENTRYPOINT ["/usr/local/bin/meanie3D"]