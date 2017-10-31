.. _install:

======================================
Installation
======================================

.. contents:: Table of Contents
   :depth: 3

FOREST is as easy to install as any of its dependencies.
The main difficulty here is to install all the dependencies,
in the right order and in the right version, before installing
FOREST.

We start with the installation of FOREST assuming that
the dependencies are already installed, but if they are not
we suggest to directly go to the installation of the external
dependencies, and only after they are installed, install FOREST.

The dependencies are PETSc, SLEPc and deal.II. They should be installed
in this specific order to satisfy their inter-dependencies. The way to
install each one of them is summarized later.

FOREST
=================

In order to install FOREST, some dependencies should be installed first.
After the dependencies have been installed we proceed to download
and install FOREST by typing

.. code-block:: console

    git clone -b maint https://bitbucket.org/dream/forest forest
    cd forest
    cmake .
    make

PETSc
================

In order to download PETSc we write

.. code-block:: console

    git clone -b maint https://bitbucket.org/petsc/petsc petsc

In order to install PETSc we write

.. code-block:: console

    ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-fblaslapack --download-mpich

which should take care of installing the mpi and configuring
the installation of petsc, or

.. code-block:: console

    ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-fblaslapack --with-mpi=0

for installing without mpi.

.. note:: For version 8.4 of dealii to link it properly, we need to add the installation of hypre as follows

    .. code-block:: console

        ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-fblaslapack --download-mpich --download-hypre=1

Then we install it by typing

.. code-block:: console

    make all test

The last thing is now to add the folders to the path. We do this by adding the following lines to the ~/.bashrc file

.. code-block:: console

    export PETSC_DIR=/folder/with/petsc
    export PETSC_ARCH=arch-linux2-c-debug

SLEPc
================

SLEPc works in a very similar way to that of PETSc.
The main point here is to tell SLEPc where is PETSc, i.e.,

* Set the environment variable SLEPC_DIR to the full path  of the SLEPc home directory, for example,
  `export SLEPC_DIR=/home/username/slepc`
* In addition to this variable, PETSC_DIR
  and PETSC_ARCH  must also be set appropriately.

We can download SLEPc in a very similar way to PETSc, by typing

.. code-block:: console

    git clone -b maint https://bitbucket.org/slepc/slepc slepc


Then we go inside the SLEPc folder and we type:

.. code-block:: console

    ./configure

The last step is to install by typing

.. code-block:: console

    make all test

deal.II
==================

For installing deal.II, we start downloading the library
to the folder that we choose, which in my case is

.. code-block:: console

    cd /home/segonpin/lib
    git clone -b dealii-8.4 https://github.com/dealii/dealii.git dealii-8.4
    cd dealii
    mkdir build
    cd build

Now we prepare the directory for the installation and
export it to the path,

.. code-block:: console

    mkdir ../../dealii-8.4.bin
    export DEAL_II_DIR=../../dealii-8.4.bin

then I run the bash script `buildingdealii.sh`
by typing

.. code-block:: console

    sh buildingdealii.sh

where the bash script contains the following lines

.. code-block:: bash

    #/bin/bash

    # here I am using my environmental variable DEAL_II_DIR
    # as the installing directory. It can be any other path.
    echo "$DEAL_II_DIR"

    # prepare the cache for the cmake
    cmake -DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR \
        -DDEAL_II_WITH_FUNCTIONPARSER=ON \
        -DDEAL_II_COMPONENT_DOCUMENTATION=ON \
        -DMPI_CXX_COMPILER=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx \
        -DDEAL_II_WITH_PETSC=ON -DPETSC_DIR=$PETSC_DIR -DPETSC_ARCH=$PETSC_ARCH \
        -DDEAL_II_WITH_SLEPC=ON -DSLEPC_DIR=$SLEPC_DIR -DSLEPC_ARCH=$PETSC_ARCH  \
        -DDEAL_II_WITH_THREADS=ON ..

    # Then install and test
    make -j7 install
    #make test
