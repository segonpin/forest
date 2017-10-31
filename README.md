What is FOREST?
===================

FOREST is an acronym for Framework Oriented for REactor SimulaTions. 
It is built on top of deal.II, a C++ program library "targeted
at the computational solution of partial differential equations using
adaptive finite elements". It makes extensive use of Object Oriented
features of C++, as well as template classes and functions.


Why FOREST?
===================

There exist many software targeted at reactor calculations. Many
of these tools are in-house tools with no available source code
(industrial software), or software developed specifically for 
the specific problem in mind (Dragon, OpenMOC). There also exist tools
intended to solve the neutron transport equation using other
Finite Element Libraries (RattleSNake, EVENT, ...), which are closer
to what this project is focused on.

Here we are interested into releasing an Open-Source library that 
can be used to perform research with the more advanced numerical 
techniques for Finite Element Methods for reactor calculations. 
For this purpose, we believe that building the library on top of
a stable and growing library, which continues developing and 
improving adds value to the present tool, in order get the latest
improvements without big efforts or modifications. This also includes
PETSc library for numerical linear algebra, SLEPc for eigenvalue
problems, and TRILINOS for nonlinear solvers. 

This library present an interface for easy usage that for the moment
account for BWR and PWR geometries, but is intended to be highly 
modular and easy to extend to more general geometries. The geometry module
is mainly intended to process the input and generate the mesh, 
that can be at the same time generated with other external tools (gmsh)
thus allowing as general geometries as possible for Finite Element
Methods. One particular restriction of the library is that the simplices
considered are only quadrangular elements, while this does not represent
are strong limitation due to the possibility of represent very different
geometries using this simplices as base elements. 

Easy Install
===================

Configure compile and install FOREST with:

```bash
$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=/path/where/deal.II/should/be/installed/to ..
$ make install    (alternatively $ make -j<N> install)
```

A detailed ReadME can be found somewhere (not yet).

Look at https://dealii.org/developer/users/cmake.html#advanced 
for more tips about how to automatize the instalation with cmake.

Tutorial
===================

The tutorial steps are located under ./examples
Information about the tutorial steps can be found at
./doc/doxygen/tutorial/index.html or at http://www.dealii.org/.

License
===================

@todo I do not know what to put here. Think about that.

Further information
===================


@todo I do not know what to put here. Think about that.
For further information have a look at "local documentation.html" or at
"global documentation.html".
  
Contact
==============

sebastian.gonzalez-pintor[at]chalmers[dot]se

Dependencies
==============

Here we have PETSc (version 3.6), SLEPc (version 3.6) and 
deal.II (version 8.4). They should be installed
in that order, to satisfy the dependencies. The way to install each
one of them is as follows. (insteresting [link](http://nicklj.com/?p=669) )


PESTc
-----

In order to download PESTc we write

~~~~.console
git clone -b maint https://bitbucket.org/petsc/petsc petsc
~~~~

In order to install PESTc we write

~~~~.console
./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-fblaslapack --download-hypre=1 --download-mpich
~~~~

which should take care of download and installing 
mpi, download and installing hypre (necessary in
dealii.8-4 and later) and configuring the installation
of petsc, or

~~~~.console
./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-fblaslapack --with-mpi=0
~~~~

for installing without mpi. Then we install it by typing

~~~~.console
make all test
~~~~

SLEPc
-----

SLEPc works in a very similar way to that of PETSc.
The main point here is to tell SLEPc where is PETSc, i.e.,

    Set the environment variable SLEPC_DIR to the full path
    of the SLEPc home directory, for example,
    `export SLEPC_DIR=/home/username/slepc`
    In addition to this variable, PETSC_DIR and PETSC_ARCH
    must also be set appropriately.

We can download SLEPc in a very similar way to PETSc, by typing

~~~~.console
git clone -b maint https://bitbucket.org/slepc/slepc slepc
~~~~

Then we go inside the SLEPc folder and we type:

~~~~.console
./configure
~~~~

The last step is to install by typing

~~~~.console
make all test
~~~~

deal.II
-------

For installing deal.II, we start downloading the library
to the folder that we choose, which in my case is

~~~~.console
cd /home/segonpin/lib
git clone -b dealii-8.4 https://github.com/dealii/dealii.git dealii
cd dealii
mkdir build
cd build
~~~~

Now we prepare the directory for the installation and
export it to the path,

~~~~.console
mkdir ../../dealii.bin
export DEAL_II_DIR=../../dealii.bin
~~~~

then I run the bash script `sebas_buildingdealii.sh`
by typing

~~~~.console
sh sebas_buildingdealii.sh
~~~~

where the bash script contains the following lines

~~~~{.bash}
#/bin/bash

# here I am using my enviromental variable DEAL_II_DIR
# as the installing directory. It can be any other path.
echo "$DEAL_II_DIR"

# prepare the cache for the cmake
cmake -DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR \
    -DDEAL_II_WITH_FUNCTIONPARSER=ON \
    -DDEAL_II_COMPONENT_DOCUMENTATION=ON \
    -DMPI_CXX_COMPILER=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx \
    -DDEAL_II_WITH_PETSC=ON -DPETSC_DIR=$PETSC_DIR -DPETSC_ARCH=$PESTC_ARCH \
    -DDEAL_II_WITH_SLEPC=ON -DSLEPC_DIR=$SLEPC_DIR -DSLEPC_ARCH=$PETSC_ARCH \
    -DDEAL_II_WITH_THREADS=ON ..

# Then install and test
make -j7 install
#make test
~~~~