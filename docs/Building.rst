.. _Building:

Building
--------

IBTFO uses executables which are customized to the case in which the user intends to run. As listed earlier, specific source code files are considered input to the program and compiled into the executable itself. Therefore, IBTFO has both compile-time and run-time inputs that the user must set for their case. The compile time inputs can be verified by running:

::
   ./IBTFO --describe

which will print out the build information (including git hashes, modules, EoS, some basic information from compiled chemistry network, etc). 

IBTFO has the ability to use two build systems. First is the GNU Make build system. This build system is best for use on large computing facility machines and production runs. With the GNU Make implementation, the build system will inspect the machine and use known compiler optimizations explicit to that machine if possible. These explicit settings are kept up-to-date by the AMReX project. The second build system implemented is CMake. This is best used for developers of IBTFO and more generalized. CMake allows for building as well as easy testing and verification of IBTFO through the use of CTest which is included in CMake.

GNU Make
~~~~~~~~

Using the GNU Make build system involves first setting environment variables for the directories of the dependencies of IBTFO which are the repositories of AMReX and PelePhysics. AMReX and PelePhysics are provided as git submodules in IBTFO and can be populated by using ``git submodule init; git submodule update`` in the IBTFO repo, or before cloning by using ``git clone --recursive <IBTFO_repo>``. Although submodules of these projects are provided, they can be placed externally as long as the ``<REPO_HOME>`` environment variables for each dependency is set correctly. An example of setting the ``<REPO_HOME>`` environment variables in the user's ``.bashrc`` is shown below:

::

   export IBTFO_HOME=${HOME}/IBTFO
   export AMREX_HOME=${IBTFO_HOME}/Submodules/AMReX
   export PELE_PHYSICS_HOME=${IBTFO_HOME}/Submodules/PelePhysics


Then one edits the ``GNUMakefile`` in any of the examples in the ``Exec`` directory and uses the ``make`` command to build the executable.

CMake
~~~~~

Using CMake involves an additional configure step before using the ``make`` command. It is also expected that the user has cloned the IBTFO repo with the ``--recursive`` option or performed ``git submodule init; git submodule update`` in the IBTFO repo to populate its submodules. 

To build with CMake, a user typically creates a ``build`` directory in the project directory and in that directory the ``cmake <options> ..`` command is used to configure the project before building it. IBTFO provides an example build directory called ``Build`` with example scripts for performing the CMake configure. Once the CMake configure step is done, then the ``make`` command will build the executable.

An example CMake configure command to build IBTFO with MPI is listed below:

::

    cmake -DCMAKE_BUILD_TYPE:STRING=Release \
          -DIBTFO_ENABLE_MPI:BOOL=ON \
          -DCMAKE_CXX_COMPILER:STRING=mpicxx \
          -DCMAKE_C_COMPILER:STRING=mpicc \
          -DCMAKE_Fortran_COMPILER:STRING=mpifort \
          .. && make

Note that CMake is able to generate makefiles for the Ninja build system as well which will allow for faster building of the executable(s).
