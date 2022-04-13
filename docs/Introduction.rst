.. highlight:: rst

.. Warning:: This documentation is a work in progress. It is reasonably complete and hopefully useful but is likely error prone and in places misleading.


Introduction
============

**IBTFO** is discrete forcing Immersed Boundary method fast Test Facility based OpenAcc. It is designed for development and validation of immersed boundary methods. IBTFO uses Built-in Euler equation solver and shared memory computation based `OpenAcc <https://www.openacc.org/>`_ . The vision of IBTFO is not to provide a robust and powerful solver, but to focus on the efficient and fast development and validation of advanced immersed boundary methods!

Immersed Boundary method
------------
The immersed boundary method (IBM) proposed by `Peskin <https://doi.org/10.1017/S0962492902000077>`_  is employed to handle elastic boundaries for simulating blood flow in the heart.  `Mittal's <10.1146/annurev.fluid.37.061903.175743>`_  review refers to this method as the ”continuous forcing method”. The forcing is incorporated into the continuous equations before discretization. This method is suitable for the case where the boundary is elastic, and the Dirac delta function used for the elastic boundary does not fit well with the rigid boundary in general. 

The discrete forcing approach is proposed to overcome the problem that the elastic force source term cannot be given precisely in the Navier-Stokes (NS) governing equations. The discrete forcing approach is implemented by modifying the computational stencil localized at the difference scheme boundary and imposing the boundary conditions on the immersed boundary. The discrete forcing approach can reconstruct the boundary accurately. Previous research has established two approaches to modify the computational stencil: cut-cell methods (**CCM**) and ghost fluid methods (**GFM**).


Dependencies
------------

Even though it has simple built-in utilities as a fast test facility, IBTFO still relies on a few libraries.

``GNU Make``-determines which files need to be compiled

``Nvidia HPC SDK``-heterogeneous parallel computing Using OpenACC. Nvidia acquired The Portland Group, Inc, the "PGI Compilers and Tools" technology is a part of the Nvidia HPC SDK product available as a free download from Nvidia.

``libpng``-libpng is the official PNG reference library

We recommend using  `Spack <https://github.com/spack/spack/>`_  to install and manage dependencies. It can easily install the required libraries and manage the version of each library. If you decide to use spack to manage dependencies, our recommended do:
::

	spack install libpng@1.6.37 nvhpc@21.9
