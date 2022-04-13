.. highlight:: rst

.. Warning:: This documentation is a work in progress. It is reasonably complete and hopefully useful but is likely error prone and in places misleading.


Introduction
============

IBTFO is discrete forcing Immersed Boundary method fast Test Facility based OpenAcc. It is designed for development and validation of immersed boundary methods. IBTFO uses Built-in Euler equation solver and shared memory computation based OpenACC. The vision of IBTFO is not to provide a robust and powerful solver, but to focus on the efficient and fast development and validation of advanced immersed boundary methods!



Dependencies
------------

Even though it has simple built-in utilities as a fast test facility, IBTFO still relies on a few libraries.

``GNU Make``-determines which files need to be compiled

``Nvidia HPC SDK``-heterogeneous parallel computing Using OpenACC. Nvidia acquired The Portland Group, Inc, the "PGI Compilers and Tools" technology is a part of the Nvidia HPC SDK product available as a free download from Nvidia.

``libpng``-libpng is the official PNG reference library
We recommend using  `Spack <https://github.com/spack/spack/>`_  to install and manage dependencies. It can easily install the required libraries and manage the version of each library. If you decide to use spack to manage dependencies, our recommended do:
::

	spack install libpng@1.6.37 nvhpc@21.9
