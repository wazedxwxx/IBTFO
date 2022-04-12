# IBTFO 
[![License](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/wazedxwxx/IBTFO/)
[![License](https://img.shields.io/github/license/wazedxwxx/IBTFO)](https://opensource.org/licenses/MIT)
**IBTFO** is **I**mmersed **B**oundary method fast **T**est **F**acility based **O**penAcc. It is designed for easy development and validation of immersed boundary methods. Built-in basic Euler equation solver and shared memory computation using OpenACC.

## Installation
Even though it has as simple utilities as possible built in as a fast test facility, `IBTFO` still relies on a few libraries.
+ `GNU Make`-determines which files need to be compiled
+ `Nvidia HPC SDK`-heterogeneous parallel computing Using OpenACC. Nvidia acquired The Portland Group, Inc, the "PGI Compilers and Tools" technology is a part of the Nvidia HPC SDK product available as a free download from Nvidia. 
+ `libpng`-libpng is the official PNG reference library

We recommend using [Spack](https://github.com/spack/spack)  to install and manage dependencies. It can easily install the required libraries and manage the version of each library.
If you decide to use spack to manage dependencies, our recommended do:
```
spack install libpng@1.6.37 nvhpc@21.9
```

## Tutorial
The tutorial of CGC can be found here.


## Contributing
If you want to contribute to the development of cgc, have a look at the contribution guidelines.

## License
**IBTFO** (Copyright (c) 2020-2022) is an open-source package made available under the MIT License.
All new contributions must be made underthe MIT license.

## Credits
The code has been developed developed by [Science and Technology on Scramjet Laboratory](https://english.nudt.edu.cn).
