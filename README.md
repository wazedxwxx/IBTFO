# <img src="https://github.com/wazedxwxx/IBTFO/blob/main/Share/Logo/logo.svg" width="64" valign="middle" alt="Spack"/> IBTFO 
[![License](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/wazedxwxx/IBTFO/)
[![License](https://img.shields.io/github/license/wazedxwxx/IBTFO)](https://opensource.org/licenses/MIT)

**IBTFO** is **discrete forcing** **I**mmersed **B**oundary method fast **T**est **F**acility based **O**penAcc. It is designed for development and validation of immersed boundary methods. IBTFO uses Built-in Euler equation solver and shared memory computation based OpenACC.
The vision of IBTFO is not to provide a robust and powerful solver, but to focus on the efficient and fast development and validation of advanced immersed boundary methods!

>Note that IBTFO currently only supports Direct BC Imposition. Such methods include ghost fluid method (GFM) and cut cell method (CCM) and their extensions.

## Installation
Even though it has simple built-in utilities as a fast test facility, `IBTFO` still relies on a few libraries.
+ `GNU Make`-determines which files need to be compiled
+ `Nvidia HPC SDK`-heterogeneous parallel computing Using OpenACC. Nvidia acquired The Portland Group, Inc, the "PGI Compilers and Tools" technology is a part of the Nvidia HPC SDK product available as a free download from Nvidia. 
+ `libpng`-libpng is the official PNG reference library

We recommend using [Spack](https://github.com/spack/spack)  to install and manage dependencies. It can easily install the required libraries and manage the version of each library.
If you decide to use spack to manage dependencies, our recommended do:
```
spack install libpng@1.6.37 nvhpc@21.9
```
## Getting Started
To build IBTFO and run a sample 2D SOD problem:
1. Import related dependencies:
```
spack load libpng@1.6.37 nvhpc@21.9
```
2. CLone the repository and enter the SOD test case
```
git clone git@github.com:wazedxwxx/IBTFO.git
cd IBTFO-main/Test/SOD/
```
3. Complie
```
make
```
4. Run
```
./REGULAR_SOD.EX sod.inp
```

### NOTE
>+  If ARCH=CPU in the compilation options, the number of cores can be selected by setting ACC_NUM_CORES
>```
>export ACC_NUM_CORES=8 
>```
>After importing the above environment variables, IBTFO will use 8 cores for calculation.
>+  If ARCH=GPU in the compilation options, the number of GPU can be selected by setting CUDA_VISIBLE_DEVICES. 
>```
> CUDA_VISIBLE_DEVICES="0" ./INCLINE_SOD.EX sod.inp 
>```
>This will only call GPU 0 for computation, otherwise all GPUs will be used.

## Tutorial
The tutorial of IBTFO can be found here.


## Contributing
If you want to contribute to the development of IBTFO, have a look at the contribution guidelines. We discourage the development of advanced advection schemes if not required for your immersed boundary method. We encourage the development and validation of your advanced immersed boundary methods. If your method is publicly published, we'd be happy to couple it into IBTFO. If you have any questions, you can discuss in the issue.

## License
**IBTFO** (Copyright (c) 2020-2022) is an open-source package made available under the MIT License.
All new contributions must be made underthe MIT license.

## Credits
The code has been developed by [Science and Technology on Scramjet Laboratory](https://english.nudt.edu.cn).
