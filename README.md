# Circle median filter toolbox (CMF)
 
This toolbox contains a fast algorithm for median filtering of signals and images with values 
    on the unit circle, for example phase or orientation data.
It is a C++ reference implementation of the algorithms described in the paper

Martin Storath, Andreas Weinmann,
"Fast median filtering for phase or orientation data",
IEEE Transactions on Pattern Analysis and Machine Intelligence, 2017

Contents:
- demos:     some demos, self explanatory (implemented in Matlab)
- auxiliary: some useful helper functions (implemented in Matlab)
- filters:   the fast algorithms for median filtering of circle valued data 
(implemented in C++ with Matlab wrappers)

Installation and usage:

- Run CMF_install.m in the Matlab console and follow the demos 
