# Circle median filter toolbox (CMF)
 
This toolbox contains a fast algorithm for median filtering of signals and images with values 
    on the unit circle, for example phase or orientation data.
It is a C++ reference implementation of the algorithms described in the paper

Martin Storath, Andreas Weinmann,
[Fast median filtering for phase or orientation data](https://hci.iwr.uni-heidelberg.de/sites/default/files/profiles/mstorath/files/storath2017fast.pdf), 
IEEE Transactions on Pattern Analysis and Machine Intelligence, 2017 (in press)

### Example


![alt tag](https://hci.iwr.uni-heidelberg.de/sites/default/files/publications/teaserimages/1908951751/mediancircularrevision_teaser_small.png)

*Left:* A circle-valued image, i.e. every pixel takes its value on the unit circle (or in angular representation a value in (-pi, pi]). The values are visualized as hue component in the HSV color space.
*Right:* Action of the circle-median filter using a filter mask of size 7 × 7. 


### Contents
- demos:     some demos, self explanatory (implemented in Matlab)
- auxiliary: some useful helper functions (implemented in Matlab)
- filters:   the fast algorithms for median filtering of circle valued data 
(implemented in C++ with Matlab wrappers)

### Installation and usage

- Run CMF_install.m in the Matlab console and follow the demos 

