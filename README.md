# func-approx-fwam2022
This repository contains demo codes for function approximation in one dimension.

To download this repository including submodules (chebfun and baobzi)  run

`git clone --recurse-submodules https://github.com/flatironinstitute/gp-shootout.git`

Run startup.m to include chebfun in the MATLAB path

To install baobzi on a linux machine, follow instructions on the
baobzi[https://github.com/flatironinstitute/baobzi.git]. Baobzi depends
on gcc, and cmake.

To install baobzi on a mac, follow the same instructions. Additionally
to run the code in matlab, you will need to copy the libbaobzi.dylib to
`/usr/local/lib` and add `${BAOBZI_ROOT}/share/baobzi/matlab` to
the matlab path as well.



