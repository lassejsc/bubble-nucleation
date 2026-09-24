# bubble-nucleation

A simple simulation used for determining the critical fractional volume of bubbles for percolation to occur in a cosmological phase transition. Details in the thesis this code was used in, which can be found on [Helda](127.0.0.1).

The condenced details of the simulation are as follows. The bubble placing algorithm is based on the simulation methods by Enqvist et al. (https://doi.org/10.1103/PhysRevD.45.3415). The percolation checking is based on finding cluster of bubbles with Depth-first search algorithm and checking whether the found cluster wraps around the periodic volume.

~**!!! Currently this repo is work in progress, the documentation and clean up of the code is incomplete slightly. !!!**~
~I really haven't had the time nor motivation to do this on my free time. ~
Changes to come, fixed the minor boundary issue and started cleaning the code


The documentation in the constant nucleation rate with expanding vacuum dominated background code is quite bad, but it almost exactly the same as the better documented fast transition case, with very minor modification to the next nucleation time for example. As such I recommend looking at the code for the fast transition.

## Usage

