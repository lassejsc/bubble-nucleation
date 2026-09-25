# bubble-nucleation

A simple simulation used for determining the critical fractional volume of bubbles for percolation to occur in a cosmological phase transition. Details in the thesis this code was used in, which can be found on [Helda](127.0.0.1).

The condenced details of the simulation are as follows. The bubble placing algorithm is based on the simulation methods by Enqvist et al. (<https://doi.org/10.1103/PhysRevD.45.3415>). The percolation checking is based on finding cluster of bubbles with Depth-first search algorithm and checking whether the found cluster wraps around the periodic volume.

~**!!! Currently this repo is work in progress, the documentation and clean up of the code is incomplete slightly. !!!**~

~I really haven't had the time nor motivation to do this on my free time.~

Changes to come, fixed the minor boundary issue and started cleaning the code

The documentation in the constant nucleation rate with expanding vacuum dominated background code is quite bad, but it almost exactly the same as the better documented fast transition case, with very minor modification to the next nucleation time for example. As such I recommend looking at the code for the fast transition.

## Usage

### Example

Example of running the simulation for equal sized spheres and computing the critical fractional volume from the finite scaling.

```
python run.py --transition="equal_size" -j 10 -n 500 -b 270 300 350 400 450 500 550 600 650 700 -r 10
```

This should produce two plots, one for the step function fits and one for finite scaling fit in the range \[270,700\]. This will do 500 repeats on 10 threads so total of 5000 trials for each boundary size with sphere radius set to 10. The finite scaling should give a result of ~0.289ish.
