

FC=gfortran
FFLAGS=-c -O3 -flto
OPTFLAGS=-O3 -flto

all: slow_transition_sim.o fast_transition_sim.o equal_sized_percolation.o

fast_transition_sim.o: percolation_tools.o
	$(FC) $(OPTFLAGS) percolation_tools.o ./Fast\ Transition/fast_transition_sim.f90 -o fast_transition_sim.o -fopenmp

slow_transition_sim.o: percolation_tools.o
	$(FC) $(OPTFLAGS) percolation_tools.o ./Slow\ Transition/slow_transition_sim.f90 -o slow_transition_sim.o -fopenmp


equal_sized_percolation.o: percolation_tools.o
	$(FC) $(OPTFLAGS) percolation_tools.o ./Equal-sized_spheres/equal_sized_percolation.f90 -o equal_sized_spheres.o -fopenmp

percolation_tools.o: percolation_tools.f90
	$(FC) $(FFLAGS) percolation_tools.f90


clean:
	rm *.o *.mod
