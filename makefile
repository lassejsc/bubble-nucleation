

FC=gfortran
FFLAGS=-c -O3 -flto
OPTFLAGS=-O3 -flto


fast_transition_sim: percolation_tools.o
	$(FC) $(OPTFLAGS) percolation_tools.o ./Fast\ Transition/fast_transition_sim.f90 -o fast_transition_sim.o -fopenmp

slow_transition_sim: percolation_tools.o
	$(FC) $(OPTFLAGS) percolation_tools.o ./Slow\ Transition/slow_transition_sim.f90 -o slow_transition_sim.o -fopenmp

percolation_tools.o: percolation_tools.f90
	$(FC) $(FFLAGS) percolation_tools.f90
	
clean:
	rm *.o *.mod
