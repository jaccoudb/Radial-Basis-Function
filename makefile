FC = gfortran
FFLAGS = -Wall -ffree-line-length-none
FLIB = -llapack -lblas
FOPT = -O3

all: rbf_program.o runcode
rbf_program.o: rbf_program.f90
	$(FC) $(FOPT) rbf_program.f90 -c
runcode: rbf_program.o
	$(FC) $(FOPT) rbf_program.o -o runcode $(FLIB)
clean:
	rm *.o *.mod run*
	
