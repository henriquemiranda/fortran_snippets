FC := gfortran
FCFLAGS := -g -Wall -O3

all:
	$(FC) $(FCFLAGS) -o opt opt.f90 $(LIBS)
	./opt

clean:
	rm opt
