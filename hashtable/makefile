FC := gfortran
FCFLAGS := -g -fcheck=bounds

all: hash

hash: hash.f90
	$(FC) $(FCFLAGS) -o hash hash.f90 $(LIBS)
