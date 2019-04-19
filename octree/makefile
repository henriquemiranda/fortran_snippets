FC := gfortran
FCFLAGS := -pg

all: octree

octree: octree.f90
	$(FC) $(FCFLAGS) -o octree octree.f90 $(LIBS)
