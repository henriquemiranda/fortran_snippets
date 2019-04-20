FC := gfortran
FCFLAGS := -g -Wall

all: octree

octree: octree.f90
	$(FC) $(FCFLAGS) -o octree octree.f90 $(LIBS)

clean:
	rm octree
