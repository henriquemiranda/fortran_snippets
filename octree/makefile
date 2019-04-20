FC := gfortran
FCFLAGS := -g

all: octree

octree: octree.f90
	$(FC) $(FCFLAGS) -o octree octree.f90 $(LIBS)

clean:
	rm octree
