FC := gfortran
FCFLAGS := -g -Wall -O3

all:
	$(FC) $(FCFLAGS) -o octree octree.f90 $(LIBS)
	./octree

clean:
	rm octree
