FC := gfortran
FCFLAGS := -g

all: stretchylist

stretchylist: stretchylist.f90
	$(FC) $(FCFLAGS) -o stretchylist stretchylist.f90 $(LIBS)

clean:
	rm stretchylist
