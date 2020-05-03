CPP = gcc -E -P
FFLAGS = -g
FC = gfortran
FFLAGS = -g -Wall -Og -fcheck=bounds --free-form
all:
	$(CPP) mastermind.F -o mastermind.f90
	$(FC) $(FFLAGS) mastermind.F -o mastermind
	./mastermind

clean:
	rm master
