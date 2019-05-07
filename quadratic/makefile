all:
	gfortran -g -Wall -ffpe-trap=zero -fcheck=bounds quadratic.f90 -lopenblas -o quadratic
	./quadratic

clean:
	rm quadratic
