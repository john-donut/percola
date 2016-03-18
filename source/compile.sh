gfortran -c random_functions.f90 -o random_functions.o
gfortran -fbounds-check -o program.x susc.f90 random_functions.o
#gfortran -fbounds-check -o program.x monte-carlo-percolation.f90 random_functions.o
