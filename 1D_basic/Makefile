### gfortran
F90 = gfortran -O2 -Wall
### Intel fortran ("-heap-arrays" option is useful)
# F90 = ifx -heap-arrays
### AMD AOCC fortran
# F90 = flang -O2

OBJ = main.o model.o set_dt.o u2v.o output.o flux_solver.o limiter.o

all: run

run : $(OBJ)
	$(F90) -o a.out $(OBJ)

%.o: %.f90 param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f a.out *.o data/*.dat
