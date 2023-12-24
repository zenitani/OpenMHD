### gfortran
F90 = gfortran -O2 -Wall
### Intel fortran ("-heap-arrays" option is useful)
# F90 = ifx -heap-arrays
### AMD AOCC fortran
# F90 = flang -O2

OBJS = model.o u2v.o output.o limiter.o flux_solver.o

run : main.o $(OBJS)
	$(F90) -o a.out main.o $(OBJS)

%.o: %.f90 param_rela.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f a.out *.o data/x-*.dat
