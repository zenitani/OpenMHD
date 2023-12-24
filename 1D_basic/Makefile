### gfortran
F90 = gfortran -O2 -Wall
### Intel fortran ("-heap-arrays" option is useful)
# F90 = ifx -heap-arrays
### AMD AOCC fortran
# F90 = flang -O2

OBJS = main.o model.o set_dt.o u2v.o output.o flux_solver.o limiter.o

all: run

run : $(OBJS)
	$(F90) -o a.out $(OBJS)

%.o: %.f90 param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f a.out *.o data/*.dat
