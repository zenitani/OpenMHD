### gfortran
F90 = gfortran -O2 -Wall
### Intel fortran ("-heap-arrays" option is useful)
# F90 = ifort -heap-arrays

OBJS = model.o set_dt.o u2v.o output.o flux_solver.o limiter.o

all: a.out

a.out : main.o $(OBJS)
	$(F90) main.o $(OBJS)

%.o: %.f90 param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f a.out *.o data/*.dat
