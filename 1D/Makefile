### gfortran
F90 = gfortran -O2 -m64 -Wall
### intel fortran ("-heap-arrays" option is useful)
#F90 = ifort -heap-arrays

OBJS = model.o set_dt.o u2v.o output.o flux_solver.o limiter.o

all: 1d_x 1d_y

1d_x : 1d_x.o $(OBJS)
	$(F90) -o 1d_x 1d_x.o $(OBJS)

1d_y : 1d_y.o $(OBJS)
	$(F90) -o 1d_y 1d_y.o $(OBJS)

%.o: %.f90 param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f 1d_x 1d_y *.o data/*.dat
