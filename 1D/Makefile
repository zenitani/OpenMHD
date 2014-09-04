### gfortran
# F90 = gfortran -O2 -m64 -Wall
### intel fortran ("-heap-arrays" option is useful)
F90 = ifort -heap-arrays

OBJS = model.o set_dt.o u2v.o output.o \
	v2u.o v2f.o v2g.o limiter_f.o limiter_g.o \
	hlld_f.o hlld_g.o hllc_f.o hllc_g.o hll_f.o hll_g.o llf_f.o llf_g.o

all: 1d_x 1d_y

1d_x : 1d_x.o $(OBJS)
	$(F90) -o 1d_x 1d_x.o $(OBJS)

1d_y : 1d_y.o $(OBJS)
	$(F90) -o 1d_y 1d_y.o $(OBJS)

%.o: %.f90 param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f 1d_x 1d_y *.o data/*.dat
