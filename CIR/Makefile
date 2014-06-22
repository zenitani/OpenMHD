#F90 = /opt/intel/Compiler/11.1/067/bin/intel64/ifort
F90 = gfortran -O2 -m64 -Wall

OBJS = model.o set_dt.o u2v.o output.o \
	v2u.o v2f.o limiter_f.o hlld_f.o hllc_f.o hll_f.o llf_f.o

all: a.out

a.out : main.o $(OBJS)
	$(F90) main.o $(OBJS)

%.o: %.f90 param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f a.out *.o data/*.dat
