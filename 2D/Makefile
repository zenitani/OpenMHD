### gfortran
# F90 = gfortran -O2 -m64 -Wall
### intel fortran ("-heap-arrays" option is useful)
F90 = ifort -heap-arrays

OBJS = model.o set_dt.o u2v.o output.o rk21.o rk22.o step1.o step2.o bc.o \
	v2u.o v2f.o v2g.o limiter_f.o limiter_g.o \
	hlld_f.o hlld_g.o hllc_f.o hllc_g.o hll_f.o hll_g.o llf_f.o llf_g.o \
	glm_f.o glm_g.o glm_ss.o

run : main.o $(OBJS)
	$(F90) -o a.out main.o $(OBJS)

%.o: %.f90 param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f a.out *.o data/*.dat
