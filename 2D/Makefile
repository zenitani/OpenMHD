### gfortran
# F90 = gfortran -O2 -m64 -Wall
# F90 = /opt/openmpi64g/bin/mpif90 -Wall -O2
### intel fortran ("-heap-arrays" option is useful, but use 'ulimit -s unlimited' for OpenMP)
# F90 = ifort -heap-arrays
# F90 = /opt/openmpi64i/bin/mpif90 -debug all -heap-arrays
# F90 = /opt/openmpi64i/bin/mpif90 -g -m64
F90 = /opt/openmpi64i/bin/mpif90 -heap-arrays
# F90 = /opt/openmpi64i/bin/mpif90 -openmp
### XC
# F90 = ftn
### JAXA/FX100
# F90 = f90sjx -Umpi

OBJS = model.o set_dt.o u2v.o input.o output.o rk21.o rk22.o step1.o step2.o bc.o \
	v2u.o v2f.o v2g.o limiter_f.o limiter_g.o \
	hlld_f.o hlld_g.o hllc_f.o hllc_g.o hll_f.o hll_g.o llf_f.o llf_g.o \
	glm_f.o glm_g.o glm_ss.o
OBJP = modelp.o mpibc.o mpibc_vlvr_f.o mpiinput.o mpioutput.o

all: run runp

run : main.o $(OBJS)
	$(F90) -o a.out main.o $(OBJS)

runp : mainp.o $(OBJS) $(OBJP)
	$(F90) -o ap.out mainp.o $(OBJS) $(OBJP)


%.o: %.f90 param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f *.out *.o data/*.dat
