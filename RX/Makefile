### gfortran
# F90 = gfortran -O2 -m64 -Wall
# F90 = /opt/openmpi64g/bin/mpif90 -Wall -O2
### intel fortran ("-heap-arrays" option is useful)
# F90 = ifort -heap-arrays
# F90 = /opt/openmpi64i/bin/mpif90 -debug all -heap-arrays
# F90 = /opt/openmpi64i/bin/mpif90 -g -m64
# F90 = /opt/openmpi64i/bin/mpif90 -heap-arrays
### NAOJ/XC
F90 = ftn
### JAXA/FX100
# F90 = f90sjx -Umpi

OBJS = model.o set_dt.o set_dt2.o set_eta.o u2v.o input.o output.o \
	rk21.o rk22.o step1.o step2.o bc.o \
	v2u.o v2f.o v2g.o limiter_f.o limiter_g.o bc_vlvr_g.o \
	hll_resistive_f.o hll_resistive_g.o \
	hlld_resistive_f.o hlld_resistive_g.o \
	glm_f.o glm_g.o glm_ss.o
OBJP = modelp.o  mpibc.o  mpibc_vlvr_f.o  mpioutput.o
OBJP2= modelp2.o mpibc2.o mpibc_vlvr_f2.o bc_vlvr_g2.o    mpiinput.o mpioutput.o

all: run runp runp2

# serial version
run : main.o $(OBJS)
	$(F90) -o a.out main.o $(OBJS)

# parallel version
runp : mainp.o $(OBJS) $(OBJP)
	$(F90) -o ap.out mainp.o $(OBJS) $(OBJP)

# parallel version (1/4 BC for reconnection papers)
runp2 : mainp2.o $(OBJS) $(OBJP2)
	$(F90) -o ap2.out mainp2.o $(OBJS) $(OBJP2)

%.o: %.f90 param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f *.out *.o data/*.dat data/*.dat.restart
