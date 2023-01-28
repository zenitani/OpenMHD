### gfortran
# F90 = gfortran -O2 -Wall
F90 = mpif90 -Wall -O2
### Intel fortran ("-heap-arrays" option is useful, but use 'ulimit -s unlimited' for OpenMP)
# F90 = ifort -ipo -heap-arrays
# F90 = mpiifort -debug all -heap-arrays
# F90 = mpiifort -g
# F90 = mpiifort -ipo -heap-arrays
# F90 = mpiifort -ipo -qopenmp
### AMD AOCC fortran
# F90 = flang -O2
### CfCA/XC, KDK-A system
# F90 = ftn -hnoomp   # OpenMP OFF (Cray compiler)
# F90 = ftn -homp     # OpenMP ON  (Cray compiler)
# F90 = ftn -ipo          # OpenMP OFF (Intel compiler)
# F90 = ftn -ipo -qopenmp # OpenMP ON  (Intel compiler)
### JAXA/FX100
# F90 = f90sjx -Umpi -Uflat
### Nagoya U/FX100
# F90 = mpifrtpx

OBJS = model.o set_dt.o set_dt2.o set_eta.o u2v.o bc.o fileio.o rk.o \
	flux_solver.o limiter.o flux_resistive.o glm_ss2.o
OBJP = modelp.o mpibc.o mpiio.o

all: run runp

# serial version
run : main.o $(OBJS)
	$(F90) -o a.out main.o $(OBJS)

# parallel version
runp : parallel.o mainp.o $(OBJS) $(OBJP)
	$(F90) -o ap.out parallel.o mainp.o $(OBJS) $(OBJP)

%.o: %.f90 param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f *.out *.o *.mod data/*.dat data/*.dat.restart
	rm -rf __pycache__/
