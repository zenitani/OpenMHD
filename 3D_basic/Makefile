### gfortran
# F90 = gfortran -O2 -Wall
F90 = mpif90 -Wall -O2
### Intel fortran ("-heap-arrays" option is useful, but use 'ulimit -s unlimited' for OpenMP)
# F90 = ifx -ipo -heap-arrays
# F90 = mpiifx -debug all -heap-arrays
# F90 = mpiifx -g
# F90 = mpiifx -ipo -heap-arrays
# F90 = mpiifx -ipo -qopenmp
### AMD AOCC fortran
# F90 = flang -O2
### CfCA/XC, KDK-A system
# F90 = ftn -hnoomp   # OpenMP OFF (Cray compiler)
# F90 = ftn -homp     # OpenMP ON  (Cray compiler)
# F90 = ftn -ipo          # OpenMP OFF (Intel compiler)
# F90 = ftn -ipo -qopenmp # OpenMP ON  (Intel compiler)
### JAXA/FX100
# F90 = f90sjx -Umpi -Uflat

OBJ   = set_dt.o u2v.o bc.o fileio.o rk.o flux_solver.o limiter.o glm_ss2.o
OBJ_S = main.o model.o
OBJ_P = parallel.o mainp.o modelp.o mpiio.o

all: run runp

# serial version
run : $(OBJ) $(OBJ_S)
	$(F90) -o a.out $(OBJ) $(OBJ_S)

# parallel version
runp : $(OBJ) $(OBJ_P)
	$(F90) -o ap.out $(OBJ) $(OBJ_P)


%.o: %.f90 param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f *.out *.o *.mod data/*.dat
	rm -rf __pycache__/
