### NVIDIA HPC SDK
# F90 = nvfortran -cuda -O2
# F90 = mpif90 -cuda -O2
F90 = mpif90 -cuda -O2 -mcmodel=medium
# F90 = mpif90 -cuda -O2 -mcmodel=medium -gpu=cuda11.6,cc80 -tp=zen2


OBJ   = u2v.o fileio.o \
	bc.o rk.o flux_solver.o limiter.o glm_ss2.o fastest_speed.o
OBJ_S = main.o model.o
OBJ_P = parallel.o mainp.o modelp.o mpiio.o

all : run runp

# serial version
run : $(OBJ) $(OBJ_S)
	$(F90) -o a.out $(OBJ) $(OBJ_S)

model.o: model.f90 param.h
	$(F90) $(FFLAGS) -c $<

fileio.o: fileio.f90 param.h
	$(F90) $(FFLAGS) -c $<

# parallel version
runp : $(OBJ) $(OBJ_P)
	$(F90) -o ap.out $(OBJ) $(OBJ_P)

modelp.o: modelp.f90 param.h
	$(F90) $(FFLAGS) -c $<

mpiio.o: mpiio.f90 param.h
	$(F90) $(FFLAGS) -c $<

%.o: %.cuf param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f *.out *.o *.mod data/*.dat
	rm -rf __pycache__/
