### NVIDIA
F90 = nvfortran -Mcuda
### NIFS
# F90 = mpif90 -Mcuda
# F90 = mpif90 -Mcuda -mcmodel=medium
### CfCA
# F90 = mpif90 -cuda -mcmodel=medium -gpu=cuda11.6,cc80 -tp=zen2

OBJS = model.o u2v.o fileio.o \
	bc.o rk.o flux_solver.o limiter.o glm_ss2.o fastest_speed.o
OBJP = modelp.o mpiio.o

all : run runp

# serial version
run : main.o $(OBJS)
	$(F90) -o a.out main.o $(OBJS)

model.o: model.f90 param.h
	$(F90) $(FFLAGS) -c $<

fileio.o: fileio.f90 param.h
	$(F90) $(FFLAGS) -c $<

# parallel version
runp : parallel.o mainp.o $(OBJS) $(OBJP)
	$(F90) -o ap.out parallel.o mainp.o $(OBJS) $(OBJP)

modelp.o: modelp.f90 param.h
	$(F90) $(FFLAGS) -c $<

mpiio.o: mpiio.f90 param.h
	$(F90) $(FFLAGS) -c $<

%.o: %.cuf param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f *.out *.o *.mod data/*.dat
	rm -rf __pycache__/
