### NVIDIA
F90 = nvfortran -Mcuda
### NIFS
# F90 = mpif90 -Mcuda
# F90 = mpif90 -Mcuda -mcmodel=medium


OBJS = model.o u2v.o fileio.o \
	bc.o rk.o flux_solver.o limiter.o flux_glm.o glm_ss2.o fastest_speed.o

all : run

# serial version
run : main.o $(OBJS)
	$(F90) -o a.out main.o $(OBJS)

model.o: model.f90 param.h
	$(F90) $(FFLAGS) -c $<

fileio.o: fileio.f90 param.h
	$(F90) $(FFLAGS) -c $<

%.o: %.cuf param.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f *.out *.o *.mod data/*.dat
