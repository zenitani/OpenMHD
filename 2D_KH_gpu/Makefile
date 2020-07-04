### PGI fortran
# F90 = pgfortran -Mcuda=cc35  # CfCA
# F90 = pgfortran -Mcuda -O2 -pg
# F90 = pgfortran -Mcuda=cc70 -O2
F90 = pgfortran -Mcuda -O2


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
