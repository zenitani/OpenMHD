#F90 = /opt/intel/Compiler/11.1/067/bin/intel64/ifort
F90 = gfortran -O2 -m64

OBJS = model.o u2v.o output.o limiter.o flux_solver.o

run : main.o $(OBJS)
	$(F90) -o a.out main.o $(OBJS)

%.o: %.f90 param_rela.h
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f a.out *.o data/x-*.dat
