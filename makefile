COMPILER= gfortran 

FLAGS = 

EXEC = planetary

SRC = $(wildcard *.f90) 

OBJ = $(SRC:.f90=.o)

$(EXEC): $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ $(LAPACK)

types.o: types.f90
	$(COMPILER) $(FLAGS) -c $<

read_write.o: read_write.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

ode_solver.o: ode_solver.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

mechanics.o: mechanics.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

main.o: main.f90 types.o read_write.o ode_solver.o mechanics.o
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf $(EXEC)
