# the compiler: gcc for C program, define as g++ for C++
CC = gcc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall -fopenmp -Ofast -lm -I.

#Dependencies
DEPS = ibm.h fuctions.h 

#Objects
OBJ = 	is_point_in_polygon.o \
		ibm_functions.o \
		generate_cartesian_mesh.o \
		place_body.o \
		ibm_mesh.o \
		ip_mesh.o \
		initialize.o \
		calculate_time_step.o \
		calculate_wall_flux.o \
		calculate_convective_flux.o \
		calculate_diffusive_flux.o \
		cell_calculations.o \
		calculate_residual.o \
		extract_flow_info.o \
		sort_data.o \
		print_data.o \
		release_memory.o \
		ibm.o

%.o: %.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS)

ibm: $(OBJ)
		gcc -o $@ $^ $(CFLAGS)


clean:
	rm -f ibm
	rm -f *.o 
	rm -f *.dat
	rm -f *.*~
	rm -rf *.dSYM
