CC = gcc 
CFLAGS = -Wall -g -fopenmp -lm -Ofast

all: is_point_in_polygon.o ibm_functions.o generate_cartesian_mesh.o place_body.o ibm_mesh.o ip_mesh.o initialize.o calculate_time_step.o calculate_wall_flux.o calculate_convective_flux.o calculate_diffusive_flux.o cell_calculations.o calculate_residual.o extract_flow_info.o sort_data.o print_data.o release_memory.o ibm.o
	$(CC) -o ibm is_point_in_polygon.o ibm_functions.o generate_cartesian_mesh.o place_body.o ibm_mesh.o ip_mesh.o initialize.o calculate_time_step.o calculate_wall_flux.o calculate_convective_flux.o calculate_diffusive_flux.o cell_calculations.o calculate_residual.o extract_flow_info.o sort_data.o print_data.o release_memory.o ibm.o $(CFLAGS)

clean:
	rm -f ibm
	rm -f *.o 
	rm -f *.dat
	rm -f *.*~
	rm -rf *.dSYM
