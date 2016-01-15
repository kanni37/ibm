/*     FUNCTION PROTOTYPES      */
#include<stdio.h>

struct points;

// File name    -   ibm_functions.c
int move(double, double, struct points*, int);
int gcell(int, int);
int gpoint(int, int);
int comp_ib_cells(const void*, const void*);

// File name    -   is_point_in_polygon.c
int is_point_in_polygon(double, double);

// File name    -   generate_cartesian_mesh.c
void generate_cartesian_mesh();

// File name    -   place_body.c
void place_body(int);
void place_circle();
void place_square();
void place_nothing();
void place_airfoil();

// File name    -   ibm_mesh.c
void ibm_mesh(int);
void ibm_FD_mesh();
void ibm_FV_mesh();

// File name    -   ip_mesh.c
void ip_mesh(int);
void johansen_colella_mesh();
void johansen_colella_mesh_modified();
void transfinite_ip_mesh();

// File name    -    initialize.c
void initialize(int, int);
void initialize_potential_flow(int);
void initialize_navier_stokes(int);

// File name    -   cell_calculations.c
void update_old_values(int);
void update_old_values_potential_flow();
void update_old_values_navier_stokes();
void cell_calculations(int);
void cell_calculations_potential_flow();
void cell_calculations_navier_stokes();

// File name    -   calculate_time_step.c
void calculate_sound_speed();
void calculate_time_step(int );
void calculate_time_step_potential_flow();
void calculate_time_step_navier_stokes();

// File name    -   calculate_convective_flux.c
void calculate_convective_flux(int, int, int, int, int);
void calculate_convective_flux_potential_flow();
void calculate_convective_flux_navier_stokes(int, int, int, int);

// File name    -   calculate_diffusive_flux.c
void calculate_diffusive_flux(int, int, int, int, int);
void calculate_diffusive_flux_potential_flow();
void calculate_diffusive_flux_navier_stokes(int, int, int, int);

// File name    -   calculate_wall_flux.c.c
void calculate_wall_flux(int, int);
void wall_flux_potential_flow_transfinite();
void wall_flux_potential_flow_johansen_colella();
void wall_flux_navier_stokes_transfinite();
void wall_flux_navier_stokes_johansen_colella();

// File name    -   calculate_residual.c
double calculate_residual(int, FILE *, int);
double calculate_residual_potential_flow();
double calculate_p_residual_navier_stokes();
double calculate_u_residual_navier_stokes();
double calculate_v_residual_navier_stokes();

// File name    -   extract_flow_info.c
void extract_ib_data_transfinite();
void extract_ib_data_johansen_colella();
void extract_cell_data();

// File name    -   sort_data.c
void sort_data();

// File name    -   print_data.c
void print_cell_data(int, int);
void print_node_data(int);
void print_ib_data();
void print_cell_data_potential_flow();
void print_cell_data_navier_stokes();
void print_node_data_potential_flow();
void print_node_data_navier_stokes();

// File name    -   release_memory.c
void release_memory();

