#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
#include<complex.h>

#define pi acos(-1)

#define diff_no 0.25

// Grid parameters navier stokes
    #define H 1.2                      //Height of the domain
    #define L 1.2                      //Length of the domain
    #define Radius 0.05                //Radius of the circle
    #define Ny 201                     //Number of points in Heigt
    #define Nx 201                     //Number of points in Length
    #define CFL 1                      //Time step 
    #define acf 10                     //Artificial Compressibility Factor

// Flow parameters
#define rho 0.4                     //Density
#define myu 1e-3                    //viscosity
#define U_ini 1                     //Inlet u elocity
#define V_ini 0                     //Inlet v velocity
#define P_out 0                     //Outlet Pressure
#define Psi0 0                      //Stream function value on cylinder 

// Common grid parameters
#define N_division 12               //Number of points to divide the immersed cells for non linear analysis
#define dist_factor 0.9               //Multiplied with the diagonal length of cell to see ...
#define d_cut_tol 1e-6

// Solver parameters
#define w 1.0                       //Parameter for Successive over-relaxation
#define residual_crit 1e-6          //Convergence criterion
                                    //how far in fluid domain to take points for transfinite interpolation
// Some useful simple functions
#define swap(x,y,z) z=x;x=y;y=z;    // To interchange x and y values
#define max(x,y) ((x-y)>0.0?x:y)    // To find minimum of two inputs
#define min(x,y)((x-y)<0.0?x:y)     // To find minimum of two inputs

/*   STRUCTURE DECLARATIONS   */

struct points
{
    int state;  
    double x, y;
    double u, v, Psi, p;
};

struct cells
{
    int state;                          // State of cell: Fluid, Solid or Immersed
    int parent;                         
    int child;
    int vertex[4];                      // 4 Vertices of each cell
    double x, y;                        // Location of cell center
    double dx, dy;                      // Dimensions of cell
    int face[4][2];                     // Keeps record of the 4 faces and their nodes
    double vol_fr;                      // Percentage of the cell volume inside the fluid domain
    double Psi[2];                      // Stream Function 
    double u[2], v[2], p[2], Cp[2];     // u, v velocity, pressure and Cp for old and new time step
    double c;                           // Speed of sound
    double U[3][2];                     // Flow variables [p, u, v]
    double F_e[3], F_w[3];              // Convective fluxes on east and west faces
    double Fv_e[3], Fv_w[3];            // Diffusive fluxes on east and west faces
    double G_n[3], G_s[3];              // Convective fluxes on north and south faces
    double Gv_n[3], Gv_s[3];            // Diffusive fluxes on north and south faces

    int ib_cell_no;                     // Some cells of the Cartesian cells are immersed cells. A separate structure for immersed cells has been made which stores extra information required to perform flow computations on immersed cells. The ib_cell_no tells the location/number of the IB cell in that structure.
};

struct ib_cells{
    int position;                       // Location inside or outside the immersed boundary
    int cell_no;                        // Links IB cell number back to Global cell number
    int cell_vol_type;                  // Types of ib cells based on boundary cut
    int vertex[4];                      // Cell vertices
    double x, y;                        // Cell center
    double dx, dy;                      // Cell dimensions
    int face[4][2];                     // Keeps record of the 4 faces and their nodes
    double vol_fr;                      // Fraction of cell volume inside fluid domain
    double alpha[4];                    // The fraction of face length in the fluid domain
    double Psi[2];                      // Stream Function 
    double u[2], v[2], p[2], Cp[2];     // Flow variables on the immersed boundary surface
    double c;                           // Speed of sound
    double U[3][2];                     // Flow variables [p, u, v] on cell center
    double F_e[3], F_w[3];              // Convective fluxes on east and west faces
    double Fv_e[3], Fv_w[3];            // Diffusive fluxes on east and west faces
    double G_n[3], G_s[3];              // Convective fluxes on north and south faces
    double Gv_n[3], Gv_s[3];            // Diffusive fluxes on north and south faces
    double wall_flux_convective[3];     // Convective flux on the boundary/wall
    double wall_flux_diffusive[3];      // Diffusive flux on the boundary/wall
    double x_cut[2], y_cut[2];          // Boundary cuts immersed cells are 2 points which are stored here
    double theta_cut[2];                // Stores the angles of the cut points
    double slope_x_cut[2], slope_y_cut[2];      // Stores the slopes of the cut points
    double d_cut;                               // Approximated linear length of boundary face inside an immersed cell
    double d_cut_NL;                            // Exact length of boundary face inside an immersed cell
    double dPsi_dx;                             
    int BI_cell_1[4];
    int BI_cell_2[4];
    int BI_cell_3[4];
    int BI_cell_4[4];
    int cut_face[2];
    double theta;                       // Angle on the mid point of the linear boundary face
    double x_m, y_m;                    // Mid point of the linear boundary face
    double n_x, n_y;                    // Components of unit normal pointing outward from the immersed boundary extended from the mid_point x_m, y_m  
    double d1, d2;                      // Distances of the 2 normal probes dropped from the mid point of the linear boundary face
    double x_1, y_1;                    // }
    double x_2, y_2;                    // }   Points in fluid domain for transfinite interpolation
    double x_3, y_3;                    // }
    double x_4, y_4;                    // }
    double Psi_1, Psi_2, Psi_3, Psi_4;
    double U0_1, U0_2, U0_3, U0_4;
    double U1_1, U1_2, U1_3, U1_4;
    double U2_1, U2_2, U2_3, U2_4;
    double dx1_1, dx2_1, dy1_1, dy2_1;
    double dx1_2, dx2_2, dy1_2, dy2_2;
    double dx1_3, dx2_3, dy1_3, dy2_3;
    double dx1_4, dx2_4, dy1_4, dy2_4;
    double dx1[4], dx2[4], dy1[4], dy2[4];
    int BI_flux_cell[4][4];
    int N_local_division;
    double slope_push;
    double eta_x[N_division], eta_y[N_division];
    double del_x[N_division], del_y[N_division];
};

struct faces
{
    int position;
    int node[2];
};


// Variables Used
extern struct points *point;                           //Global
extern struct points *point_Body;                      //Global
extern struct cells *cell;                             //Global
extern struct ib_cells *cell_ib;                       //Global
extern struct faces *face_Body;                        //Global

extern int Total_Grid_Points;                          //Global
extern int Total_Body_Points;                          //Global
extern int N_fl, N_sl, N_ib;                           //Global
extern int cell_Nx, cell_Ny, Total_cells;              //Global
extern double dx, dy, r, dmin;                         //Global

extern double xo, yo;                                  //Global
extern double xmax, xmin, ymax, ymin;                  //Global

extern double dt;
extern double Re;
