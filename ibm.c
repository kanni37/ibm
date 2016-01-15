#include "ibm.h"
#include "functions.h"

// Initialize or Restart ?
//          1   :   Initialize      // Start from initialized values from inlet
//          2   :   Restart         // Start from previous data
int init_type = 1;

// Boundary Conditions for Domain
//
//                      (3)
//     _______________________________________
//    |                                       |  
//    |                                       |
//    |                  __                   |
//    |                 /  \                  |   
// (1)|                 \__/                  |(2)    
//    |                                       |    
//    |                                       |
//    |_______________________________________|
//                      (4)
//
//
//  (1) Front
//          1   :   Velocity inlet
int Front = 1;
//  (2) Back
//          1   :   Pressure Outlet
int Back = 1;
//  (3) Top
//          1   :   Farfield      
//          2   :   Wall
//          3   :   Symmetry
int Top = 1;
//  (4) Bottom  
//          1   :   Farfield     
//          2   :   Wall
//          3   :   Symmetry 
int Bottom = 1;

// Equation type
//          1   :   Potential Flow
//          2   :   Incompressible Navier Stokes
int eq_type = 2;

// Method_type 
//          1   :   Transfinite Interpolation   // (refer thesis)
//          2   :   IBM1                        // Original Embedded Boundary Method. (refer thesis)
//          3   :   IBM2                        // 2 points normal to boundary at equal distance (refer thesis).
int M_type = 2;

// Discretization_type
//          1   :   Finite Volume
//          2   :   Finite Difference
int D_type = 1;

// Body to place
//          1   :   No body
//          2   :   2D Circular Cylinder
//          3   :   2D Square Cylinder
//          4   :   2D Airfoil
int Immersed_body_type = 2;

// Variables Used
struct points *point = NULL;                    //Global        // Cartesian Grid points/nodes
struct points *point_Body = NULL;               //Global        // Points which define the geometry of the body
struct cells *cell = NULL;                      //Global        // Caresian grid cells
struct ib_cells *cell_ib = NULL;                //Global        // Immersed cells which contain the boundary
struct faces *face_Body = NULL;                 //Global        // Faces or lines which make the boundary of immersed body.

int Total_Grid_Points;                          //Global        // Total number of grid points
int Total_Body_Points;                          //Global        // Total number of body points
int N_fl, N_sl, N_ib;                           //Global        // N_fl = Number of fluid cells, N_sl = Number of solid cells, N_ib = Number of ib cells
int cell_Nx, cell_Ny, Total_cells;              //Global        // cell_Nx = Cells in x direction, cell_Ny = Cells in y direction
double dx, dy, r, dmin;                         //Global        // dx, dy = grid spacing in x and y // r = ratio of grid spacing // dmin = minimum distance of placing probe points for immersed cells

double xo, yo;                                  //Global        // Position of body in the grid
double xmax, xmin, ymax, ymin;                  //Global        // Bounding box for the body required for faster ray-casting
double dt = 1e-6;
double Re;

int main(int argc, char *argv[])
{
    double residual;
    int count_residual, count_iter, count_print;                // Counters
    struct timeval start_time, stop_time, elapsed_time;         // Timers

    gettimeofday(&start_time, NULL);    // Unix timer
    
    //------------------------------ Setting up Mesh ---------------------------

    // Generate Cartesian Mesh.
    // generate_FVM_mesh for a finite volume mesh // Grid points and cells constructed
    // generate_FDM_mesh for a fintite difference mesh // Only grid points constructed
    generate_cartesian_mesh();

    // Place Geometry over Cartesian Mesh.
    // To place any other geometry: In the place_body.c file make a new function which defines the shape of that geometry.
    place_body(Immersed_body_type);

    // Identify Cartesin grid points inside and outside the body
    // Used Ray casting algorithm (finding point in a polygon)
    // ibm_FD_mesh() //Identifies Fluid, Solid and Ghost Nodes and generates ibm mesh with all required information
    // ibm_FV_mesh() //Identifies Fluid, Solid and Immersed Cells and generates ibm mesh with all required information
    // Information stored in state of points or cells
    // state: 0 = Fluid, 1 = Solid, 2 = Immersed/Ghost
    ibm_mesh(D_type);

    // To enforce boundary condition on the immersed boundary, some information is needed
    // transfinite_ip_mesh() generates points for Immersed cells which will be used for transfinite interpolation
    // Program can be automated based on number of layers required ! or if some other scheme is to be used !
    // --------------
    // 2nd part: Finds the information required to interpolate the data on the transfinite grid points
    // Information includes: Surrounding cells; Distances for bilinear interpolation;
    ip_mesh(M_type);
    
    //-------------------------------- Flow Solver -----------------------------

    // Initilizing the flow variables
    initialize(eq_type, init_type);

    // Gauss Siedel with successive over relaxation
    count_residual = 0;
    count_iter = 0;
    count_print = 0;

    FILE *fp_residual = NULL;
    fp_residual = fopen("Residual.dat", "w");
    fprintf(fp_residual,"VARIABLES = \"Iteration Count\", \"Residual_p\", \"Residual_u\", \"Residual_v\"\n");
    
    do{

        // Calculating time step for stable time marching
        // calculate_time_step(eq_type);

        // Calculating Diffusive Fluxes 
        // Implements boundary condition within the fluxes calculated
        // Boundary condition here include: Farfield, Wall , Inlet and Outlet
        calculate_convective_flux(eq_type, Front, Back, Top, Bottom);
       
        // Calculating Convective Fluxes
        // Implements boundary condition within the fluxes calculated
        // Boundary condition here include: Farfield, Wall , Inlet and Outlet
        calculate_diffusive_flux(eq_type, Front, Back, Top, Bottom);

        // Finding wall flux for Immersed Cells (ib cells)
        // All other cells apart from immersed cells will have zero wall ... 
        // flux since they do not contain the immersed boundary
        calculate_wall_flux(eq_type, M_type);

        // Calculations/Iteration for fluid and solid cells
        // Its not necessary to perform calculations for solid cells
        // Can put all fluid and solid cells in their own structures ...
        // which would help removing the if condition in following function
        cell_calculations(eq_type);

        // Calculating residual using Root Mean Square
        residual = calculate_residual(eq_type, fp_residual, count_iter);    
         
        if(count_iter == 200*count_residual){
            printf("\nIteration: %d\t\tresidual = %.7lf",count_iter,  residual);
            count_residual++;
        }
        if(count_iter == 2000*count_print){
            printf("\nPrinting Files.....\nPrint count = %d", count_print);
            // Printing information extracted to files to be read by tecplot
            print_cell_data(eq_type, M_type);
            print_node_data(eq_type);
            print_ib_data();
            count_print++;
        }
        
        // Update old values for Streamfunction.
        // Used for calculating residual after each iteration
        update_old_values(eq_type);
        
        // Update iteration counter
        count_iter++;
    } while(residual > residual_crit);
    
    fclose(fp_residual);

    //-------------------------- POST PROCESSING ----------------------------
    
    sort_data();
    
    int n;

    FILE *fp = NULL;
    fp = fopen("Transfinite_Coordinates.dat", "w");
    fprintf(fp,"VARIABLES = \"x\", \"y\"\n");
    fprintf(fp, "ZONE I = %d, \t J = %d, DATAPACKING = POINT\n", 3, N_ib+1);
    for(n=1; n<=N_ib; n++){
        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_cut[0], cell_ib[n].y_cut[0]);
        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_1, cell_ib[n].y_1);
//        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_2, cell_ib[n].y_2);
        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_3, cell_ib[n].y_3);
//        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_4, cell_ib[n].y_4);
    }
    
    fprintf(fp,"%lf\t%lf\n", cell_ib[1].x_cut[0], cell_ib[1].y_cut[0]);
    fprintf(fp,"%lf\t%lf\n", cell_ib[1].x_1, cell_ib[1].y_1);
    fprintf(fp,"%lf\t%lf\n", cell_ib[1].x_3, cell_ib[1].y_3);
    
    fclose(fp);

    // Printing information extracted to files to be read by tecplot
    print_cell_data(eq_type, M_type);
    print_node_data(eq_type);
    print_ib_data();

    // Release the memory blocks allocated for Computation
    release_memory();

    // Time taken for program execution
    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine
     
    printf("\n\n --------COMPUTATION ENDS-------- \n\n\a");
    printf("\nTotal run time was %f seconds.\n\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
    
    return 0;
}
