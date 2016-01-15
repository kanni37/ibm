#include "functions.h"
#include "ibm.h"

void print_cell_data(int eq_type, int M_type)
{
    // Extract dependent variables from the solution of stream function.
    // Dependent variables: u, v, p, Cp
    switch (eq_type){
        case 1:
            switch (M_type){
                case 1:
                    extract_cell_data();
                    extract_ib_data_transfinite();
                    break;
                case 2:
                    extract_cell_data();
                    extract_ib_data_johansen_colella();
                    break;
                case 3:
                    extract_cell_data();
                    extract_ib_data_johansen_colella();
                    break;
                default:
                    printf("Error, bad input, quitting");
                    break;
            }
            print_cell_data_potential_flow();
            break;
        case 2:
            print_cell_data_navier_stokes();
            break;
        default:
            printf("Error, bad eq_type input, quitting");
            break;
    }
}

void print_node_data(int eq_type)
{
    switch (eq_type){
        case 1:
            print_node_data_potential_flow();
            break;
        case 2:
            print_node_data_navier_stokes();
            break;
        default:
            printf("Error, bad eq_type input, quitting");
            break;
    }
}

// Prints data on all cells used in flow simulation
void print_cell_data_potential_flow()
{
    FILE *fp = NULL;
    int i;
    fp = fopen("Cell_Data.dat", "w");
    fprintf(fp,"VARIABLES = \"x\", \"y\", \"Stream Function\", \"u\", \"v\", \"p\", \"state\"\n");
    fprintf(fp, "ZONE I = %d, \t J = %d, DATAPACKING = POINT\n", cell_Nx, cell_Ny);
    for(i=1; i<=Total_cells; i++){
        fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n", 
                cell[i].x, cell[i].y, cell[i].Psi[1], 
                cell[i].u[0], cell[i].v[0], cell[i].p[0], cell[i].state);
    }
    fclose(fp);
}

// Prints data on all nodes used in flow simulation
void print_node_data_potential_flow()
{
    FILE *fp = NULL;
    int i, j, k, n;

    fp = fopen("Node_Data.dat", "w");
    fprintf(fp,"VARIABLES = \"x\", \"y\", \"Stream Function\", \"u\", \"v\", \"p\", \"state\"\n");
    fprintf(fp, "ZONE I = %d, \t J = %d, DATAPACKING = POINT\n", Nx, Ny);

    for(j=1; j<=Ny; j++){
        for(i=1; i<=Nx; i++){
            k = (j-1)*Nx + i;
            n = (j-1)*cell_Nx + i;
            if(j==1){
                if(i==1){
                    point[k].Psi = cell[n].Psi[1];
                    point[k].u = cell[n].u[0];
                    point[k].v = cell[n].v[0];
                    point[k].p = cell[n].p[0];
                }else if(i==Nx){
                    point[k].Psi = cell[n-1].Psi[1];
                    point[k].u = cell[n-1].u[0];
                    point[k].v = cell[n-1].v[0];
                    point[k].p = cell[n-1].p[0];
                }else{
                    point[k].Psi = (cell[n-1].Psi[1] + cell[n].Psi[1])/2;
                    point[k].u = (cell[n-1].u[0] + cell[n].u[0])/2;
                    point[k].v = (cell[n-1].v[0] + cell[n].v[0])/2;
                    point[k].p = (cell[n-1].p[0] + cell[n].p[0])/2;
                }
            }else if(j==Ny){
                if(i==1){
                    point[k].Psi = cell[n-cell_Nx].Psi[1];
                    point[k].u = cell[n-cell_Nx].u[0];
                    point[k].v = cell[n-cell_Nx].v[0];
                    point[k].p = cell[n-cell_Nx].p[0];
                }else if(i==Nx){
                    point[k].Psi = cell[n-cell_Nx-1].Psi[1];
                    point[k].u = cell[n-cell_Nx-1].u[0];
                    point[k].v = cell[n-cell_Nx-1].v[0];
                    point[k].p = cell[n-cell_Nx-1].p[0];
                }else{
                    point[k].Psi = (cell[n-cell_Nx-1].Psi[1] + cell[n-cell_Nx].Psi[1])/2;
                    point[k].u = (cell[n-cell_Nx-1].u[0] + cell[n-cell_Nx].u[0])/2;
                    point[k].v = (cell[n-cell_Nx-1].v[0] + cell[n-cell_Nx].v[0])/2;
                    point[k].p = (cell[n-cell_Nx-1].p[0] + cell[n-cell_Nx].p[0])/2;
                }
            }else if(i==1){
                if(j==1){
                    point[k].Psi = cell[n].Psi[1];
                    point[k].u = cell[n].u[0];
                    point[k].v = cell[n].v[0];
                    point[k].p = cell[n].p[0];
                }else if(j==Ny){
                    point[k].Psi = cell[n-cell_Nx].Psi[1];
                    point[k].u = cell[n-cell_Nx].u[0];
                    point[k].v = cell[n-cell_Nx].v[0];
                    point[k].p = cell[n-cell_Nx].p[0];
                }else{
                    point[k].Psi = (cell[n].Psi[1] + cell[n-cell_Nx].Psi[1])/2;
                    point[k].u = (cell[n-cell_Nx].u[0] + cell[n].u[0])/2;
                    point[k].v = (cell[n-cell_Nx].v[0] + cell[n].v[0])/2;
                    point[k].p = (cell[n-cell_Nx].p[0] + cell[n].p[0])/2;
                }
            }else if(i==Nx){
                if(j==1){
                    point[k].Psi = cell[n-1].Psi[1];
                    point[k].u = cell[n-1].u[0];
                    point[k].v = cell[n-1].v[0];
                    point[k].p = cell[n-1].p[0];
                }else if(j==Ny){
                    point[k].Psi = cell[n-cell_Nx-1].Psi[1];
                    point[k].u = cell[n-cell_Nx-1].u[0];
                    point[k].v = cell[n-cell_Nx-1].v[0];
                    point[k].p = cell[n-cell_Nx-1].p[0];
                }else{
                    point[k].Psi = (cell[n-1].Psi[1] + cell[n-1-cell_Nx].Psi[1])/2;
                    point[k].u = (cell[n-1].u[0] + cell[n-1-cell_Nx].u[0])/2;
                    point[k].v = (cell[n-1].v[0] + cell[n-1-cell_Nx].v[0])/2;
                    point[k].p = (cell[n-1].p[0] + cell[n-1-cell_Nx].p[0])/2;
                }
            }else{
                point[k].Psi = (cell[n-1].Psi[1] + cell[n].Psi[1] + cell[n-cell_Nx-1].Psi[1] + cell[n-cell_Nx].Psi[1])/4;
                point[k].u = (cell[n-1].u[0] + cell[n].u[0] + cell[n-cell_Nx-1].u[0] + cell[n-cell_Nx].u[0])/4;
                point[k].v = (cell[n-1].v[0] + cell[n].v[0] + cell[n-cell_Nx-1].v[0] + cell[n-cell_Nx].v[0])/4;
                point[k].p = (cell[n-1].p[0] + cell[n].p[0] + cell[n-cell_Nx-1].p[0] + cell[n-cell_Nx].p[0])/4;
            }
            fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n", 
                    point[k].x, point[k].y, point[k].Psi, 
                    point[k].u, point[k].v, point[k].p, point[i].state);
        }
    }
    fclose(fp);
}


// Prints data on all cells used in flow simulation
void print_cell_data_navier_stokes()
{
    FILE *fp = NULL;
    int i;
    fp = fopen("Cell_Data.dat", "w");
    fprintf(fp,"VARIABLES = \"x\", \"y\", \"u\", \"v\", \"p\", \"state\"\n");
    fprintf(fp, "ZONE I = %d, \t J = %d, DATAPACKING = POINT\n", cell_Nx, cell_Ny);
    for(i=1; i<=Total_cells; i++){
        fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%d\n", 
                cell[i].x, cell[i].y, 
                cell[i].U[1][0], cell[i].U[2][0], cell[i].U[0][0], cell[i].state);
    }
    fclose(fp);
}

// Prints data on all nodes used in flow simulation
void print_node_data_navier_stokes()
{
    FILE *fp = NULL;
    int i, j, k, n;

    fp = fopen("Node_Data.dat", "w");
    fprintf(fp,"VARIABLES = \"x\", \"y\", \"u\", \"v\", \"p\", \"state\"\n");
    fprintf(fp, "ZONE I = %d, \t J = %d, DATAPACKING = POINT\n", Nx, Ny);

    for(j=1; j<=Ny; j++){
        for(i=1; i<=Nx; i++){
            k = (j-1)*Nx + i;
            n = (j-1)*cell_Nx + i;
            if(j==1){
                if(i==1){
                    point[k].u = cell[n].U[1][0];
                    point[k].v = cell[n].U[2][0];
                    point[k].p = cell[n].U[0][0];
                }else if(i==Nx){
                    point[k].u = cell[n-1].U[1][0];
                    point[k].v = cell[n-1].U[2][0];
                    point[k].p = cell[n-1].U[0][0];
                }else{
                    point[k].u = (cell[n-1].U[1][0] + cell[n].U[1][0])/2;
                    point[k].v = (cell[n-1].U[2][0] + cell[n].U[2][0])/2;
                    point[k].p = (cell[n-1].U[0][0] + cell[n].U[0][0])/2;
                }
            }else if(j==Ny){
                if(i==1){
                    point[k].u = cell[n-cell_Nx].U[1][0];
                    point[k].v = cell[n-cell_Nx].U[2][0];
                    point[k].p = cell[n-cell_Nx].U[0][0];
                }else if(i==Nx){
                    point[k].u = cell[n-cell_Nx-1].U[1][0];
                    point[k].v = cell[n-cell_Nx-1].U[2][0];
                    point[k].p = cell[n-cell_Nx-1].U[0][0];
                }else{
                    point[k].u = (cell[n-cell_Nx-1].U[1][0] + cell[n-cell_Nx].U[1][0])/2;
                    point[k].v = (cell[n-cell_Nx-1].U[2][0] + cell[n-cell_Nx].U[2][0])/2;
                    point[k].p = (cell[n-cell_Nx-1].U[0][0] + cell[n-cell_Nx].U[0][0])/2;
                }
            }else if(i==1){
                if(j==1){
                    point[k].u = cell[n].U[1][0];
                    point[k].v = cell[n].U[2][0];
                    point[k].p = cell[n].U[0][0];
                }else if(j==Ny){
                    point[k].u = cell[n-cell_Nx].U[1][0];
                    point[k].v = cell[n-cell_Nx].U[2][0];
                    point[k].p = cell[n-cell_Nx].U[0][0];
                }else{
                    point[k].u = (cell[n-cell_Nx].U[1][0] + cell[n].U[1][0])/2;
                    point[k].v = (cell[n-cell_Nx].U[2][0] + cell[n].U[2][0])/2;
                    point[k].p = (cell[n-cell_Nx].U[0][0] + cell[n].U[0][0])/2;
                }
            }else if(i==Nx){
                if(j==1){
                    point[k].u = cell[n-1].U[1][0];
                    point[k].v = cell[n-1].U[2][0];
                    point[k].p = cell[n-1].U[0][0];
                }else if(j==Ny){
                    point[k].u = cell[n-cell_Nx-1].U[1][0];
                    point[k].v = cell[n-cell_Nx-1].U[2][0];
                    point[k].p = cell[n-cell_Nx-1].U[0][0];
                }else{
                    point[k].u = (cell[n-1].U[1][0] + cell[n-1-cell_Nx].U[1][0])/2;
                    point[k].v = (cell[n-1].U[2][0] + cell[n-1-cell_Nx].U[2][0])/2;
                    point[k].p = (cell[n-1].U[0][0] + cell[n-1-cell_Nx].U[0][0])/2;
                }
            }else{
                point[k].u = (cell[n-1].U[1][0] + cell[n].U[1][0] + cell[n-cell_Nx-1].U[1][0] + cell[n-cell_Nx].U[1][0])/4;
                point[k].v = (cell[n-1].U[2][0] + cell[n].U[2][0] + cell[n-cell_Nx-1].U[2][0] + cell[n-cell_Nx].U[2][0])/4;
                point[k].p = (cell[n-1].U[0][0] + cell[n].U[0][0] + cell[n-cell_Nx-1].U[0][0] + cell[n-cell_Nx].U[0][0])/4;
            }
            fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%d\n", 
                    point[k].x, point[k].y, 
                    point[k].u, point[k].v, point[k].p, point[k].state);
        }
    }
    fclose(fp);
}

// Prints data on immersed boundary/surface
void print_ib_data()
{
    int n;
    FILE *fp_1, *fp_2, *fp_3, *fp_4, *fp_5;
    double Cp;
    double u, v, p;
    double theta;

    char filename1[sizeof "IB_Data_99.dat"];
    char filename2[sizeof "IB_Data_analytical.dat"];

    sprintf(filename1, "IB_Data_%02d.dat", N_division);
    sprintf(filename2, "IB_Data_analytical.dat");

    fp_1 = fopen(filename1, "w");
    fp_2 = fopen(filename2, "w");
    fp_3 = fopen("Absolute_Error.dat", "w");
    fp_4 = fopen("Relative_Error.dat", "w");
    fp_5 = fopen("Relative_Error_max.dat", "w");

    fprintf(fp_1, "VARIABLES = \"x\", \"y\", \"theta\", \"u\", \"v\", \"p\", \"Cp\"\n");
    fprintf(fp_2, "VARIABLES = \"x\", \"y\", \"theta\", \"u\", \"v\", \"p\", \"Cp\"\n");
    fprintf(fp_3, "VARIABLES = \"x\", \"y\", \"theta\", \"Error_u\", \"Error_v\", \"Error_p\", \"Error_Cp\"\n");
    fprintf(fp_4, "VARIABLES = \"x\", \"y\", \"theta\", \"Error_u\", \"Error_v\", \"Error_p\", \"Error_Cp\"\n");
    fprintf(fp_5, "VARIABLES = \"x\", \"y\", \"theta\", \"Error_u\", \"Error_v\", \"Error_p\", \"Error_Cp\"\n");

    for(n=1; n<=N_ib; n++){
        if(cell_ib[n].d_cut > d_cut_tol){
            theta = cell_ib[n].theta;
            Cp = 1 - 4*sin(theta)*sin(theta);
            u = 2*U_ini*sin(theta)*sin(theta); 
            v = -2*U_ini*sin(theta)*cos(theta);
            p = 0.5*(U_ini*U_ini - u*u -v*v);
            fprintf(fp_1, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
                    cell_ib[n].x_m, cell_ib[n].y_m, theta*180/pi, cell_ib[n].u[0], cell_ib[n].v[0], cell_ib[n].p[0], cell_ib[n].Cp[0]);
            fprintf(fp_2, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
                    cell_ib[n].x_m, cell_ib[n].y_m, theta*180/pi, u, v ,p, Cp);
            fprintf(fp_3, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
                    cell_ib[n].x_m, cell_ib[n].y_m, theta*180/pi, fabs(cell_ib[n].u[0] - u), fabs(cell_ib[n].v[0] - v), fabs(cell_ib[n].p[0] - p), fabs(cell_ib[n].Cp[0] - Cp));
            fprintf(fp_4, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
                    cell_ib[n].x_m, cell_ib[n].y_m, theta*180/pi, fabs((cell_ib[n].u[0] - u)/u*100), fabs((cell_ib[n].v[0] - v)/v*100), fabs((cell_ib[n].p[0] - p)/p*100), fabs((cell_ib[n].Cp[0] - Cp)/Cp*100));
            fprintf(fp_5, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
                    cell_ib[n].x_m, cell_ib[n].y_m, theta*180/pi, fabs((cell_ib[n].u[0] - u)/(2*U_ini)*100), 
                    fabs((cell_ib[n].v[0] - v)/(U_ini)*100), fabs((cell_ib[n].p[0] - p)/(1.5*U_ini*U_ini)*100), fabs((cell_ib[n].Cp[0] - Cp)/3*100));
        }
    }

    fclose(fp_1);
    fclose(fp_2);
    fclose(fp_3);
    fclose(fp_4);
    fclose(fp_5);
}
