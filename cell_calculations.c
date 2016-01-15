#include "functions.h"
#include "ibm.h"

void update_old_values(int eq_type)
{
    switch (eq_type){
        case 1:
            update_old_values_potential_flow();
            break;
        case 2:
            update_old_values_navier_stokes();
            break;
        default:
            printf("Error, bad eq_type input, quitting");
            break;
    }
}

void cell_calculations(int eq_type)
{
    switch (eq_type){
        case 1:
            cell_calculations_potential_flow();
            break;
        case 2:
            cell_calculations_navier_stokes();
            break;
        default:
            printf("Error, bad eq_type input, quitting");
            break;
    }
}

void update_old_values_potential_flow()
{
    int n;

    //Updating Old values
#pragma omp parallel for private(n) schedule(auto)
    for(n=1; n<=Total_cells; n++){
        cell[n].Psi[0] = cell[n].Psi[1];
    }
    for(n=1; n<=N_ib; n++){
        cell_ib[n].Psi[0] = cell_ib[n].Psi[1];
    }
}

void cell_calculations_potential_flow()
{
    int n, n_ib;

#pragma omp parallel for private(n_ib) schedule(auto)
    for(n=1; n<=Total_cells; n++){
        if(cell[n].state != 2){
            cell[n].Psi[1] = cell[n].Psi[0] + dt*acf/(dx*dy)*((cell[n].Fv_e[0] - cell[n].Fv_w[0])*dy + (cell[n].Gv_n[0] - cell[n].Gv_s[0])*dx);
            cell[n].Psi[1] = w*cell[n].Psi[1] + (1-w)*cell[n].Psi[0];
        }else if(cell[n].state == 2){
            n_ib = cell[n].ib_cell_no;
            cell[n].Psi[1] = cell[n].Psi[0] + dt*acf/(dx*dy)*((cell[n].Fv_e[0]*cell_ib[n_ib].alpha[1] - cell[n].Fv_w[0]*cell_ib[n_ib].alpha[3])*dy + 
                    (cell[n].Gv_n[0]*cell_ib[n_ib].alpha[2] - cell[n].Gv_s[0]*cell_ib[n_ib].alpha[0])*dx + cell_ib[n_ib].wall_flux_diffusive[0]);
            cell[n].Psi[1] = w*cell[n].Psi[1] + (1-w)*cell[n].Psi[0];
            // updating value in cell_ib cell structure simultaneously
            cell_ib[n_ib].Psi[1] = cell[n].Psi[1];
        }
    }
}

void update_old_values_navier_stokes()
{
    int n;

#pragma omp parallel for private(n) schedule(auto)
    //Updating Old values
    for(n=1; n<=Total_cells; n++){
        cell[n].U[0][0] = cell[n].U[0][1];
        cell[n].U[1][0] = cell[n].U[1][1];
        cell[n].U[2][0] = cell[n].U[2][1];
    }
    for(n=1; n<=N_ib; n++){
        cell_ib[n].U[0][0] = cell_ib[n].U[0][1];
        cell_ib[n].U[1][0] = cell_ib[n].U[1][1];
        cell_ib[n].U[2][0] = cell_ib[n].U[2][1];
    }
}

void cell_calculations_navier_stokes()
{
    int n, n_ib;
    //int t_division = 50;
    //int t_count;

#pragma omp parallel for private(n_ib) schedule(auto)
    // Cell calculations
    for(n=1; n<=Total_cells; n++){
        if(cell[n].state != 2){
            cell[n].U[0][1] = cell[n].U[0][0] + acf*dt/(dx*dy)*(-(cell[n].F_e[0] - cell[n].F_w[0])*dy - (cell[n].G_n[0] - cell[n].G_s[0])*dx +
                    (cell[n].Fv_e[0] - cell[n].Fv_w[0])*dy + (cell[n].Gv_n[0] - cell[n].Gv_s[0])*dx );
            cell[n].U[1][1] = cell[n].U[1][0] + dt/(dx*dy)*(-(cell[n].F_e[1] - cell[n].F_w[1])*dy - (cell[n].G_n[1] - cell[n].G_s[1])*dx + 
                    (cell[n].Fv_e[1] - cell[n].Fv_w[1])*dy + (cell[n].Gv_n[1] - cell[n].Gv_s[1])*dx);
            cell[n].U[2][1] = cell[n].U[2][0] + dt/(dx*dy)*(-(cell[n].F_e[2] - cell[n].F_w[2])*dy - (cell[n].G_n[2] - cell[n].G_s[2])*dx + 
                    (cell[n].Fv_e[2] - cell[n].Fv_w[2])*dy + (cell[n].Gv_n[2] - cell[n].Gv_s[2])*dx);
        }else if(cell[n].state == 2){
           // for(t_count = 1; t_count<=t_division; t_count++){
                //Insert IB cell Calculations here
                n_ib = cell[n].ib_cell_no;

                cell[n].U[0][1] = cell[n].U[0][0] + 
                    acf*dt/(dx*dy)*(-(cell[n].F_e[0]*cell_ib[n_ib].alpha[1] - cell[n].F_w[0]*cell_ib[n_ib].alpha[3])*dy - 
                            (cell[n].G_n[0]*cell_ib[n_ib].alpha[2] - cell[n].G_s[0]*cell_ib[n_ib].alpha[0])*dx -
                            cell_ib[n_ib].wall_flux_convective[0] +
                            (cell[n].Fv_e[0]*cell_ib[n_ib].alpha[1] - cell[n].Fv_w[0]*cell_ib[n_ib].alpha[3])*dy + 
                            (cell[n].Gv_n[0]*cell_ib[n_ib].alpha[2] - cell[n].Gv_s[0]*cell_ib[n_ib].alpha[0])*dx +
                            cell_ib[n_ib].wall_flux_diffusive[0]);
                cell[n].U[1][1] = cell[n].U[1][0] + 
                    dt/(dx*dy)*(-(cell[n].F_e[1]*cell_ib[n_ib].alpha[1] - cell[n].F_w[1]*cell_ib[n_ib].alpha[3])*dy - 
                            (cell[n].G_n[1]*cell_ib[n_ib].alpha[2] - cell[n].G_s[1]*cell_ib[n_ib].alpha[0])*dx -
                            cell_ib[n_ib].wall_flux_convective[1] +
                            (cell[n].Fv_e[1]*cell_ib[n_ib].alpha[1] - cell[n].Fv_w[1]*cell_ib[n_ib].alpha[3])*dy + 
                            (cell[n].Gv_n[1]*cell_ib[n_ib].alpha[2] - cell[n].Gv_s[1]*cell_ib[n_ib].alpha[0])*dx +
                            cell_ib[n_ib].wall_flux_diffusive[1]);
                cell[n].U[2][1] = cell[n].U[2][0] + 
                    dt/(dx*dy)*(-(cell[n].F_e[2]*cell_ib[n_ib].alpha[1] - cell[n].F_w[2]*cell_ib[n_ib].alpha[3])*dy - 
                            (cell[n].G_n[2]*cell_ib[n_ib].alpha[2] - cell[n].G_s[2]*cell_ib[n_ib].alpha[0])*dx -
                            cell_ib[n_ib].wall_flux_convective[2] +
                            (cell[n].Fv_e[2]*cell_ib[n_ib].alpha[1] - cell[n].Fv_w[2]*cell_ib[n_ib].alpha[3])*dy + 
                            (cell[n].Gv_n[2]*cell_ib[n_ib].alpha[2] - cell[n].Gv_s[2]*cell_ib[n_ib].alpha[0])*dx +
                            cell_ib[n_ib].wall_flux_diffusive[2]);

                //Updating value in cell_ib cell structure simulatneously
                cell_ib[n_ib].U[0][1] = cell[n].U[0][1]; 
                cell_ib[n_ib].U[1][1] = cell[n].U[1][1]; 
                cell_ib[n_ib].U[2][1] = cell[n].U[2][1]; 
           // }
        }
    }



#pragma omp parallel for schedule(auto)
    // Successive Over-Relaxation
    for(n=1; n<=Total_cells; n++){
        cell[n].U[0][1] = w*cell[n].U[0][1] + (1-w)*cell[n].U[0][0];
        cell[n].U[1][1] = w*cell[n].U[1][1] + (1-w)*cell[n].U[1][0];
        cell[n].U[2][1] = w*cell[n].U[2][1] + (1-w)*cell[n].U[2][0];
    }

}
