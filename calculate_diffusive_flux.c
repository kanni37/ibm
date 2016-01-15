#include "functions.h"
#include "ibm.h"

void calculate_diffusive_flux(int eq_type, int Front, int Back, int Top, int Bottom)
{
    switch (eq_type){
        case 1:
            calculate_diffusive_flux_potential_flow();
            break;
        case 2:
            calculate_diffusive_flux_navier_stokes(Front, Back, Top, Bottom);
            break;
        default:
            printf("Error, bad eq_type input, quitting");
            break;
    }
}


void calculate_diffusive_flux_potential_flow()
{
    int i, j, n, k, n_ib;

    #pragma omp parallel for private(i, j, k, n, n_ib) schedule(auto)
    //Interior Fluid and Solid Cells
    for(j=2; j<cell_Ny;j++){
        for(i=2; i<cell_Nx; i++){
            n = (j-1)*cell_Nx + i;
            if(cell[n].state != 2){
                // East Face
                cell[n].Fv_e[0] = (cell[n+1].Psi[0] - cell[n].Psi[0])/dx;
                // West Face
                cell[n].Fv_w[0] = (cell[n].Psi[0] - cell[n-1].Psi[0])/dx;
                // North Face 
                cell[n].Gv_n[0] = (cell[n+cell_Nx].Psi[0] - cell[n].Psi[0])/dy;
                // South Face
                cell[n].Gv_s[0] = (cell[n].Psi[0] - cell[n-cell_Nx].Psi[0])/dy;
            }else if(cell[n].state == 2){
                n_ib = cell[n].ib_cell_no;
                // East Face
                k = 1;
                cell[n].Fv_e[0] = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].Psi[0] -  cell[cell_ib[n_ib].BI_flux_cell[k][3]].Psi[0])/dx*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].Psi[0] -  cell[cell_ib[n_ib].BI_flux_cell[k][0]].Psi[0])/dx*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                // West Face
                k = 3;
                cell[n].Fv_w[0] = -((cell[cell_ib[n_ib].BI_flux_cell[k][2]].Psi[0] -  cell[cell_ib[n_ib].BI_flux_cell[k][3]].Psi[0])/dx*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].Psi[0] -  cell[cell_ib[n_ib].BI_flux_cell[k][0]].Psi[0])/dx*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                // North Face 
                k = 2;
                cell[n].Gv_n[0] = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].Psi[0] -  cell[cell_ib[n_ib].BI_flux_cell[k][1]].Psi[0])/dy*(cell_ib[n_ib].dx2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][3]].Psi[0] -  cell[cell_ib[n_ib].BI_flux_cell[k][0]].Psi[0])/dy*(cell_ib[n_ib].dx1[k]))/
                    (cell_ib[n_ib].dx1[k] + cell_ib[n_ib].dx2[k]); 
                // South Face
                k = 0;
                cell[n].Gv_s[0] = -((cell[cell_ib[n_ib].BI_flux_cell[k][2]].Psi[0] -  cell[cell_ib[n_ib].BI_flux_cell[k][1]].Psi[0])/dy*(cell_ib[n_ib].dx2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][3]].Psi[0] -  cell[cell_ib[n_ib].BI_flux_cell[k][0]].Psi[0])/dy*(cell_ib[n_ib].dx1[k]))/
                    (cell_ib[n_ib].dx1[k] + cell_ib[n_ib].dx2[k]); 
            }
        }
    }
    //Inlet
    for(j=2; j<=cell_Ny-1; j++){
        n = (j-1)*cell_Nx+1;
        // East Face
        cell[n].Fv_e[0] = (cell[n+1].Psi[0] - cell[n].Psi[0])/dx;
        // West Face
        cell[n].Fv_w[0] = 0; 
        // North Face 
        cell[n].Gv_n[0] = (cell[n+cell_Nx].Psi[0] - cell[n].Psi[0])/dy;
        // South Face
        cell[n].Gv_s[0] = (cell[n].Psi[0] - cell[n-cell_Nx].Psi[0])/dy;
    }
    //Outlet
    for(j=2; j<=cell_Ny-1; j++){
        n = j*cell_Nx;
        // East Face
        cell[n].Fv_e[0] = 0; 
        // West Face
        cell[n].Fv_w[0] = (cell[n].Psi[0] - cell[n-1].Psi[0])/dx;
        // North Face 
        cell[n].Gv_n[0] = (cell[n+cell_Nx].Psi[0] - cell[n].Psi[0])/dy;
        // South Face
        cell[n].Gv_s[0] = (cell[n].Psi[0] - cell[n-cell_Nx].Psi[0])/dy;
    }
    //Farfield Bottom
    for(i=1; i<=cell_Nx; i++){
        n = i;
        if (i == 1) {
            // East Face
            cell[n].Fv_e[0] = (cell[n+1].Psi[0] - cell[n].Psi[0])/dx;
            // West Face
            cell[n].Fv_w[0] = 0; 
            // North Face 
            cell[n].Gv_n[0] = (cell[n+cell_Nx].Psi[0] - cell[n].Psi[0])/dy;
            // South Face
            cell[n].Gv_s[0] = U_ini; 
        }else if (i == cell_Nx){
            // East Face
            cell[n].Fv_e[0] = 0;
            // West Face
            cell[n].Fv_w[0] = (cell[n].Psi[0] - cell[n-1].Psi[0])/dx;
            // North Face 
            cell[n].Gv_n[0] = (cell[n+cell_Nx].Psi[0] - cell[n].Psi[0])/dy;
            // South Face
            cell[n].Gv_s[0] = U_ini; 
        }else{
            // East Face
            cell[n].Fv_e[0] = (cell[n+1].Psi[0] - cell[n].Psi[0])/dx;
            // West Face
            cell[n].Fv_w[0] = (cell[n].Psi[0] - cell[n-1].Psi[0])/dx;
            // North Face 
            cell[n].Gv_n[0] = (cell[n+cell_Nx].Psi[0] - cell[n].Psi[0])/dy;
            // South Face
            cell[n].Gv_s[0] = U_ini; 
        }
    }
    //Farfield Top
    for(i=1; i<=cell_Nx; i++){
        n=(cell_Ny-1)*cell_Nx + i;
        if (i == 1) {
            // East Face
            cell[n].Fv_e[0] = (cell[n+1].Psi[0] - cell[n].Psi[0])/dx;
            // West Face
            cell[n].Fv_w[0] = 0; 
            // North Face 
            cell[n].Gv_n[0] = U_ini; 
            // South Face
            cell[n].Gv_s[0] = (cell[n].Psi[0] - cell[n-cell_Nx].Psi[0])/dy;
        }else if (i == cell_Nx){
            // East Face
            cell[n].Fv_e[0] = 0; 
            // West Face
            cell[n].Fv_w[0] = (cell[n].Psi[0] - cell[n-1].Psi[0])/dx;
            // North Face 
            cell[n].Gv_n[0] = U_ini; 
            // South Face
            cell[n].Gv_s[0] = (cell[n].Psi[0] - cell[n-cell_Nx].Psi[0])/dy;
        }else{
            // East Face
            cell[n].Fv_e[0] = (cell[n+1].Psi[0] - cell[n].Psi[0])/dx;
            // West Face
            cell[n].Fv_w[0] = (cell[n].Psi[0] - cell[n-1].Psi[0])/dx;
            // North Face 
            cell[n].Gv_n[0] = U_ini; 
            // South Face
            cell[n].Gv_s[0] = (cell[n].Psi[0] - cell[n-cell_Nx].Psi[0])/dy;
        }
    }
}

void calculate_diffusive_flux_navier_stokes(int Front, int Back, int Top, int Bottom)
{
    int i, j, n, k, n_ib;

    #pragma omp parallel for private(i, j, n, k, n_ib) schedule(auto)
    //Interior Fluid and Solid Cells
    for(j=2; j<cell_Ny;j++){
        for(i=2; i<cell_Nx; i++){
            n = (j-1)*cell_Nx + i;
            if(cell[n].state != 2){
                // East Face
                cell[n].Fv_e[0] = 0; 
                cell[n].Fv_e[1] = myu/rho*(cell[n+1].U[1][0] - cell[n].U[1][0])/dx;
                cell[n].Fv_e[2] = myu/rho*(cell[n+1].U[2][0] - cell[n].U[2][0])/dx;
                // West Face
                cell[n].Fv_w[0] = 0; 
                cell[n].Fv_w[1] = myu/rho*(cell[n].U[1][0] - cell[n-1].U[1][0])/dx;
                cell[n].Fv_w[2] = myu/rho*(cell[n].U[2][0] - cell[n-1].U[2][0])/dx;
                // North Face
                cell[n].Gv_n[0] = 0;
                cell[n].Gv_n[1] = myu/rho*(cell[n+cell_Nx].U[1][0] - cell[n].U[1][0])/dy;
                cell[n].Gv_n[2] = myu/rho*(cell[n+cell_Nx].U[2][0] - cell[n].U[2][0])/dy;
                // South Face
                cell[n].Gv_s[0] = 0;
                cell[n].Gv_s[1] = myu/rho*(cell[n].U[1][0] - cell[n-cell_Nx].U[1][0])/dy;
                cell[n].Gv_s[2] = myu/rho*(cell[n].U[2][0] - cell[n-cell_Nx].U[2][0])/dy;
            }else if(cell[n].state == 2){
                n_ib = cell[n].ib_cell_no;
                // East Face
                k = 1;
                cell[n].Fv_e[0] = 0;
                cell[n].Fv_e[1] = myu/rho*((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[1][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[1][0])/dx*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[1][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[1][0])/dx*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                cell[n].Fv_e[2] = myu/rho*((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[2][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[2][0])/dx*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[2][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[2][0])/dx*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                // West Face
                k = 3;
                cell[n].Fv_w[0] = 0; 
                cell[n].Fv_w[1] = -myu/rho*((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[1][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[1][0])/dx*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[1][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[1][0])/dx*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                cell[n].Fv_w[2] = -myu/rho*((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[2][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[2][0])/dx*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[2][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[2][0])/dx*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                // North Face 
                k = 2;
                cell[n].Gv_n[0] = 0; 
                cell[n].Gv_n[1] = myu/rho*((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[1][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[1][0])/dy*(cell_ib[n_ib].dx2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[1][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[1][0])/dy*(cell_ib[n_ib].dx1[k]))/
                    (cell_ib[n_ib].dx1[k] + cell_ib[n_ib].dx2[k]); 
                cell[n].Gv_n[2] = myu/rho*((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[2][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[2][0])/dy*(cell_ib[n_ib].dx2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[2][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[2][0])/dy*(cell_ib[n_ib].dx1[k]))/
                    (cell_ib[n_ib].dx1[k] + cell_ib[n_ib].dx2[k]); 
                // South Face
                k = 0;
                cell[n].Gv_s[0] = 0; 
                cell[n].Gv_s[1] = -myu/rho*((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[1][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[1][0])/dy*(cell_ib[n_ib].dx2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[1][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[1][0])/dy*(cell_ib[n_ib].dx1[k]))/
                    (cell_ib[n_ib].dx1[k] + cell_ib[n_ib].dx2[k]); 
                cell[n].Gv_s[2] = -myu/rho*((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[2][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[2][0])/dy*(cell_ib[n_ib].dx2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[2][0] -  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[2][0])/dy*(cell_ib[n_ib].dx1[k]))/
                    (cell_ib[n_ib].dx1[k] + cell_ib[n_ib].dx2[k]); 
            }
        }
    }
    //Inlet
    for(j=2; j<cell_Ny; j++){
        n = 1 + (j-1)*cell_Nx;
        // East Face
        cell[n].Fv_e[0] = 0; 
        cell[n].Fv_e[1] = myu/rho*(cell[n+1].U[1][0] - cell[n].U[1][0])/dx;
        cell[n].Fv_e[2] = myu/rho*(cell[n+1].U[2][0] - cell[n].U[2][0])/dx;
        // West Face
        cell[n].Fv_w[0] = 0; 
        cell[n].Fv_w[1] = myu/rho*(cell[n].U[1][0] - U_ini)/(dx/2);
        cell[n].Fv_w[2] = myu/rho*(cell[n].U[2][0] - V_ini)/(dx/2);
        // North Face
        cell[n].Gv_n[0] = 0;
        cell[n].Gv_n[1] = myu/rho*(cell[n+cell_Nx].U[1][0] - cell[n].U[1][0])/dy;
        cell[n].Gv_n[2] = myu/rho*(cell[n+cell_Nx].U[2][0] - cell[n].U[2][0])/dy;
        // South Face
        cell[n].Gv_s[0] = 0;
        cell[n].Gv_s[1] = myu/rho*(cell[n].U[1][0] - cell[n-cell_Nx].U[1][0])/dy;
        cell[n].Gv_s[2] = myu/rho*(cell[n].U[2][0] - cell[n-cell_Nx].U[2][0])/dy;
    } 
    //Outlet
    for(j=2; j<cell_Ny; j++){
        n = j*cell_Nx;
        // East Face
        cell[n].Fv_e[0] = 0; 
        cell[n].Fv_e[1] = 0;
        cell[n].Fv_e[2] = 0;
        // West Face
        cell[n].Fv_w[0] = 0; 
        cell[n].Fv_w[1] = myu/rho*(cell[n].U[1][0] - cell[n-1].U[1][0])/dx;
        cell[n].Fv_w[2] = myu/rho*(cell[n].U[2][0] - cell[n-1].U[2][0])/dx;
        // North Face
        cell[n].Gv_n[0] = 0;
        cell[n].Gv_n[1] = myu/rho*(cell[n+cell_Nx].U[1][0] - cell[n].U[1][0])/dy;
        cell[n].Gv_n[2] = myu/rho*(cell[n+cell_Nx].U[2][0] - cell[n].U[2][0])/dy;
        // South Face
        cell[n].Gv_s[0] = 0;
        cell[n].Gv_s[1] = myu/rho*(cell[n].U[1][0] - cell[n-cell_Nx].U[1][0])/dy;
        cell[n].Gv_s[2] = myu/rho*(cell[n].U[2][0] - cell[n-cell_Nx].U[2][0])/dy;
    } 
    //Farfield/Wall Bottom
    for(i=1; i<=cell_Nx; i++){
        n = i;
        if(i == 1){
            // East Face
            cell[n].Fv_e[0] = 0; 
            cell[n].Fv_e[1] = myu/rho*(cell[n+1].U[1][0] - cell[n].U[1][0])/dx;
            cell[n].Fv_e[2] = myu/rho*(cell[n+1].U[2][0] - cell[n].U[2][0])/dx;
            // West Face
            cell[n].Fv_w[0] = 0; 
            cell[n].Fv_w[1] = myu/rho*(cell[n].U[1][0] - U_ini)/(dx/2);
            cell[n].Fv_w[2] = myu/rho*(cell[n].U[2][0] - V_ini)/(dx/2);
            // North Face
            cell[n].Gv_n[0] = 0;
            cell[n].Gv_n[1] = myu/rho*(cell[n+cell_Nx].U[1][0] - cell[n].U[1][0])/dy;
            cell[n].Gv_n[2] = myu/rho*(cell[n+cell_Nx].U[2][0] - cell[n].U[2][0])/dy;
            // South Face
            if(Bottom == 1){            // Farfield
                cell[n].Gv_s[0] = 0;
                cell[n].Gv_s[1] = myu/rho*(cell[n].U[1][0] - U_ini)/(dy/2);
                cell[n].Gv_s[2] = myu/rho*(cell[n].U[2][0] - V_ini)/(dy/2);
            }else if(Bottom == 2){      // Wall
                cell[n].Gv_s[0] = 0;
                cell[n].Gv_s[1] = myu/rho*(cell[n].U[1][0] - 0)/(dy/2);
                cell[n].Gv_s[2] = myu/rho*(cell[n].U[2][0] - 0)/(dy/2);
            }else if(Bottom == 3){      // Symmetry
                cell[n].Gv_s[0] = 0;
                cell[n].Gv_s[1] = 0; 
                cell[n].Gv_s[2] = 0;
            }
        }else if(i == cell_Nx){
            // East Face
            cell[n].Fv_e[0] = 0; 
            cell[n].Fv_e[1] = 0; 
            cell[n].Fv_e[2] = 0;
            // West Face
            cell[n].Fv_w[0] = 0; 
            cell[n].Fv_w[1] = myu/rho*(cell[n].U[1][0] - cell[n-1].U[1][0])/dx;
            cell[n].Fv_w[2] = myu/rho*(cell[n].U[2][0] - cell[n-1].U[2][0])/dx;
            // North Face
            cell[n].Gv_n[0] = 0;
            cell[n].Gv_n[1] = myu/rho*(cell[n+cell_Nx].U[1][0] - cell[n].U[1][0])/dy;
            cell[n].Gv_n[2] = myu/rho*(cell[n+cell_Nx].U[2][0] - cell[n].U[2][0])/dy;
            // South Face
            if(Bottom == 1){            // Farfield
                cell[n].Gv_s[0] = 0;
                cell[n].Gv_s[1] = myu/rho*(cell[n].U[1][0] - U_ini)/(dy/2);
                cell[n].Gv_s[2] = myu/rho*(cell[n].U[2][0] - V_ini)/(dy/2);
            }else if(Bottom == 2){      // Wall
                cell[n].Gv_s[0] = 0;
                cell[n].Gv_s[1] = myu/rho*(cell[n].U[1][0] - 0)/(dy/2);
                cell[n].Gv_s[2] = myu/rho*(cell[n].U[2][0] - 0)/(dy/2);
            }else if(Bottom == 3){  // Symmetry
                cell[n].Gv_s[0] = 0;
                cell[n].Gv_s[1] = 0; 
                cell[n].Gv_s[2] = 0;
            }
        }else{
            // East Face
            cell[n].Fv_e[0] = 0; 
            cell[n].Fv_e[1] = myu/rho*(cell[n+1].U[1][0] - cell[n].U[1][0])/dx;
            cell[n].Fv_e[2] = myu/rho*(cell[n+1].U[2][0] - cell[n].U[2][0])/dx;
            // West Face
            cell[n].Fv_w[0] = myu/rho*0; 
            cell[n].Fv_w[1] = myu/rho*(cell[n].U[1][0] - cell[n-1].U[1][0])/dx;
            cell[n].Fv_w[2] = myu/rho*(cell[n].U[2][0] - cell[n-1].U[2][0])/dx;
            // North Face
            cell[n].Gv_n[0] = 0;
            cell[n].Gv_n[1] = myu/rho*(cell[n+cell_Nx].U[1][0] - cell[n].U[1][0])/dy;
            cell[n].Gv_n[2] = myu/rho*(cell[n+cell_Nx].U[2][0] - cell[n].U[2][0])/dy;
            // South Face
            if(Bottom == 1){            // Farfield
                cell[n].Gv_s[0] = 0;
                cell[n].Gv_s[1] = myu/rho*(cell[n].U[1][0] - U_ini)/(dy/2);
                cell[n].Gv_s[2] = myu/rho*(cell[n].U[2][0] - V_ini)/(dy/2);
            }else if(Bottom == 2){      // Wall
                cell[n].Gv_s[0] = 0;
                cell[n].Gv_s[1] = myu/rho*(cell[n].U[1][0] - 0)/(dy/2);
                cell[n].Gv_s[2] = myu/rho*(cell[n].U[2][0] - 0)/(dy/2);
            }else if(Bottom == 3){  // Symmetry
                cell[n].Gv_s[0] = 0;
                cell[n].Gv_s[1] = 0; 
                cell[n].Gv_s[2] = 0;
            }
        } 
    }
    //Farfield/Wall Top
    for(i=1; i<=cell_Nx; i++){
        n = (cell_Ny-1)*cell_Nx+i;
        if(i == 1){
            // East Face
            cell[n].Fv_e[0] = 0; 
            cell[n].Fv_e[1] = myu/rho*(cell[n+1].U[1][0] - cell[n].U[1][0])/dx;
            cell[n].Fv_e[2] = myu/rho*(cell[n+1].U[2][0] - cell[n].U[2][0])/dx;
            // West Face
            cell[n].Fv_w[0] = 0; 
            cell[n].Fv_w[1] = myu/rho*(cell[n].U[1][0] - U_ini)/(dx/2);
            cell[n].Fv_w[2] = myu/rho*(cell[n].U[2][0] - V_ini)/(dx/2);
            // North Face
            if(Top == 1){            // Farfield
                cell[n].Gv_n[0] = 0;
                cell[n].Gv_n[1] = myu/rho*(U_ini - cell[n].U[1][0])/(dy/2);
                cell[n].Gv_n[2] = myu/rho*(V_ini - cell[n].U[2][0])/(dy/2);
            }else if(Top == 2){      // Wall
                cell[n].Gv_n[0] = 0;
                cell[n].Gv_n[1] = myu/rho*(0 - cell[n].U[1][0])/(dy/2);
                cell[n].Gv_n[2] = myu/rho*(0 - cell[n].U[2][0])/(dy/2);
            }else if(Top == 3){
                cell[n].Gv_n[0] = 0;
                cell[n].Gv_n[1] = 0; 
                cell[n].Gv_n[2] = 0;
            }
            // South Face
            cell[n].Gv_s[0] = 0;
            cell[n].Gv_s[1] = myu/rho*(cell[n].U[1][0] - cell[n-cell_Nx].U[1][0])/dy;
            cell[n].Gv_s[2] = myu/rho*(cell[n].U[2][0] - cell[n-cell_Nx].U[2][0])/dy;
        }else if(i == cell_Nx){
            // East Face
            cell[n].Fv_e[0] = 0; 
            cell[n].Fv_e[1] = 0; 
            cell[n].Fv_e[2] = 0; 
            // West Face
            cell[n].Fv_w[0] = 0; 
            cell[n].Fv_w[1] = myu/rho*(cell[n].U[1][0] - cell[n-1].U[1][0])/dx;
            cell[n].Fv_w[2] = myu/rho*(cell[n].U[2][0] - cell[n-1].U[2][0])/dx;
            // North Face
            if(Top == 1){            // Farfield
                cell[n].Gv_n[0] = 0;
                cell[n].Gv_n[1] = myu/rho*(U_ini - cell[n].U[1][0])/(dy/2);
                cell[n].Gv_n[2] = myu/rho*(V_ini - cell[n].U[2][0])/(dy/2);
            }else if(Top == 2){      // Wall
                cell[n].Gv_n[0] = 0;
                cell[n].Gv_n[1] = myu/rho*(0 - cell[n].U[1][0])/(dy/2);
                cell[n].Gv_n[2] = myu/rho*(0 - cell[n].U[2][0])/(dy/2);
            }else if(Top == 3){
                cell[n].Gv_n[0] = 0;
                cell[n].Gv_n[1] = 0; 
                cell[n].Gv_n[2] = 0;
            }
            // South Face
            cell[n].Gv_s[0] = 0;
            cell[n].Gv_s[1] = myu/rho*(cell[n].U[1][0] - cell[n-cell_Nx].U[1][0])/dy;
            cell[n].Gv_s[2] = myu/rho*(cell[n].U[2][0] - cell[n-cell_Nx].U[2][0])/dy;
        }else{
            // East Face
            cell[n].Fv_e[0] = 0; 
            cell[n].Fv_e[1] = myu/rho*(cell[n+1].U[1][0] - cell[n].U[1][0])/dx;
            cell[n].Fv_e[2] = myu/rho*(cell[n+1].U[2][0] - cell[n].U[2][0])/dx;
            // West Face
            cell[n].Fv_w[0] = 0; 
            cell[n].Fv_w[1] = myu/rho*(cell[n].U[1][0] - cell[n-1].U[1][0])/dx;
            cell[n].Fv_w[2] = myu/rho*(cell[n].U[2][0] - cell[n-1].U[2][0])/dx;
            // North Face
            if(Top == 1){            // Farfield
                cell[n].Gv_n[0] = 0;
                cell[n].Gv_n[1] = myu/rho*(U_ini - cell[n].U[1][0])/(dy/2);
                cell[n].Gv_n[2] = myu/rho*(V_ini - cell[n].U[2][0])/(dy/2);
            }else if(Top == 2){      // Wall
                cell[n].Gv_n[0] = 0;
                cell[n].Gv_n[1] = myu/rho*(0 - cell[n].U[1][0])/(dy/2);
                cell[n].Gv_n[2] = myu/rho*(0 - cell[n].U[2][0])/(dy/2);
            }else if(Top == 3){
                cell[n].Gv_n[0] = 0;
                cell[n].Gv_n[1] = 0; 
                cell[n].Gv_n[2] = 0;
            }
            // South Face
            cell[n].Gv_s[0] = 0;
            cell[n].Gv_s[1] = myu/rho*(cell[n].U[1][0] - cell[n-cell_Nx].U[1][0])/dy;
            cell[n].Gv_s[2] = myu/rho*(cell[n].U[2][0] - cell[n-cell_Nx].U[2][0])/dy;
        } 
    } 
}


