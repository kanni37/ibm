#include "functions.h"
#include "ibm.h"

void calculate_convective_flux(int eq_type, int Front, int Back, int Top, int Bottom)
{
    switch (eq_type){
        case 1:
            calculate_convective_flux_potential_flow();
            break;
        case 2:
            calculate_convective_flux_navier_stokes(Front, Back, Top, Bottom);
            break;
        default:
            printf("Error, bad eq_type input, quitting");
            break;
    }
}

void calculate_convective_flux_navier_stokes(int Front, int Back, int Top, int Bottom)
{
    int i, j, k, n, n_ib;
    double u_face, v_face, p_face;
    double pressure_term = 0;
    double pr_coeff = 0;

#pragma omp parallel for private(i, j, k, n, n_ib, u_face, v_face, p_face, pressure_term, pr_coeff) schedule(auto)
    //Interior Fluid and Solid Cells
    for(j=2; j<cell_Ny;j++){
        for(i=2; i<cell_Nx; i++){
            n = (j-1)*cell_Nx + i;
            pr_coeff = dx*dy/(4*myu);//10*pow(cell[n].dx*cell[n].dx + cell[n].dy*cell[n].dy, 0.5);
            if(cell[n].state != 2){
                // East Face 
                if(i == cell_Nx-1){
                    pressure_term = 0;
                }else{
                    pressure_term = pr_coeff*(cell[n+2].U[0][0] - 3*cell[n+1].U[0][0] + 3*cell[n].U[0][0] - cell[n-1].U[0][0]);
                }
                u_face = (cell[n+1].U[1][0] + cell[n].U[1][0])/2 + pressure_term; 
                v_face = (cell[n+1].U[2][0] + cell[n].U[2][0])/2;
                p_face = (cell[n+1].U[0][0] + cell[n].U[0][0])/2;
                cell[n].F_e[0] = u_face; 
                cell[n].F_e[1] = u_face*u_face + p_face/rho; 
                cell[n].F_e[2] = u_face*v_face; 
                // West Face
                if(i == 2){
                    pressure_term = 0;
                }else{
                    pressure_term = pr_coeff*(cell[n+1].U[0][0] - 3*cell[n].U[0][0] + 3*cell[n-1].U[0][0] - cell[n-2].U[0][0]);
                }
                u_face = (cell[n-1].U[1][0] + cell[n].U[1][0])/2 + pressure_term;
                v_face = (cell[n-1].U[2][0] + cell[n].U[2][0])/2;
                p_face = (cell[n-1].U[0][0] + cell[n].U[0][0])/2;
                cell[n].F_w[0] = u_face; 
                cell[n].F_w[1] = u_face*u_face + p_face/rho; 
                cell[n].F_w[2] = u_face*v_face; 
                // North Face
                if(j == cell_Ny-1){
                    pressure_term = 0;
                }else{
                    pressure_term = pr_coeff*(cell[n+2*cell_Nx].U[0][0] - 3*cell[n+1*cell_Nx].U[0][0] + 3*cell[n].U[0][0] - cell[n-1*cell_Nx].U[0][0]);
                }
                u_face = (cell[n+cell_Nx].U[1][0] + cell[n].U[1][0])/2;
                v_face = (cell[n+cell_Nx].U[2][0] + cell[n].U[2][0])/2 + pressure_term; 
                p_face = (cell[n+cell_Nx].U[0][0] + cell[n].U[0][0])/2;
                cell[n].G_n[0] = v_face; 
                cell[n].G_n[1] = u_face*v_face; 
                cell[n].G_n[2] = v_face*v_face + p_face/rho; 
                // South Face
                if(j == 2){
                    pressure_term = 0;
                }else{
                    pressure_term = pr_coeff*(cell[n+1*cell_Nx].U[0][0] - 3*cell[n].U[0][0] + 3*cell[n-1*cell_Nx].U[0][0] - cell[n-2*cell_Nx].U[0][0]);
                }
                u_face = (cell[n-cell_Nx].U[1][0] + cell[n].U[1][0])/2; 
                v_face = (cell[n-cell_Nx].U[2][0] + cell[n].U[2][0])/2 + pressure_term;
                p_face = (cell[n-cell_Nx].U[0][0] + cell[n].U[0][0])/2;
                cell[n].G_s[0] = v_face; 
                cell[n].G_s[1] = u_face*v_face; 
                cell[n].G_s[2] = v_face*v_face + p_face/rho; 
            }else if(cell[n].state == 2){
                n_ib = cell[n].ib_cell_no;
                // East Face
                k = 1;
                pressure_term = pr_coeff*(cell[n+2].U[0][0] - 3*cell[n+1].U[0][0] + 3*cell[n].U[0][0] - cell[n-1].U[0][0]);
                u_face = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[1][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[1][0])/2*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[1][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[1][0])/2*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]) + pressure_term; 
                v_face = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[2][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[2][0])/2*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[2][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[2][0])/2*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                p_face = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[0][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[0][0])/2*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[0][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[0][0])/2*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                cell[n].F_e[0] = u_face; 
                cell[n].F_e[1] = u_face*u_face + p_face/rho; 
                cell[n].F_e[2] = u_face*v_face; 
                // West Face
                k = 3;
                pressure_term = pr_coeff*(cell[n+1].U[0][0] - 3*cell[n].U[0][0] + 3*cell[n-1].U[0][0] - cell[n-2].U[0][0]);
                u_face = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[1][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[1][0])/2*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[1][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[1][0])/2*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]) + pressure_term; 
                v_face = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[2][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[2][0])/2*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[2][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[2][0])/2*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                p_face = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[0][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[0][0])/2*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[0][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[0][0])/2*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                cell[n].F_w[0] = u_face; 
                cell[n].F_w[1] = u_face*u_face + p_face/rho;
                cell[n].F_w[2] = u_face*v_face;
                // North Face 
                k = 2;
                pressure_term = pr_coeff*(cell[n+2*cell_Nx].U[0][0] - 3*cell[n+1*cell_Nx].U[0][0] + 3*cell[n].U[0][0] - cell[n-1*cell_Nx].U[0][0]);
                u_face = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[1][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[1][0])/2*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[1][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[1][0])/2*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                v_face = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[2][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[2][0])/2*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[2][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[2][0])/2*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]) + pressure_term; 
                p_face = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[0][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[0][0])/2*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[0][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[0][0])/2*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                cell[n].G_n[0] = v_face; 
                cell[n].G_n[1] = u_face*v_face; 
                cell[n].G_n[2] = v_face*v_face + p_face/rho; 
                // South Face
                k = 0;
                pressure_term = pr_coeff*(cell[n+1*cell_Nx].U[0][0] - 3*cell[n].U[0][0] + 3*cell[n-1*cell_Nx].U[0][0] - cell[n-2*cell_Nx].U[0][0]);
                u_face = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[1][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[1][0])/2*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[1][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[1][0])/2*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                v_face = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[2][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[2][0])/2*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[2][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[2][0])/2*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]) + pressure_term; 
                p_face = ((cell[cell_ib[n_ib].BI_flux_cell[k][2]].U[0][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][3]].U[0][0])/2*(cell_ib[n_ib].dy2[k]) + 
                        (cell[cell_ib[n_ib].BI_flux_cell[k][1]].U[0][0] +  cell[cell_ib[n_ib].BI_flux_cell[k][0]].U[0][0])/2*(cell_ib[n_ib].dy1[k]))/
                    (cell_ib[n_ib].dy1[k] + cell_ib[n_ib].dy2[k]); 
                cell[n].G_s[0] = v_face;
                cell[n].G_s[1] = u_face*v_face; 
                cell[n].G_s[2] = v_face*v_face + p_face/rho; 
            }
        }
    }

    //Inlet
    for(j=2; j<cell_Ny; j++){
        n = 1 + (j-1)*cell_Nx;
        // East Face 
        u_face = (cell[n+1].U[1][0] + cell[n].U[1][0])/2;
        v_face = (cell[n+1].U[2][0] + cell[n].U[2][0])/2;
        p_face = (cell[n+1].U[0][0] + cell[n].U[0][0])/2;
        cell[n].F_e[0] = u_face; 
        cell[n].F_e[1] = u_face*u_face + p_face/rho; 
        cell[n].F_e[2] = u_face*v_face; 
        // West Face
        u_face = U_ini;
        v_face = V_ini;
        p_face = cell[n].U[0][0];
        cell[n].F_w[0] = u_face; 
        cell[n].F_w[1] = u_face*u_face + p_face/rho; 
        cell[n].F_w[2] = u_face*v_face; 
        // North Face
        if(j == cell_Ny-1){
            pressure_term = 0;
        }else{
            pressure_term = pr_coeff*(cell[n+2*cell_Nx].U[0][0] - 3*cell[n+1*cell_Nx].U[0][0] + 3*cell[n].U[0][0] - cell[n-1*cell_Nx].U[0][0]);
        }
        u_face = (cell[n+cell_Nx].U[1][0] + cell[n].U[1][0])/2;
        v_face = (cell[n+cell_Nx].U[2][0] + cell[n].U[2][0])/2 + pressure_term; 
        p_face = (cell[n+cell_Nx].U[0][0] + cell[n].U[0][0])/2;
        cell[n].G_n[0] = v_face; 
        cell[n].G_n[1] = u_face*v_face; 
        cell[n].G_n[2] = v_face*v_face + p_face/rho; 
        // South Face
        if(j == 2){
            pressure_term = 0;
        }else{
            pressure_term = pr_coeff*(cell[n+1*cell_Nx].U[0][0] - 3*cell[n].U[0][0] + 3*cell[n-1*cell_Nx].U[0][0] - cell[n-2*cell_Nx].U[0][0]);
        }
        u_face = (cell[n-cell_Nx].U[1][0] + cell[n].U[1][0])/2; 
        v_face = (cell[n-cell_Nx].U[2][0] + cell[n].U[2][0])/2 + pressure_term;
        p_face = (cell[n-cell_Nx].U[0][0] + cell[n].U[0][0])/2;
        cell[n].G_s[0] = v_face; 
        cell[n].G_s[1] = u_face*v_face; 
        cell[n].G_s[2] = v_face*v_face + p_face/rho; 
    } 
    //Outlet
    for(j=2; j<cell_Ny; j++){
        n = j*cell_Nx;
        // East Face 
        u_face = cell[n].U[1][0]; 
        v_face = cell[n].U[2][0];
        p_face = P_out;
        cell[n].F_e[0] = u_face; 
        cell[n].F_e[1] = u_face*u_face + p_face/rho; 
        cell[n].F_e[2] = u_face*v_face; 
        // West Face
        u_face = (cell[n-1].U[1][0] + cell[n].U[1][0])/2;
        v_face = (cell[n-1].U[2][0] + cell[n].U[2][0])/2;
        p_face = (cell[n-1].U[0][0] + cell[n].U[0][0])/2;
        cell[n].F_w[0] = u_face; 
        cell[n].F_w[1] = u_face*u_face + p_face/rho; 
        cell[n].F_w[2] = u_face*v_face; 
        // North Face
        if(j == cell_Ny-1){
            pressure_term = 0;
        }else{
            pressure_term = pr_coeff*(cell[n+2*cell_Nx].U[0][0] - 3*cell[n+1*cell_Nx].U[0][0] + 3*cell[n].U[0][0] - cell[n-1*cell_Nx].U[0][0]);
        }
        u_face = (cell[n+cell_Nx].U[1][0] + cell[n].U[1][0])/2;
        v_face = (cell[n+cell_Nx].U[2][0] + cell[n].U[2][0])/2 + pressure_term; 
        p_face = (cell[n+cell_Nx].U[0][0] + cell[n].U[0][0])/2;
        cell[n].G_n[0] = v_face; 
        cell[n].G_n[1] = u_face*v_face; 
        cell[n].G_n[2] = v_face*v_face + p_face/rho; 
        // South Face
        if(j == 2){
            pressure_term = 0;
        }else{
            pressure_term = pr_coeff*(cell[n+1*cell_Nx].U[0][0] - 3*cell[n].U[0][0] + 3*cell[n-1*cell_Nx].U[0][0] - cell[n-2*cell_Nx].U[0][0]);
        }
        u_face = (cell[n-cell_Nx].U[1][0] + cell[n].U[1][0])/2; 
        v_face = (cell[n-cell_Nx].U[2][0] + cell[n].U[2][0])/2 + pressure_term;
        p_face = (cell[n-cell_Nx].U[0][0] + cell[n].U[0][0])/2;
        cell[n].G_s[0] = v_face; 
        cell[n].G_s[1] = u_face*v_face; 
        cell[n].G_s[2] = v_face*v_face + p_face/rho; 
    } 
    //Farfield/Wall Bottom
    for(i=1; i<=cell_Nx; i++){
        n = i;
        if(i == 1){
            // East Face
            u_face = (cell[n+1].U[1][0] + cell[n].U[1][0])/2;
            v_face = (cell[n+1].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n+1].U[0][0] + cell[n].U[0][0])/2;
            cell[n].F_e[0] = u_face; 
            cell[n].F_e[1] = u_face*u_face + p_face/rho; 
            cell[n].F_e[2] = u_face*v_face; 
            // West Face
            u_face = U_ini; 
            v_face = V_ini;
            p_face = cell[n].U[0][0];
            cell[n].F_w[0] = u_face; 
            cell[n].F_w[1] = u_face*u_face + p_face/rho; 
            cell[n].F_w[2] = u_face*v_face; 
            // North Face
            u_face = (cell[n+cell_Nx].U[1][0] + cell[n].U[1][0])/2;
            v_face = (cell[n+cell_Nx].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n+cell_Nx].U[0][0] + cell[n].U[0][0])/2;
            cell[n].G_n[0] = v_face; 
            cell[n].G_n[1] = u_face*v_face; 
            cell[n].G_n[2] = v_face*v_face + p_face/rho; 
            // South Face
            if(Bottom == 1){        // Farfield
                u_face = U_ini; 
                v_face = V_ini;
                p_face = P_out;
            }else if(Bottom == 2){    // Wall
                u_face = 0; 
                v_face = 0;
                p_face = cell[n].U[0][0];
            }else if(Bottom == 3){  // Symmetry
                u_face = cell[n].U[1][0]; 
                v_face = cell[n].U[2][0];
                p_face = cell[n].U[0][0];
            }
            cell[n].G_s[0] = v_face; 
            cell[n].G_s[1] = u_face*v_face; 
            cell[n].G_s[2] = v_face*v_face + p_face/rho; 
        }else if(i == cell_Nx){
            // East Face
            u_face = cell[n].U[1][0];
            v_face = cell[n].U[2][0];
            p_face = P_out;
            cell[n].F_e[0] = u_face; 
            cell[n].F_e[1] = u_face*u_face + p_face/rho; 
            cell[n].F_e[2] = u_face*v_face; 
            // West Face
            u_face = (cell[n-1].U[1][0] + cell[n].U[1][0])/2;
            v_face = (cell[n-1].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n-1].U[0][0] + cell[n].U[0][0])/2;
            cell[n].F_w[0] = u_face; 
            cell[n].F_w[1] = u_face*u_face + p_face/rho; 
            cell[n].F_w[2] = u_face*v_face; 
            // North Face
            u_face = (cell[n+cell_Nx].U[1][0] + cell[n].U[1][0])/2;
            v_face = (cell[n+cell_Nx].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n+cell_Nx].U[0][0] + cell[n].U[0][0])/2;
            cell[n].G_n[0] = v_face; 
            cell[n].G_n[1] = u_face*v_face; 
            cell[n].G_n[2] = v_face*v_face + p_face/rho; 
            // South Face
            if(Bottom == 1){        // Farfield
                u_face = U_ini; 
                v_face = V_ini;
                p_face = P_out;
            }else if(Bottom == 2){    // Wall
                u_face = 0; 
                v_face = 0;
                p_face = cell[n].U[0][0];
            }else if(Bottom == 3){  // Symmetry
                u_face = cell[n].U[1][0]; 
                v_face = cell[n].U[2][0];
                p_face = cell[n].U[0][0];
            }
            cell[n].G_s[0] = v_face; 
            cell[n].G_s[1] = u_face*v_face; 
            cell[n].G_s[2] = v_face*v_face + p_face/rho; 
        }else{
            // East Face 
            if(i == cell_Nx-1){
                pressure_term = 0;
            }else{
                pressure_term = pr_coeff*(cell[n+2].U[0][0] - 3*cell[n+1].U[0][0] + 3*cell[n].U[0][0] - cell[n-1].U[0][0]);
            }
            u_face = (cell[n+1].U[1][0] + cell[n].U[1][0])/2 + pressure_term; 
            v_face = (cell[n+1].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n+1].U[0][0] + cell[n].U[0][0])/2;
            cell[n].F_e[0] = u_face; 
            cell[n].F_e[1] = u_face*u_face + p_face/rho; 
            cell[n].F_e[2] = u_face*v_face; 
            // West Face
            if(i == 2){
                pressure_term = 0;
            }else{
                pressure_term = pr_coeff*(cell[n+1].U[0][0] - 3*cell[n].U[0][0] + 3*cell[n-1].U[0][0] - cell[n-2].U[0][0]);
            }
            u_face = (cell[n-1].U[1][0] + cell[n].U[1][0])/2 + pressure_term;
            v_face = (cell[n-1].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n-1].U[0][0] + cell[n].U[0][0])/2;
            cell[n].F_w[0] = u_face; 
            cell[n].F_w[1] = u_face*u_face + p_face/rho; 
            cell[n].F_w[2] = u_face*v_face; 
            // North Face
            u_face = (cell[n+cell_Nx].U[1][0] + cell[n].U[1][0])/2;
            v_face = (cell[n+cell_Nx].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n+cell_Nx].U[0][0] + cell[n].U[0][0])/2;
            cell[n].G_n[0] = v_face; 
            cell[n].G_n[1] = u_face*v_face; 
            cell[n].G_n[2] = v_face*v_face + p_face/rho; 
            // South Face
            if(Bottom == 1){    //Farfield
                u_face = U_ini; 
                v_face = V_ini;
                p_face = P_out;
            }else if(Bottom == 2){    // Wall
                u_face = 0; 
                v_face = 0;
                p_face = cell[n].U[0][0];
            }else if(Bottom == 3){  // Symmetry
                u_face = cell[n].U[1][0]; 
                v_face = cell[n].U[2][0];
                p_face = cell[n].U[0][0];
            }
            cell[n].G_s[0] = v_face; 
            cell[n].G_s[1] = u_face*v_face; 
            cell[n].G_s[2] = v_face*v_face + p_face/rho; 
        } 
    }
    //Farfield/Wall Top
    for(i=1; i<=cell_Nx; i++){
        n = (cell_Ny-1)*cell_Nx+i;
        if(i == 1){
            // East Face
            u_face = (cell[n+1].U[1][0] + cell[n].U[1][0])/2;
            v_face = (cell[n+1].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n+1].U[0][0] + cell[n].U[0][0])/2;
            cell[n].F_e[0] = u_face; 
            cell[n].F_e[1] = u_face*u_face + p_face/rho; 
            cell[n].F_e[2] = u_face*v_face; 
            // West Face
            u_face = U_ini;
            v_face = V_ini;
            p_face = cell[n].U[0][0];
            cell[n].F_w[0] = u_face; 
            cell[n].F_w[1] = u_face*u_face + p_face/rho; 
            cell[n].F_w[2] = u_face*v_face; 
            if(Top == 1){        // Farfield
                u_face = U_ini; 
                v_face = V_ini;
                p_face = P_out;
            }else if(Top == 2){    // Wall
                u_face = 0; 
                v_face = 0;
                p_face = cell[n].U[0][0];
            }else if(Top == 3){     //Symmetry
                u_face = cell[n].U[1][0]; 
                v_face = cell[n].U[2][0];
                p_face = cell[n].U[0][0];
            }
            cell[n].G_n[0] = v_face; 
            cell[n].G_n[1] = u_face*v_face; 
            cell[n].G_n[2] = v_face*v_face + p_face/rho; 
            // South Face
            u_face = (cell[n-cell_Nx].U[1][0] + cell[n].U[1][0])/2;
            v_face = (cell[n-cell_Nx].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n-cell_Nx].U[0][0] + cell[n].U[0][0])/2;
            cell[n].G_s[0] = v_face; 
            cell[n].G_s[1] = u_face*v_face; 
            cell[n].G_s[2] = v_face*v_face + p_face/rho; 
        }else if(i == cell_Nx){
            // East Face
            u_face = cell[n].U[1][0]; 
            v_face = cell[n].U[2][0];
            p_face = P_out;
            cell[n].F_e[0] = u_face; 
            cell[n].F_e[1] = u_face*u_face + p_face/rho; 
            cell[n].F_e[2] = u_face*v_face; 
            // West Face
            u_face = (cell[n-1].U[1][0] + cell[n].U[1][0])/2;
            v_face = (cell[n-1].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n-1].U[0][0] + cell[n].U[0][0])/2;
            cell[n].F_w[0] = u_face; 
            cell[n].F_w[1] = u_face*u_face + p_face/rho; 
            cell[n].F_w[2] = u_face*v_face; 
            // North Face
            if(Top == 1){        // Farfield
                u_face = U_ini; 
                v_face = V_ini;
                p_face = P_out;
            }else if(Top == 2){    // Wall
                u_face = 0; 
                v_face = 0;
                p_face = cell[n].U[0][0];
            }else if(Top == 3){     //Symmetry
                u_face = cell[n].U[1][0]; 
                v_face = cell[n].U[2][0];
                p_face = cell[n].U[0][0];
            }
            cell[n].G_n[0] = v_face; 
            cell[n].G_n[1] = u_face*v_face; 
            cell[n].G_n[2] = v_face*v_face + p_face/rho; 
            // South Face
            u_face = (cell[n-cell_Nx].U[1][0] + cell[n].U[1][0])/2;
            v_face = (cell[n-cell_Nx].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n-cell_Nx].U[0][0] + cell[n].U[0][0])/2;
            cell[n].G_s[0] = v_face; 
            cell[n].G_s[1] = u_face*v_face; 
            cell[n].G_s[2] = v_face*v_face + p_face/rho; 
        }else{
            // East Face 
            if(i == cell_Nx-1){
                pressure_term = 0;
            }else{
                pressure_term = pr_coeff*(cell[n+2].U[0][0] - 3*cell[n+1].U[0][0] + 3*cell[n].U[0][0] - cell[n-1].U[0][0]);
            }
            u_face = (cell[n+1].U[1][0] + cell[n].U[1][0])/2 + pressure_term; 
            v_face = (cell[n+1].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n+1].U[0][0] + cell[n].U[0][0])/2;
            cell[n].F_e[0] = u_face; 
            cell[n].F_e[1] = u_face*u_face + p_face/rho; 
            cell[n].F_e[2] = u_face*v_face; 
            // West Face
            if(i == 2){
                pressure_term = 0;
            }else{
                pressure_term = pr_coeff*(cell[n+1].U[0][0] - 3*cell[n].U[0][0] + 3*cell[n-1].U[0][0] - cell[n-2].U[0][0]);
            }
            u_face = (cell[n-1].U[1][0] + cell[n].U[1][0])/2 + pressure_term;
            v_face = (cell[n-1].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n-1].U[0][0] + cell[n].U[0][0])/2;
            cell[n].F_w[0] = u_face; 
            cell[n].F_w[1] = u_face*u_face + p_face/rho; 
            cell[n].F_w[2] = u_face*v_face; 
            // North Face
            if(Top == 1){        // Farfield
                u_face = U_ini; 
                v_face = V_ini;
                p_face = P_out;
            }else if(Top == 2){    // Wall
                u_face = 0; 
                v_face = 0;
                p_face = cell[n].U[0][0];
            }else if(Top == 3){     //Symmetry
                u_face = cell[n].U[1][0]; 
                v_face = cell[n].U[2][0];
                p_face = cell[n].U[0][0];
            }
            cell[n].G_n[0] = v_face; 
            cell[n].G_n[1] = u_face*v_face; 
            cell[n].G_n[2] = v_face*v_face + p_face/rho; 
            // South Face
            u_face = (cell[n-cell_Nx].U[1][0] + cell[n].U[1][0])/2;
            v_face = (cell[n-cell_Nx].U[2][0] + cell[n].U[2][0])/2;
            p_face = (cell[n-cell_Nx].U[0][0] + cell[n].U[0][0])/2;
            cell[n].G_s[0] = v_face; 
            cell[n].G_s[1] = u_face*v_face; 
            cell[n].G_s[2] = v_face*v_face + p_face/rho; 
        } 
    } 
}


void calculate_convective_flux_potential_flow()
{
    // No calculation required since no convective flux for potential flow
    /*    int i, j, n;

#pragma omp parallel for private(i, j, n) schedule(auto)
    //Convective fluxes are all zero for potential flow
    for(j=1; j<=cell_Ny;j++){
    for(i=1; i<=cell_Nx; i++){
    n = (j-1)*cell_Nx + i;
    // East Face`
    cell[n].F_e[0] = 0;
    // West Face
    cell[n].F_w[0] = 0;
    // North Face
    cell[n].G_n[0] = 0;
    // South Face
    cell[n].G_s[0] = 0;
    }
    }
    */
}
