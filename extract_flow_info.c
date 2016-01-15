#include "ibm.h"

void extract_cell_data()
{
    int i, j, n;

    // Flow info on Cartesian grid
    // Calculating u velocity
    for(j=1; j<=cell_Ny; j++){
        for(i=1; i<=cell_Nx; i++){
            n = (j-1)*cell_Nx + i;
            if(j==1){
                cell[n].u[0] = (cell[n+cell_Nx].Psi[1] - cell[n].Psi[1])/(cell[n+cell_Nx].y-cell[n].y);
            }else if(j==cell_Ny){
                cell[n].u[0] = (cell[n-cell_Nx].Psi[1] - cell[n].Psi[1])/(cell[n-cell_Nx].y-cell[n].y);
            }else{
                cell[n].u[0] = (cell[n+cell_Nx].Psi[1] - cell[n-cell_Nx].Psi[1])/(cell[n+cell_Nx].y-cell[n-cell_Nx].y);
            }
        }
    }
    //Calculating v velocity
    for(j=1; j<=cell_Ny; j++){
        for(i=1; i<=cell_Nx; i++){
            n = (j-1)*cell_Nx + i;
            if(i==1){
                cell[n].v[0] = -1*(cell[n+1].Psi[1] - cell[n].Psi[1])/(cell[n+1].x-cell[n].x);
            }else if(i==cell_Nx){
                cell[n].v[0] = -1*(cell[n-1].Psi[1] - cell[n].Psi[1])/(cell[n-1].x-cell[n].x);
            }else{
                cell[n].v[0] = -1*(cell[n+1].Psi[1] - cell[n-1].Psi[1])/(cell[n+1].x-cell[n-1].x);
            }
        }
    }
    //Calculating pressure values
    for(n=1; n<=Total_cells; n++){
        cell[n].p[0] = 0.5*(U_ini*U_ini - cell[n].u[0]*cell[n].u[0] - cell[n].v[0]*cell[n].v[0]);
    }
}

void extract_ib_data_johansen_colella()
{
    int n;

    // Flow info on Immersed Boundary surface
    //Calculating Cp values for immersed boundary surface 
    for(n=1; n<=N_ib; n++){
        if(cell_ib[n].d_cut != 0){
            cell_ib[n].Cp[0] = 1 - pow(cell_ib[n].wall_flux_diffusive[0]/(cell_ib[n].d_cut),2)/(U_ini*U_ini);
        }
    }

    //Calculating u and v values for immersed boundary surface
    for(n=1; n<=N_ib; n++){
        if(cell_ib[n].d_cut != 0){
            cell_ib[n].u[0] = -(cell_ib[n].wall_flux_diffusive[0]/cell_ib[n].d_cut)*sin(cell_ib[n].theta);
            cell_ib[n].v[0] = (cell_ib[n].wall_flux_diffusive[0]/cell_ib[n].d_cut)*cos(cell_ib[n].theta);
            cell_ib[n].p[0] = 0.5*(U_ini*U_ini - cell_ib[n].u[0]*cell_ib[n].u[0] - cell_ib[n].v[0]*cell_ib[n].v[0]);
        }
    }
}

void extract_ib_data_transfinite()
{
    int n;

    // Flow info on Immersed Boundary surface
    //Calculating Cp values for immersed boundary surface 
    for(n=1; n<=N_ib; n++){
        if(cell_ib[n].d_cut != 0){
            cell_ib[n].Cp[0] = 1 - pow(cell_ib[n].wall_flux_diffusive[0]/((cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])*Radius),2)/(U_ini*U_ini);
        }
    }

    //Calculating u and v values for immersed boundary surface
    for(n=1; n<=N_ib; n++){
        if(cell_ib[n].d_cut != 0){
            cell_ib[n].u[0] = -cell_ib[n].wall_flux_diffusive[0]/((cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])*Radius)*sin(cell_ib[n].theta);
            cell_ib[n].v[0] = cell_ib[n].wall_flux_diffusive[0]/((cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])*Radius)*cos(cell_ib[n].theta);
            cell_ib[n].p[0] = .5*(U_ini*U_ini - cell_ib[n].u[0]*cell_ib[n].u[0] - cell_ib[n].v[0]*cell_ib[n].v[0]);
        }
    }
}
