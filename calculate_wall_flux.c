#include "ibm.h"
#include "functions.h"

// Boundary conditions for Potential Flow
void calculate_wall_flux(int eq_type, int M_type) 
{
    switch(eq_type){
        case 1:
            switch (M_type){
                case 1:
                    wall_flux_potential_flow_transfinite();
                    break;
                case 2:
                    wall_flux_potential_flow_johansen_colella();
                    break;
                case 3:
                    wall_flux_potential_flow_johansen_colella();
                    break;
                default:
                    printf("Error, bad input, quitting");
                    break;
            }
            break;
        case 2:
            switch (M_type){
                case 1:
                    wall_flux_navier_stokes_transfinite();
                    break;
                case 2:
                    wall_flux_navier_stokes_johansen_colella();
                    break;
                case 3:
                    wall_flux_navier_stokes_johansen_colella();
                    break;
                default:
                    printf("Error, bad input, quitting");
                    break;
            }
            break;
        default:
            printf("Error, bad input, quitting");
            break;
    }
}

void wall_flux_potential_flow_transfinite()
{
    int k, n;
    double zi_1, zi_2, zi_m;
    double Psi_x, Psi_y;

    //For IB cells
    for(n=1; n<=N_ib; n++){
        if(cell_ib[n].d_cut == 0){
            cell_ib[n].wall_flux_diffusive[0] = 0;
        }else{
            cell_ib[n].wall_flux_diffusive[0] = 0;
            
            cell_ib[n].Psi_1 =  ((cell[cell_ib[n].BI_cell_1[3]].Psi[0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[2]].Psi[0]*cell_ib[n].dx2_1)*cell_ib[n].dy2_1 + 
                    (cell[cell_ib[n].BI_cell_1[0]].Psi[0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[1]].Psi[0]*cell_ib[n].dx2_1)*cell_ib[n].dy1_1)/
                ((cell_ib[n].dx1_1 + cell_ib[n].dx2_1)*(cell_ib[n].dy1_1 + cell_ib[n].dy2_1));
            cell_ib[n].Psi_2 =  ((cell[cell_ib[n].BI_cell_2[3]].Psi[0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[2]].Psi[0]*cell_ib[n].dx2_2)*cell_ib[n].dy2_2 + 
                    (cell[cell_ib[n].BI_cell_2[0]].Psi[0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[1]].Psi[0]*cell_ib[n].dx2_2)*cell_ib[n].dy1_2)/
                ((cell_ib[n].dx1_2 + cell_ib[n].dx2_2)*(cell_ib[n].dy1_2 + cell_ib[n].dy2_2));
            cell_ib[n].Psi_3 =  ((cell[cell_ib[n].BI_cell_3[3]].Psi[0]*cell_ib[n].dx1_3 + cell[cell_ib[n].BI_cell_3[2]].Psi[0]*cell_ib[n].dx2_3)*cell_ib[n].dy2_3 + 
                    (cell[cell_ib[n].BI_cell_3[0]].Psi[0]*cell_ib[n].dx1_3 + cell[cell_ib[n].BI_cell_3[1]].Psi[0]*cell_ib[n].dx2_3)*cell_ib[n].dy1_3)/
                ((cell_ib[n].dx1_3 + cell_ib[n].dx2_3)*(cell_ib[n].dy1_3 + cell_ib[n].dy2_3));
            cell_ib[n].Psi_4 =  ((cell[cell_ib[n].BI_cell_4[3]].Psi[0]*cell_ib[n].dx1_4 + cell[cell_ib[n].BI_cell_4[2]].Psi[0]*cell_ib[n].dx2_4)*cell_ib[n].dy2_4 + 
                    (cell[cell_ib[n].BI_cell_4[0]].Psi[0]*cell_ib[n].dx1_4 + cell[cell_ib[n].BI_cell_4[1]].Psi[0]*cell_ib[n].dx2_4)*cell_ib[n].dy1_4)/
                ((cell_ib[n].dx1_4 + cell_ib[n].dx2_4)*(cell_ib[n].dy1_4 + cell_ib[n].dy2_4));

            //printf("Psi_1 = %lf, Psi_2 = %lf, Psi_3 = %lf, Psi_4 = %lf\n",cell_ib[n].Psi_1, cell_ib[n].Psi_2, cell_ib[n].Psi_3, cell_ib[n].Psi_4);

            //-------------------------------Calculating flux through wall----------------------------------
            for(k=0; k<cell_ib[n].N_local_division; k++){
                zi_1 = 1.0*k/cell_ib[n].N_local_division;
                zi_2 = zi_1 + 1.0/cell_ib[n].N_local_division;
                zi_m = (zi_1 + zi_2)/2;

                Psi_x = cell_ib[n].eta_x[k]*(-Psi0*1.5 + 2.0*(cell_ib[n].Psi_1 + zi_m*(cell_ib[n].Psi_2 - cell_ib[n].Psi_1)) - 
                        0.5*(cell_ib[n].Psi_3 + zi_m*(cell_ib[n].Psi_4 - cell_ib[n].Psi_3)));
                Psi_y = cell_ib[n].eta_y[k]*(-Psi0*1.5 + 2.0*(cell_ib[n].Psi_1 + zi_m*(cell_ib[n].Psi_2 - cell_ib[n].Psi_1)) - 
                        0.5*(cell_ib[n].Psi_3 + zi_m*(cell_ib[n].Psi_4 - cell_ib[n].Psi_3)));
                //printf("Psi_x = %lf, Psi_y = %lf\n",Psi_x, Psi_y);

                cell_ib[n].wall_flux_diffusive[0] = cell_ib[n].wall_flux_diffusive[0] + (Psi_x*(-cell_ib[n].del_y[k]) + Psi_y*(cell_ib[n].del_x[k]));
            }
        }
    }
}

void wall_flux_potential_flow_johansen_colella()
{
    int n;
    double Psi_1, Psi_2;
    double d1, d2;

    //For IB cells
    for(n=1; n<=N_ib; n++){
        if(cell_ib[n].d_cut == 0){
            cell_ib[n].wall_flux_diffusive[0] = 0;
        }else{
            cell_ib[n].wall_flux_diffusive[0] = 0;
            
            Psi_1 =    ((cell[cell_ib[n].BI_cell_1[3]].Psi[0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[2]].Psi[0]*cell_ib[n].dx2_1)*cell_ib[n].dy2_1 + 
                    (cell[cell_ib[n].BI_cell_1[0]].Psi[0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[1]].Psi[0]*cell_ib[n].dx2_1)*cell_ib[n].dy1_1)/
                ((cell_ib[n].dy1_1 + cell_ib[n].dy2_1)*(cell_ib[n].dx1_1 + cell_ib[n].dx2_1));

            Psi_2 =    ((cell[cell_ib[n].BI_cell_2[3]].Psi[0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[2]].Psi[0]*cell_ib[n].dx2_2)*cell_ib[n].dy2_2 + 
                    (cell[cell_ib[n].BI_cell_2[0]].Psi[0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[1]].Psi[0]*cell_ib[n].dx2_2)*cell_ib[n].dy1_2)/
                ((cell_ib[n].dy1_2 + cell_ib[n].dy2_2)*(cell_ib[n].dx1_2 + cell_ib[n].dx2_2));

            d1 = cell_ib[n].d1; 
            d2 = cell_ib[n].d2;

            cell_ib[n].wall_flux_diffusive[0] = (d2/d1*(Psi0 - Psi_1) - d1/d2*(Psi0 - Psi_2))/(d2-d1)*cell_ib[n].d_cut;                  // For non-constant grid spacing    
        }
    }
}

void wall_flux_navier_stokes_transfinite()
{
    int k, n;
    double zi_1, zi_2, zi_m;
    double del_x, del_y;
    double U1_x, U1_y;
    double U2_x, U2_y;
    double P_wall;
    double x_zi, y_zi, x_eta, y_eta;
    double Jacob, eta_x, eta_y;

    //Convective Wall Fluxes
    for(n=1; n<=N_ib; n++){
        cell_ib[n].wall_flux_convective[0] = 0;
        cell_ib[n].wall_flux_convective[1] = 0;
        cell_ib[n].wall_flux_convective[2] = 0;

        cell_ib[n].U0_1 =  ((cell[cell_ib[n].BI_cell_1[3]].U[0][0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[2]].U[0][0]*cell_ib[n].dx2_1)*cell_ib[n].dy2_1 + 
                (cell[cell_ib[n].BI_cell_1[0]].U[0][0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[1]].U[0][0]*cell_ib[n].dx2_1)*cell_ib[n].dy1_1)/
            ((cell_ib[n].dx1_1 + cell_ib[n].dx2_1)*(cell_ib[n].dy1_1 + cell_ib[n].dy2_1));
        cell_ib[n].U0_2 =  ((cell[cell_ib[n].BI_cell_2[3]].U[0][0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[2]].U[0][0]*cell_ib[n].dx2_2)*cell_ib[n].dy2_2 + 
                (cell[cell_ib[n].BI_cell_2[0]].U[0][0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[1]].U[0][0]*cell_ib[n].dx2_2)*cell_ib[n].dy1_2)/
            ((cell_ib[n].dx1_2 + cell_ib[n].dx2_2)*(cell_ib[n].dy1_2 + cell_ib[n].dy2_2));
        cell_ib[n].U0_3 =  ((cell[cell_ib[n].BI_cell_3[3]].U[0][0]*cell_ib[n].dx1_3 + cell[cell_ib[n].BI_cell_3[2]].U[0][0]*cell_ib[n].dx2_3)*cell_ib[n].dy2_3 + 
                (cell[cell_ib[n].BI_cell_3[0]].U[0][0]*cell_ib[n].dx1_3 + cell[cell_ib[n].BI_cell_3[1]].U[0][0]*cell_ib[n].dx2_3)*cell_ib[n].dy1_3)/
            ((cell_ib[n].dx1_3 + cell_ib[n].dx2_3)*(cell_ib[n].dy1_3 + cell_ib[n].dy2_3));
        cell_ib[n].U0_4 =  ((cell[cell_ib[n].BI_cell_4[3]].U[0][0]*cell_ib[n].dx1_4 + cell[cell_ib[n].BI_cell_4[2]].U[0][0]*cell_ib[n].dx2_4)*cell_ib[n].dy2_4 + 
                (cell[cell_ib[n].BI_cell_4[0]].U[0][0]*cell_ib[n].dx1_4 + cell[cell_ib[n].BI_cell_4[1]].U[0][0]*cell_ib[n].dx2_4)*cell_ib[n].dy1_4)/
            ((cell_ib[n].dx1_4 + cell_ib[n].dx2_4)*(cell_ib[n].dy1_4 + cell_ib[n].dy2_4));

        //printf("U[0]_1 = %lf, U[0]_2 = %lf, U[0]_3 = %lf, U[0]_4 = %lf\n",cell_ib[n].U[0]_1, cell_ib[n].U[0]_2, cell_ib[n].U[0]_3, cell_ib[n].U[0]_4);

        //-------------------------------Calculating flux through wall----------------------------------
        for(k=0; k<cell_ib[n].N_local_division; k++){
            zi_1 = 1.0*k/cell_ib[n].N_local_division;
            zi_2 = zi_1 + 1.0/cell_ib[n].N_local_division;
            zi_m = (zi_1 + zi_2)/2;

            del_x = Radius*(cos(cell_ib[n].theta_cut[0] + zi_2*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])) 
                    - cos(cell_ib[n].theta_cut[0] + zi_1*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])));
            del_y = Radius*(sin(cell_ib[n].theta_cut[0] + zi_2*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])) 
                    - sin(cell_ib[n].theta_cut[0] + zi_1*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])));
            //printf("del_x = %lf, del_y = %lf\n",del_x, del_y);

            x_zi  = - Radius*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])*sin(cell_ib[n].theta_cut[0] + 
                    zi_m*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0]));
            x_eta = - 1.5*(xo + Radius*cos(cell_ib[n].theta_cut[0] + zi_m*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0]))) + 
                2.0*(cell_ib[n].x_1 + zi_m*(cell_ib[n].x_2 - cell_ib[n].x_1)) - 
                0.5*(cell_ib[n].x_3 + zi_m*(cell_ib[n].x_4 - cell_ib[n].x_3));
            y_zi  =   Radius*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])*cos(cell_ib[n].theta_cut[0] + 
                    zi_m*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0]));
            y_eta = - 1.5*(yo + Radius*sin(cell_ib[n].theta_cut[0] + zi_m*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0]))) + 
                2.0*(cell_ib[n].y_1 + zi_m*(cell_ib[n].y_2 - cell_ib[n].y_1)) - 
                0.5*(cell_ib[n].y_3 + zi_m*(cell_ib[n].y_4 - cell_ib[n].y_3));
            //printf("x_zi = %lf, x_eta = %lf, y_zi = %lf, y_eta = %lf\n",x_zi, x_eta, y_zi, y_eta);

            Jacob = x_zi*y_eta - y_zi*x_eta;
            eta_x = -y_zi/Jacob;
            eta_y = x_zi/Jacob;

            // Need normal pressure derivative to be zero. Finding Normal pressure by multiplying by surface area vector with the pressure gradient.
            // Finding pressure on surface by equating the normal pressure to zero.
            // Psi_x*-dy + Psi_y*dx = 0

            P_wall = 1.0/((-del_y)*eta_x*1.5 + del_x*eta_y*1.5)*(-del_y)*eta_x*(2.0*(cell_ib[n].U0_1 + zi_m*(cell_ib[n].U0_2 - cell_ib[n].U0_1)) - 
                    0.5*(cell_ib[n].U0_3 + zi_m*(cell_ib[n].U0_4 - cell_ib[n].U0_3))) + del_x*eta_y*(2.0*(cell_ib[n].U0_1 + zi_m*(cell_ib[n].U0_2 - cell_ib[n].U0_1)) - 
                    0.5*(cell_ib[n].U0_3 + zi_m*(cell_ib[n].U0_4 - cell_ib[n].U0_3)));

            //printf("Psi_x = %lf, Psi_y = %lf\n",Psi_x, Psi_y);
            cell_ib[n].wall_flux_convective[0] = 0; 
            cell_ib[n].wall_flux_convective[1] = cell_ib[n].wall_flux_convective[1] + P_wall/rho*(-del_y); 
            cell_ib[n].wall_flux_convective[2] = cell_ib[n].wall_flux_convective[2] + P_wall/rho*(del_x);
        }
    }

    // Diffusive Wall Fluxes
    for(n=1; n<=N_ib; n++){
        cell_ib[n].wall_flux_diffusive[0] = 0;
        cell_ib[n].wall_flux_diffusive[1] = 0;
        cell_ib[n].wall_flux_diffusive[2] = 0;

        cell_ib[n].U1_1 =  ((cell[cell_ib[n].BI_cell_1[3]].U[1][0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[2]].U[1][0]*cell_ib[n].dx2_1)*cell_ib[n].dy2_1 + 
                (cell[cell_ib[n].BI_cell_1[0]].U[1][0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[1]].U[1][0]*cell_ib[n].dx2_1)*cell_ib[n].dy1_1)/
            ((cell_ib[n].dx1_1 + cell_ib[n].dx2_1)*(cell_ib[n].dy1_1 + cell_ib[n].dy2_1));
        cell_ib[n].U1_2 =  ((cell[cell_ib[n].BI_cell_2[3]].U[1][0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[2]].U[1][0]*cell_ib[n].dx2_2)*cell_ib[n].dy2_2 + 
                (cell[cell_ib[n].BI_cell_2[0]].U[1][0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[1]].U[1][0]*cell_ib[n].dx2_2)*cell_ib[n].dy1_2)/
            ((cell_ib[n].dx1_2 + cell_ib[n].dx2_2)*(cell_ib[n].dy1_2 + cell_ib[n].dy2_2));
        cell_ib[n].U1_3 =  ((cell[cell_ib[n].BI_cell_3[3]].U[1][0]*cell_ib[n].dx1_3 + cell[cell_ib[n].BI_cell_3[2]].U[1][0]*cell_ib[n].dx2_3)*cell_ib[n].dy2_3 + 
                (cell[cell_ib[n].BI_cell_3[0]].U[1][0]*cell_ib[n].dx1_3 + cell[cell_ib[n].BI_cell_3[1]].U[1][0]*cell_ib[n].dx2_3)*cell_ib[n].dy1_3)/
            ((cell_ib[n].dx1_3 + cell_ib[n].dx2_3)*(cell_ib[n].dy1_3 + cell_ib[n].dy2_3));
        cell_ib[n].U1_4 =  ((cell[cell_ib[n].BI_cell_4[3]].U[1][0]*cell_ib[n].dx1_4 + cell[cell_ib[n].BI_cell_4[2]].U[1][0]*cell_ib[n].dx2_4)*cell_ib[n].dy2_4 + 
                (cell[cell_ib[n].BI_cell_4[0]].U[1][0]*cell_ib[n].dx1_4 + cell[cell_ib[n].BI_cell_4[1]].U[1][0]*cell_ib[n].dx2_4)*cell_ib[n].dy1_4)/
            ((cell_ib[n].dx1_4 + cell_ib[n].dx2_4)*(cell_ib[n].dy1_4 + cell_ib[n].dy2_4));

        cell_ib[n].U2_1 =  ((cell[cell_ib[n].BI_cell_1[3]].U[2][0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[2]].U[2][0]*cell_ib[n].dx2_1)*cell_ib[n].dy2_1 + 
                (cell[cell_ib[n].BI_cell_1[0]].U[2][0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[1]].U[2][0]*cell_ib[n].dx2_1)*cell_ib[n].dy1_1)/
            ((cell_ib[n].dx1_1 + cell_ib[n].dx2_1)*(cell_ib[n].dy1_1 + cell_ib[n].dy2_1));
        cell_ib[n].U2_2 =  ((cell[cell_ib[n].BI_cell_2[3]].U[2][0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[2]].U[2][0]*cell_ib[n].dx2_2)*cell_ib[n].dy2_2 + 
                (cell[cell_ib[n].BI_cell_2[0]].U[2][0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[1]].U[2][0]*cell_ib[n].dx2_2)*cell_ib[n].dy1_2)/
            ((cell_ib[n].dx1_2 + cell_ib[n].dx2_2)*(cell_ib[n].dy1_2 + cell_ib[n].dy2_2));
        cell_ib[n].U2_3 =  ((cell[cell_ib[n].BI_cell_3[3]].U[2][0]*cell_ib[n].dx1_3 + cell[cell_ib[n].BI_cell_3[2]].U[2][0]*cell_ib[n].dx2_3)*cell_ib[n].dy2_3 + 
                (cell[cell_ib[n].BI_cell_3[0]].U[2][0]*cell_ib[n].dx1_3 + cell[cell_ib[n].BI_cell_3[1]].U[2][0]*cell_ib[n].dx2_3)*cell_ib[n].dy1_3)/
            ((cell_ib[n].dx1_3 + cell_ib[n].dx2_3)*(cell_ib[n].dy1_3 + cell_ib[n].dy2_3));
        cell_ib[n].U2_4 =  ((cell[cell_ib[n].BI_cell_4[3]].U[2][0]*cell_ib[n].dx1_4 + cell[cell_ib[n].BI_cell_4[2]].U[2][0]*cell_ib[n].dx2_4)*cell_ib[n].dy2_4 + 
                (cell[cell_ib[n].BI_cell_4[0]].U[2][0]*cell_ib[n].dx1_4 + cell[cell_ib[n].BI_cell_4[1]].U[2][0]*cell_ib[n].dx2_4)*cell_ib[n].dy1_4)/
            ((cell_ib[n].dx1_4 + cell_ib[n].dx2_4)*(cell_ib[n].dy1_4 + cell_ib[n].dy2_4));

        //printf("U[0]_1 = %lf, U[0]_2 = %lf, U[0]_3 = %lf, U[0]_4 = %lf\n",cell_ib[n].U[0]_1, cell_ib[n].U[0]_2, cell_ib[n].U[0]_3, cell_ib[n].U[0]_4);

        //-------------------------------Calculating flux through wall----------------------------------
        for(k=0; k<cell_ib[n].N_local_division; k++){
            zi_1 = 1.0*k/cell_ib[n].N_local_division;
            zi_2 = zi_1 + 1.0/cell_ib[n].N_local_division;
            zi_m = (zi_1 + zi_2)/2;

            del_x = Radius*(cos(cell_ib[n].theta_cut[0] + zi_2*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])) 
                    - cos(cell_ib[n].theta_cut[0] + zi_1*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])));
            del_y = Radius*(sin(cell_ib[n].theta_cut[0] + zi_2*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])) 
                    - sin(cell_ib[n].theta_cut[0] + zi_1*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])));
            //printf("del_x = %lf, del_y = %lf\n",del_x, del_y);

            x_zi  = - Radius*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])*sin(cell_ib[n].theta_cut[0] + 
                    zi_m*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0]));
            x_eta = - 1.5*(xo + Radius*cos(cell_ib[n].theta_cut[0] + zi_m*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0]))) + 
                2.0*(cell_ib[n].x_1 + zi_m*(cell_ib[n].x_2 - cell_ib[n].x_1)) - 
                0.5*(cell_ib[n].x_3 + zi_m*(cell_ib[n].x_4 - cell_ib[n].x_3));
            y_zi  =   Radius*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])*cos(cell_ib[n].theta_cut[0] + 
                    zi_m*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0]));
            y_eta = - 1.5*(yo + Radius*sin(cell_ib[n].theta_cut[0] + zi_m*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0]))) + 
                2.0*(cell_ib[n].y_1 + zi_m*(cell_ib[n].y_2 - cell_ib[n].y_1)) - 
                0.5*(cell_ib[n].y_3 + zi_m*(cell_ib[n].y_4 - cell_ib[n].y_3));
            //printf("x_zi = %lf, x_eta = %lf, y_zi = %lf, y_eta = %lf\n",x_zi, x_eta, y_zi, y_eta);

            Jacob = x_zi*y_eta - y_zi*x_eta;
            eta_x = -y_zi/Jacob;
            eta_y = x_zi/Jacob;

            // Need normal pressure derivative to be zero. Finding Normal pressure by multiplying by surface area vector with the pressure gradient.
            // Finding pressure on surface by equating the normal pressure to zero.
            // Psi_x*-dy + Psi_y*dx = 0

            U1_x = eta_x*(0*1.5 + 2.0*(cell_ib[n].U1_1 + zi_m*(cell_ib[n].U1_2 - cell_ib[n].U1_1)) - 
                    0.5*(cell_ib[n].U1_3 + zi_m*(cell_ib[n].U1_4 - cell_ib[n].U1_3)));
            U1_y = eta_y*(0*1.5 + 2.0*(cell_ib[n].U1_1 + zi_m*(cell_ib[n].U1_2 - cell_ib[n].U1_1)) - 
                    0.5*(cell_ib[n].U1_3 + zi_m*(cell_ib[n].U1_4 - cell_ib[n].U1_3)));
            U2_x = eta_x*(0*1.5 + 2.0*(cell_ib[n].U2_1 + zi_m*(cell_ib[n].U2_2 - cell_ib[n].U2_1)) - 
                    0.5*(cell_ib[n].U2_3 + zi_m*(cell_ib[n].U2_4 - cell_ib[n].U2_3)));
            U2_y = eta_y*(0*1.5 + 2.0*(cell_ib[n].U2_1 + zi_m*(cell_ib[n].U2_2 - cell_ib[n].U2_1)) - 
                    0.5*(cell_ib[n].U2_3 + zi_m*(cell_ib[n].U2_4 - cell_ib[n].U2_3)));

            cell_ib[n].wall_flux_diffusive[0] = 0; 
            cell_ib[n].wall_flux_diffusive[1] = cell_ib[n].wall_flux_diffusive[1] + myu/rho*(U1_x*(-del_y) + U1_y*(del_x));
            cell_ib[n].wall_flux_diffusive[2] = cell_ib[n].wall_flux_diffusive[2] + myu/rho*(U2_x*(-del_y) + U2_y*(del_x));
        }
    }
}

void wall_flux_navier_stokes_johansen_colella()
{
    int n;
    double U0_1, U0_2;
    double U1_1, U1_2;
    double U2_1, U2_2;
    double d1, d2;


    //Convective Wall Fluxes
    for(n=1; n<=N_ib; n++){
        U0_1 =    ((cell[cell_ib[n].BI_cell_1[3]].U[0][1]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[2]].U[0][1]*cell_ib[n].dx2_1)*cell_ib[n].dy2_1 + 
                (cell[cell_ib[n].BI_cell_1[0]].U[0][1]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[1]].U[0][1]*cell_ib[n].dx2_1)*cell_ib[n].dy1_1)/
            ((cell_ib[n].dy1_1 + cell_ib[n].dy2_1)*(cell_ib[n].dx1_1 + cell_ib[n].dx2_1));

        U0_2 =    ((cell[cell_ib[n].BI_cell_2[3]].U[0][1]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[2]].U[0][1]*cell_ib[n].dx2_2)*cell_ib[n].dy2_2 + 
                (cell[cell_ib[n].BI_cell_2[0]].U[0][1]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[1]].U[0][1]*cell_ib[n].dx2_2)*cell_ib[n].dy1_2)/
            ((cell_ib[n].dy1_2 + cell_ib[n].dy2_2)*(cell_ib[n].dx1_2 + cell_ib[n].dx2_2));

        d1 = cell_ib[n].d1; 
        d2 = cell_ib[n].d2;

        cell_ib[n].wall_flux_convective[0] = 0;                                                         // u and v velocity on wall zero 
        cell_ib[n].wall_flux_convective[1] = (d2/d1*U0_1 - d1/d2*U0_2)/(d2/d1 - d1/d2)/rho*(cell_ib[n].y_cut[1] - cell_ib[n].y_cut[0]);        // p/rho or u0 found by equating normal derivative on wall to zero 
        cell_ib[n].wall_flux_convective[2] = -(d2/d1*U0_1 - d1/d2*U0_2)/(d2/d1 - d1/d2)/rho*(cell_ib[n].x_cut[1] - cell_ib[n].x_cut[0]);        // p/rho or u0 found by equating normal derivative on wall to zero 
    }

    //Diffusive Wall Fluxes
    for(n=1; n<=N_ib; n++){
        U1_1 =    ((cell[cell_ib[n].BI_cell_1[3]].U[1][0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[2]].U[1][0]*cell_ib[n].dx2_1)*cell_ib[n].dy2_1 + 
                (cell[cell_ib[n].BI_cell_1[0]].U[1][0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[1]].U[1][0]*cell_ib[n].dx2_1)*cell_ib[n].dy1_1)/
            ((cell_ib[n].dy1_1 + cell_ib[n].dy2_1)*(cell_ib[n].dx1_1 + cell_ib[n].dx2_1));

        U1_2 =    ((cell[cell_ib[n].BI_cell_2[3]].U[1][0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[2]].U[1][0]*cell_ib[n].dx2_2)*cell_ib[n].dy2_2 + 
                (cell[cell_ib[n].BI_cell_2[0]].U[1][0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[1]].U[1][0]*cell_ib[n].dx2_2)*cell_ib[n].dy1_2)/
            ((cell_ib[n].dy1_2 + cell_ib[n].dy2_2)*(cell_ib[n].dx1_2 + cell_ib[n].dx2_2));

        U2_1 =    ((cell[cell_ib[n].BI_cell_1[3]].U[2][0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[2]].U[2][0]*cell_ib[n].dx2_1)*cell_ib[n].dy2_1 + 
                (cell[cell_ib[n].BI_cell_1[0]].U[2][0]*cell_ib[n].dx1_1 + cell[cell_ib[n].BI_cell_1[1]].U[2][0]*cell_ib[n].dx2_1)*cell_ib[n].dy1_1)/
            ((cell_ib[n].dy1_1 + cell_ib[n].dy2_1)*(cell_ib[n].dx1_1 + cell_ib[n].dx2_1));

        U2_2 =    ((cell[cell_ib[n].BI_cell_2[3]].U[2][0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[2]].U[2][0]*cell_ib[n].dx2_2)*cell_ib[n].dy2_2 + 
                (cell[cell_ib[n].BI_cell_2[0]].U[2][0]*cell_ib[n].dx1_2 + cell[cell_ib[n].BI_cell_2[1]].U[2][0]*cell_ib[n].dx2_2)*cell_ib[n].dy1_2)/
            ((cell_ib[n].dy1_2 + cell_ib[n].dy2_2)*(cell_ib[n].dx1_2 + cell_ib[n].dx2_2));

        d1 = cell_ib[n].d1; 
        d2 = cell_ib[n].d2;

        cell_ib[n].wall_flux_diffusive[0] = 0; 
        cell_ib[n].wall_flux_diffusive[1] = myu/rho*(d2/d1*(0 - U1_1) - d1/d2*(0 - U1_2))/(d2-d1)*cell_ib[n].d_cut;                  // For non-constant grid spacing    
        cell_ib[n].wall_flux_diffusive[2] = myu/rho*(d2/d1*(0 - U2_1) - d1/d2*(0 - U2_2))/(d2-d1)*cell_ib[n].d_cut;                  // For non-constant grid spacing    
        //        cell_ib[n].wall_flux_diffusive[1] = myu/rho*(0 - U1_1)/d1*cell_ib[n].d_cut;                  // For non-constant grid spacing    
        //        cell_ib[n].wall_flux_diffusive[2] = myu/rho*(0 - U2_1)/d1*cell_ib[n].d_cut;                  // For non-constant grid spacing    

        //printf("\nwall_flux1 = %lf, wall_flux2 = %lf",cell_ib[n].wall_flux_diffusive[1], cell_ib[n].wall_flux_diffusive[2]);
    }
}

