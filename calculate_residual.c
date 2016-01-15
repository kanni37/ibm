#include "functions.h"
#include "ibm.h"

double calculate_residual(int eq_type, FILE* fp_residual, int count_iter)
{   
    double residual;
    double residual_u, residual_v, residual_p;

    switch (eq_type){
        case 1:
            residual = calculate_residual_potential_flow();
            break;
        case 2:
            residual_u = calculate_u_residual_navier_stokes();
            residual_v = calculate_v_residual_navier_stokes();
            residual_p = calculate_p_residual_navier_stokes();
            fprintf(fp_residual, "%d\t%lf\t%lf\t%lf\n", count_iter, log10(residual_p), log10(residual_u), log10(residual_v));  
            residual = max(max(residual_u, residual_v), max(residual_v, residual_p));
            break;
        default:
            residual = -1;
            printf("Error, bad eq_type input, quitting");
            break;
    }

    return residual;
}

double calculate_residual_potential_flow()
{
    int n;
    double residual = 0;        // residual after an iteration
    int cell_count = 0;         // cells participating in residual calculation

    #pragma omp parallel for reduction(+: residual, cell_count)
    for(n=1; n<=Total_cells; n++){
        if(cell[n].state != 1){
            if(fabs(cell[n].Psi[1]) > 0.0000001){
                residual += pow((cell[n].Psi[1] - cell[n].Psi[0])/(cell[n].Psi[1]),2);
                cell_count++;
            }
        }
    }
    if(cell_count == 0){
        residual = 1;
    }else{
        residual = pow(residual/cell_count,0.5);
    }

    return residual;
}


double calculate_u_residual_navier_stokes()
{
    int n;
    double residual_u = 0;        // residual after an iteration
    int cell_count_u = 0;         // cells participating in residual calculation

    #pragma omp parallel for reduction(+: residual_u, cell_count_u)
    for(n=1; n<=Total_cells; n++){
        if(cell[n].state == 0){
            if(fabs(cell[n].U[1][1]) > 0.0000001){
                residual_u += pow((cell[n].U[1][1] - cell[n].U[1][0])/(cell[n].U[1][1]),2);
                cell_count_u++;
            }
        }
    }

    if(cell_count_u == 0){
        residual_u = 1;
    }else{
        residual_u = pow(residual_u/cell_count_u,0.5);
    }

    return residual_u;
}


double calculate_v_residual_navier_stokes()
{
    int n;
    double residual_v = 0;        // residual after an iteration
    int cell_count_v = 0;         // cells participating in residual calculation

    #pragma omp parallel for reduction(+: residual_v, cell_count_v)
    for(n=1; n<=Total_cells; n++){
        if(cell[n].state == 0){
            if(fabs(cell[n].U[2][1]) > 0.0000001){
                residual_v += pow((cell[n].U[2][1] - cell[n].U[2][0])/(cell[n].U[2][1]),2);
                cell_count_v++;
            }
        }
    }

    if(cell_count_v == 0){
        residual_v = 1;
    }else{
        residual_v = pow(residual_v/cell_count_v,0.5);
    }

    return residual_v;
}

double calculate_p_residual_navier_stokes()
{
    int n;
    double residual_p = 0;        // residual after an iteration
    int cell_count_p = 0;         // cells participating in residual calculation

    #pragma omp parallel for reduction(+: residual_p, cell_count_p)
    for(n=1; n<=Total_cells; n++){
        if(cell[n].state == 0){
            if(fabs(cell[n].U[0][1]) > 0.0000001){
                residual_p += pow((cell[n].U[0][1] - cell[n].U[0][0])/(cell[n].U[0][1]), 2);
                cell_count_p++;
            }
        }
    }

    if(cell_count_p == 0){
        residual_p = 1;
    }else{
        residual_p = pow(residual_p/cell_count_p,0.5);
    }

    return residual_p;
}
