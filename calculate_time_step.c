#include "ibm.h"
#include "functions.h"

void calculate_time_step(int eq_type)
{   
    switch (eq_type){
        case 1:
            calculate_time_step_potential_flow();
            break;
        case 2:
            calculate_time_step_navier_stokes();
            break;
        default:
            printf("Error, bad eq_type input, quitting");
            break;
    }
}

void calculate_sound_speed()
{
    int n;

    for(n=1; n<=Total_cells; n++){
        cell[n].c = pow(pow(cell[n].U[1][0],2) + pow(cell[n].U[2][0],2) + acf, 0.5);
    }

    for(n=1; n<=N_ib; n++){
        cell_ib[n].c = pow(pow(cell_ib[n].U[1][0],2) + pow(cell_ib[n].U[2][0],2) + acf, 0.5);
    }

}

void calculate_time_step_potential_flow()
{
    double dt_min = 1000;
    double dt_local;
    int n;

    calculate_sound_speed();

    for(n=1; n<=Total_cells; n++){
        if(cell[n].vol_fr > 0.99){
            dt_local = diff_no*cell[n].vol_fr*(dx*dy)/acf;
            dt_min = min(dt_min, dt_local);        
        }
    }

    dt = dt_min;
    
//   printf("dt = %lf\n", dt);
}



// Still to add
void calculate_time_step_navier_stokes()
{

}
