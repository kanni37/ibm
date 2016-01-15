#include "ibm.h"

void generate_cartesian_mesh()
{
    int i, j, n;
    double Cell_Re;
    //-------------------Defining Cartesian Background Grid------------
    dx = 1.0*L/(Nx-1);
    dy = 1.0*H/(Ny-1);
    
    dmin = dist_factor*pow((pow(dx,2) + pow(dy,2)),0.5);       //Perpendicular distance from body (used when calculating first order accurate gradients)
    
    r = dy/dx;
    
    Re = rho*U_ini*2*Radius/myu;
    
    Cell_Re = rho*U_ini*min(dx, dy)/myu;

    printf("\ndx = %f,dy = %f", dx, dy);
    printf("\ndmin = %f", dmin);
    printf("\n\nFlow Reynolds Number = %f", Re);
    printf("\n\nCell Reynolds Number = %f", Cell_Re);

    Total_Grid_Points = Nx*Ny;
    
    point = malloc((Total_Grid_Points+1)*sizeof(struct points));
    
    for(j=1; j<=Ny; j++){
        for(i=1; i<=Nx; i++){
            n = (j-1)*Nx + i;
            point[n].x = (i-1)*dx;
            point[n].y = (j-1)*dy;
            //printf("point[%d].x = %lf\n", n, point[n].x);
        }
    }
}

