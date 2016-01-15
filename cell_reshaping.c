#include "ibm.h"
#include "functions.h"

void cell_reshaping()
{
    double m_x, m_y;
    int n;

    // Finding cell centers inside and outside the boundary
    for(n=1;n<=N_ib;n++){
        if(is_point_in_polygon(cell_ib[n].x, cell_ib[n].y) == 1){
            cell_ib[n].position = 1;           // Solid
        }else{
            cell_ib[n].position = 0;           // Fluid
        }
    }

    // Finding Parent and child cells
    // Parent   :   Larger cell with which the smaller cell is merged
    // Child    :   Smaller cell which is merged with larger/pared cell
    int k;

    for(n=1; n<=N_ib; n++){
        if(cell_ib[n].position == 1){
            k = cell_ib[n].cell_no;
            cell[k].child = NULL;                     
            if(fabs(cell_ib[n].n_x) > fabs(cell_ib[n].n_y)){
                if(cell_ib[n].n_x >= 0){
                    cell[k].parent = k + 1;
                    cell[k+1].child = k;
                    cell[k+1].parent = NULL;
                }else{
                    cell[k].parent = k - 1;
                    cell[k-1].child = k;
                    cell[k-1].parent = NULL;
                }    
            }else{
                if(cell_ib[n].n_y >= 0){
                    cell[k].parent = k + cell_Nx;
                    cell[k+cell_Nx].child = k;
                    cell[k+cell_Nx].parent = NULL;
                }else{
                    cell[k].parent = k - cell_Nx;
                    cell[k-cell_Nx].child = k;
                    cell[k-cell_Nx].parent = NULL;
                }    
            }
        }
    } 
}
