#include "functions.h"
#include "ibm.h"

void initialize(int eq_type, int init_type)
{
    switch (eq_type){
        case 1:
            initialize_potential_flow(init_type);
            break;
        case 2:
            initialize_navier_stokes(init_type);
            break;
        default:
            printf("Error, bad eq_type input, quitting");
            break;
    }
}

void initialize_potential_flow(int init_type)
{
    int i;
    if (init_type == 1){
        //Initial Values for cells
        for(i=1; i<=Total_cells; i++){
            cell[i].Psi[0] = Psi0;
            cell[i].Psi[1] = Psi0;
        }

        //Initial values for immersed cells
        for(i = 1; i<= N_ib; i++){
            cell_ib[i].Psi[0] = Psi0;
            cell_ib[i].Psi[1] = Psi0;
        }
    }else if(init_type == 2){
        FILE *fp = NULL;
        char buffer[100];

        fp = fopen("Cell_Data.dat","r"); 
        // First 2 reads correspond to tecplot syntax
        fgets(buffer, 100, fp);
        fgets(buffer, 100, fp);
        // Read the data from here on
        for(i=1; i<=Total_cells; i++){
            fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d", 
                    &cell[i].x, &cell[i].y, &cell[i].Psi[0], 
                    &cell[i].u[0], &cell[i].v[0], &cell[i].p[0], &cell[i].state);
        }
        fclose(fp);

        for(i = 1; i<= N_ib; i++){
            cell_ib[i].Psi[0] = cell[cell_ib[i].cell_no].Psi[0];
        }
    }else{
        printf("Init_type you specified does not exist in code !\n");
        getchar();
    }
}

void initialize_navier_stokes(int init_type)
{
    int i;
    if(init_type == 1){
        //Initial Values for cells
        for(i=1; i<=Total_cells; i++){
            cell[i].U[0][0] = P_out;
            cell[i].U[1][0] = U_ini;
            cell[i].U[2][0] = V_ini;
            cell[i].U[0][1] = P_out;
            cell[i].U[1][1] = U_ini;
            cell[i].U[2][1] = V_ini;
        }

        //Initial values for immersed cells
        for(i = 1; i<= N_ib; i++){
            cell_ib[i].U[0][0] = P_out;
            cell_ib[i].U[1][0] = U_ini;
            cell_ib[i].U[2][0] = V_ini;
        }
    }else if(init_type == 2){
        FILE *fp = NULL;
        char buffer[100];
        
        fp = fopen("Cell_Data.dat","r"); 
        
        if(fp == NULL){
            printf("File does not exist. Change init_type in ibm.c");
            getchar();
        }

        // First 2 reads correspond to tecplot syntax
        fgets(buffer, 100, fp);
        fgets(buffer, 100, fp);
        // Read the data from here on
        for(i=1; i<=Total_cells; i++){
            fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%d", 
                    &cell[i].x, &cell[i].y, 
                    &cell[i].U[1][0], &cell[i].U[2][0], &cell[i].U[0][0], &cell[i].state);
        }
        fclose(fp);

        //Initial values for immersed cells
        for(i = 1; i<= N_ib; i++){
            cell_ib[i].U[0][0] = cell[cell_ib[i].cell_no].U[0][0];
            cell_ib[i].U[1][0] = cell[cell_ib[i].cell_no].U[1][0];
            cell_ib[i].U[2][0] =  cell[cell_ib[i].cell_no].U[2][0];
        }
    }else{
        printf("\n----------------------------------------------\n");
        printf("Init_type you specified does not exist in code !\n");
        printf("----------------------------------------------\n");
        getchar();
    }
}
