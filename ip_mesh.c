#include "ibm.h"
#include "functions.h"

void ip_mesh(int M_type)
{
    switch (M_type){
        case 1:
            transfinite_ip_mesh();
            break;
        case 2:
            johansen_colella_mesh();
            break;
        case 3:
            johansen_colella_mesh_modified();
            break;
        default:
            printf("Error, bad input, quitting");
            break;
    }
}

void transfinite_ip_mesh()
{
    int i, j, n;
    double dxp, dyp;
    double m_x, m_y;
    int n_ip, I_ip = 0, J_ip = 0;    


    // Finding transfinite grid points to be used for enforcing boundary condition
    // 2 layers used in current case
    for (n=1; n<=N_ib; n++) {

        cell_ib[n].x_m = (cell_ib[n].x_cut[0] + cell_ib[n].x_cut[1])/2;
        cell_ib[n].y_m = (cell_ib[n].y_cut[0] + cell_ib[n].y_cut[1])/2;
        //printf("\n xm,ym = %lf,%lf ", cell_ib[n].x_m,  cell_ib[n].y_m);

        if(cell_ib[n].d_cut == 0){
            m_x = 0;
            m_y = 0;
        }else{
            m_x =  cell_ib[n].slope_x_cut[0];
            m_y =  cell_ib[n].slope_y_cut[0];
        }

        dxp = -m_y*dmin;
        dyp = m_x*dmin;

        cell_ib[n].x_1 = cell_ib[n].x_cut[0]  + dxp;
        cell_ib[n].y_1 = cell_ib[n].y_cut[0]  + dyp;
        cell_ib[n].x_3 = cell_ib[n].x_cut[0]  + 2*dxp;
        cell_ib[n].y_3 = cell_ib[n].y_cut[0]  + 2*dyp;

        if(is_point_in_polygon(cell_ib[n].x_1, cell_ib[n].y_1) == 1){
            cell_ib[n].x_1 = cell_ib[n].x_cut[0] - dxp;
            cell_ib[n].y_1 = cell_ib[n].y_cut[0] - dyp;
            cell_ib[n].x_3 = cell_ib[n].x_cut[0] - 2*dxp;
            cell_ib[n].y_3 = cell_ib[n].y_cut[0] - 2*dyp;
        }

        if(cell_ib[n].d_cut == 0){
            m_x = 0;
            m_y = 0;
        }else{
            m_x =  cell_ib[n].slope_x_cut[1];
            m_y =  cell_ib[n].slope_y_cut[1];
        }

        dxp = -m_y*dmin;
        dyp = m_x*dmin;

        cell_ib[n].x_2 = cell_ib[n].x_cut[1]  + dxp;
        cell_ib[n].y_2 = cell_ib[n].y_cut[1]  + dyp;
        cell_ib[n].x_4 = cell_ib[n].x_cut[1]  + 2*dxp;
        cell_ib[n].y_4 = cell_ib[n].y_cut[1]  + 2*dyp;

        if(is_point_in_polygon(cell_ib[n].x_2, cell_ib[n].y_2) == 1){
            cell_ib[n].x_2 = cell_ib[n].x_cut[1] - dxp;
            cell_ib[n].y_2 = cell_ib[n].y_cut[1] - dyp;
            cell_ib[n].x_4 = cell_ib[n].x_cut[1] - 2*dxp;
            cell_ib[n].y_4 = cell_ib[n].y_cut[1] - 2*dyp;
        }
    }

    /*
       for(n=1; n<=N_ib; n++){
       printf("\nIb_cell %d , x = %lf, y = %lf, x_p = %lf, y_p = %lf", 
       n, cell_ib[n].x_cut[0], cell_ib[n].y_cut[0], cell_ib[n].xT[0], cell_ib[n].yT[0]);
       }
       */

    FILE *fp = NULL;
    fp = fopen("Transfinite_Coordinates.dat", "w");
    fprintf(fp,"VARIABLES = \"x\", \"y\"\n");
    for(n=1; n<=N_ib; n++){
        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_cut[0], cell_ib[n].y_cut[0]);
        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_1, cell_ib[n].y_1);
        //        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_2, cell_ib[n].y_2);
        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_3, cell_ib[n].y_3);
        //        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_4, cell_ib[n].y_4);
    }
    fclose(fp);


    // cell_ib[n].BI_cell: Finding Cells from which to interpolate the solution on transfinite grid
    // cell_ib[n].dx1,dx2,dy1,dy2: Finding information required to do bilinear interpolation to find ...
    // data on transfinite grid points
    for (n=1; n<=N_ib; n++) {
        for (i=1; i<=cell_Nx; i++){
            n_ip = i;
            if (cell_ib[n].x_1<=cell[n_ip].x){
                I_ip = i;
                break;
            }
        }
        for (j=1; j<=cell_Ny; j++){
            n_ip = (j-1)*cell_Nx + 1;
            if (cell_ib[n].y_1<=cell[n_ip].y){
                J_ip = j;
                break;
            }
        }

        //printf("\ncellIB = %d I_ip = %d, J_ip = %d", n, I_ip, J_ip);

        cell_ib[n].BI_cell_1[0] = gcell(I_ip-1,J_ip-1);
        cell_ib[n].BI_cell_1[1] = gcell(I_ip,J_ip-1);
        cell_ib[n].BI_cell_1[2] = gcell(I_ip,J_ip);
        cell_ib[n].BI_cell_1[3] = gcell(I_ip-1,J_ip);

        cell_ib[n].dx1_1 = fabs(cell_ib[n].x_1 - cell[cell_ib[n].BI_cell_1[2]].x);
        cell_ib[n].dx2_1 = fabs(cell_ib[n].x_1 - cell[cell_ib[n].BI_cell_1[3]].x);
        cell_ib[n].dy1_1 = fabs(cell_ib[n].y_1 - cell[cell_ib[n].BI_cell_1[2]].y);
        cell_ib[n].dy2_1 = fabs(cell_ib[n].y_1 - cell[cell_ib[n].BI_cell_1[1]].y);

        // printf("\nBI_cellx for 0,1,2,3 = %lf, %lf, %lf, %lf",
        //         cell[cell_ib[n].BI_cell[k][0]].x, cell[cell_ib[n].BI_cell[k][1]].x,
        //         cell[cell_ib[n].BI_cell[k][2]].x, cell[cell_ib[n].BI_cell[k][3]].x);
        // printf("\nBI_cellx for 0,1,2,3 = %lf, %lf, %lf, %lf",
        //           cell_ib[n].dxT1[k], cell_ib[n].dxT2[k],
        //           cell_ib[n].dyT1[k], cell_ib[n].dyT2[k]);
    }

    for (n=1; n<=N_ib; n++) {
        for (i=1; i<=cell_Nx; i++){
            n_ip = i;
            if (cell_ib[n].x_2<=cell[n_ip].x){
                I_ip = i;
                break;
            }
        }

        for (j=1; j<=cell_Ny; j++){
            n_ip = (j-1)*cell_Nx + 1;
            if (cell_ib[n].y_2<=cell[n_ip].y){
                J_ip = j;
                break;
            }
        }
        //printf("\ncellIB = %d I_ip = %d, J_ip = %d", n, I_ip, J_ip);

        cell_ib[n].BI_cell_2[0] = gcell(I_ip-1,J_ip-1);
        cell_ib[n].BI_cell_2[1] = gcell(I_ip,J_ip-1);
        cell_ib[n].BI_cell_2[2] = gcell(I_ip,J_ip);
        cell_ib[n].BI_cell_2[3] = gcell(I_ip-1,J_ip);

        cell_ib[n].dx1_2 = fabs(cell_ib[n].x_2 - cell[cell_ib[n].BI_cell_2[2]].x);
        cell_ib[n].dx2_2 = fabs(cell_ib[n].x_2 - cell[cell_ib[n].BI_cell_2[3]].x);
        cell_ib[n].dy1_2 = fabs(cell_ib[n].y_2 - cell[cell_ib[n].BI_cell_2[2]].y);
        cell_ib[n].dy2_2 = fabs(cell_ib[n].y_2 - cell[cell_ib[n].BI_cell_2[1]].y);
    }

    for (n=1; n<=N_ib; n++) {
        for (i=1; i<=cell_Nx; i++){
            n_ip = i;
            if (cell_ib[n].x_3<=cell[n_ip].x){
                I_ip = i;
                break;
            }
        }
        for (j=1; j<=cell_Ny; j++){
            n_ip = (j-1)*cell_Nx + 1;
            if (cell_ib[n].y_3<=cell[n_ip].y){
                J_ip = j;
                break;
            }
        }
        //printf("\ncellIB = %d I_ip = %d, J_ip = %d", n, I_ip, J_ip);

        cell_ib[n].BI_cell_3[0] = gcell(I_ip-1,J_ip-1);
        cell_ib[n].BI_cell_3[1] = gcell(I_ip,J_ip-1);
        cell_ib[n].BI_cell_3[2] = gcell(I_ip,J_ip);
        cell_ib[n].BI_cell_3[3] = gcell(I_ip-1,J_ip);

        cell_ib[n].dx1_3 = fabs(cell_ib[n].x_3 - cell[cell_ib[n].BI_cell_3[2]].x);
        cell_ib[n].dx2_3 = fabs(cell_ib[n].x_3 - cell[cell_ib[n].BI_cell_3[3]].x);
        cell_ib[n].dy1_3 = fabs(cell_ib[n].y_3 - cell[cell_ib[n].BI_cell_3[2]].y);
        cell_ib[n].dy2_3 = fabs(cell_ib[n].y_3 - cell[cell_ib[n].BI_cell_3[1]].y);

        // printf("\nBI_cellx for 0,1,2,3 = %lf, %lf, %lf, %lf",
        //            cell[cell_ib[n].BI_cell[k][0]].x, cell[cell_ib[n].BI_cell[k][1]].x,
        //            cell[cell_ib[n].BI_cell[k][2]].x, cell[cell_ib[n].BI_cell[k][3]].x);
        // printf("\nBI_cellx for 0,1,2,3 = %lf, %lf, %lf, %lf",
        //           cell_ib[n].dxT1[k], cell_ib[n].dxT2[k],
        //           cell_ib[n].dyT1[k], cell_ib[n].dyT2[k]);
    }

    for (n=1; n<=N_ib; n++) {
        for (i=1; i<=cell_Nx; i++){
            n_ip = i;
            if (cell_ib[n].x_4<=cell[n_ip].x){
                I_ip = i;
                break;
            }
        }
        for (j=1; j<=cell_Ny; j++){
            n_ip = (j-1)*cell_Nx + 1;
            if (cell_ib[n].y_4<=cell[n_ip].y){
                J_ip = j;
                break;
            }
        }
        //printf("\ncellIB = %d I_ip = %d, J_ip = %d", n, I_ip, J_ip);

        cell_ib[n].BI_cell_4[0] = gcell(I_ip-1,J_ip-1);
        cell_ib[n].BI_cell_4[1] = gcell(I_ip,J_ip-1);
        cell_ib[n].BI_cell_4[2] = gcell(I_ip,J_ip);
        cell_ib[n].BI_cell_4[3] = gcell(I_ip-1,J_ip);

        cell_ib[n].dx1_4 = fabs(cell_ib[n].x_4 - cell[cell_ib[n].BI_cell_4[2]].x);
        cell_ib[n].dx2_4 = fabs(cell_ib[n].x_4 - cell[cell_ib[n].BI_cell_4[3]].x);
        cell_ib[n].dy1_4 = fabs(cell_ib[n].y_4 - cell[cell_ib[n].BI_cell_4[2]].y);
        cell_ib[n].dy2_4 = fabs(cell_ib[n].y_4 - cell[cell_ib[n].BI_cell_4[1]].y);

        // printf("\nBI_cellx for 0,1,2,3 = %lf, %lf, %lf, %lf",
        //             cell[cell_ib[n].BI_cell[k][0]].x, cell[cell_ib[n].BI_cell[k][1]].x,
        //             cell[cell_ib[n].BI_cell[k][2]].x, cell[cell_ib[n].BI_cell[k][3]].x);
        // printf("\nBI_cellx for 0,1,2,3 = %lf, %lf, %lf, %lf",
        //           cell_ib[n].dxT1[k], cell_ib[n].dxT2[k], cell_ib[n].dyT1[k], cell_ib[n].dyT2[k]);
    }

    int k;
    double zi_1, zi_2, zi_m;
    double x_zi, y_zi, x_eta, y_eta;
    double Jacob;

    for (n=1; n<=N_ib; n++) {
        for(k=0; k<cell_ib[n].N_local_division; k++){
            zi_1 = 1.0*k/cell_ib[n].N_local_division;
            zi_2 = zi_1 + 1.0/cell_ib[n].N_local_division;
            zi_m = (zi_1 + zi_2)/2;

            cell_ib[n].del_x[k] = Radius*(cos(cell_ib[n].theta_cut[0] + zi_2*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])) 
                    - cos(cell_ib[n].theta_cut[0] + zi_1*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])));
            cell_ib[n].del_y[k] = Radius*(sin(cell_ib[n].theta_cut[0] + zi_2*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])) 
                    - sin(cell_ib[n].theta_cut[0] + zi_1*(cell_ib[n].theta_cut[1] - cell_ib[n].theta_cut[0])));
            //    printf("del_x = %lf, del_y = %lf\n",del_x, del_y);

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
            cell_ib[n].eta_x[k] = -y_zi/Jacob;
            cell_ib[n].eta_y[k] = x_zi/Jacob;
        }
    }
}

void johansen_colella_mesh()
{
    int i, j, n;
    double dxp, dyp;
    double m_x, m_y;
    double n_x, n_y;
    int n_ip, I_ip = 0, J_ip = 0;    
    double t1, t2;
    int N;

    // Finding normals
    for (n=1; n<=N_ib; n++) {

        cell_ib[n].x_1 = xo;
        cell_ib[n].x_2 = xo;
        cell_ib[n].y_1 = yo;
        cell_ib[n].y_2 = yo;
        
        cell_ib[n].x_m = (cell_ib[n].x_cut[0] + cell_ib[n].x_cut[1])/2;
        cell_ib[n].y_m = (cell_ib[n].y_cut[0] + cell_ib[n].y_cut[1])/2;

        if(cell_ib[n].d_cut == 0){
            cell_ib[n].d1 = 0; 
            cell_ib[n].d2 = 0; 
        }
        else{
            m_x = (cell_ib[n].x_cut[1] - cell_ib[n].x_cut[0])/cell_ib[n].d_cut;
            m_y = (cell_ib[n].y_cut[1] - cell_ib[n].y_cut[0])/cell_ib[n].d_cut;
            n_x = -m_y;
            n_y = m_x;

            dxp = dmin*n_x;
            dyp = dmin*n_y;

            cell_ib[n].x_1 = cell_ib[n].x_m  + dxp;
            cell_ib[n].y_1 = cell_ib[n].y_m  + dyp;

            if(is_point_in_polygon(cell_ib[n].x_1, cell_ib[n].y_1) == 1){
                n_x = -n_x;
                n_y = -n_y;
            }
            cell_ib[n].n_x = n_x;
            cell_ib[n].n_y = n_y;

            // Finding the points in flow domain for boundary flux calculations ...
            // as used by Johansen and Colella
            N = cell_ib[n].cell_no;

            if (fabs(n_x) > fabs(n_y)){
                if (n_x >= 0){
                    cell_ib[n].x_1 = cell[N+1].x;
                    cell_ib[n].x_2 = cell[N+2].x;

                    t1 = (cell_ib[n].x_1 - cell_ib[n].x_m)/n_x;
                    t2 = (cell_ib[n].x_2 - cell_ib[n].x_m)/n_x;

                    cell_ib[n].y_1 = cell_ib[n].y_m + n_y*t1;
                    cell_ib[n].y_2 = cell_ib[n].y_m + n_y*t2;
                }else{
                    cell_ib[n].x_1 = cell[N-1].x;
                    cell_ib[n].x_2 = cell[N-2].x;

                    t1 = (cell_ib[n].x_1 - cell_ib[n].x_m)/n_x;
                    t2 = (cell_ib[n].x_2 - cell_ib[n].x_m)/n_x;

                    cell_ib[n].y_1 = cell_ib[n].y_m + n_y*t1;
                    cell_ib[n].y_2 = cell_ib[n].y_m + n_y*t2;
                }
            }else{ 
                if (n_y >= 0){
                    cell_ib[n].y_1 = cell[N+cell_Nx].y;
                    cell_ib[n].y_2 = cell[N+2*cell_Nx].y;

                    t1 = (cell_ib[n].y_1 - cell_ib[n].y_m)/n_y;
                    t2 = (cell_ib[n].y_2 - cell_ib[n].y_m)/n_y;

                    cell_ib[n].x_1 = cell_ib[n].x_m + n_x*t1;
                    cell_ib[n].x_2 = cell_ib[n].x_m + n_x*t2;
                }else{
                    cell_ib[n].y_1 = cell[N-cell_Nx].y;
                    cell_ib[n].y_2 = cell[N-2*cell_Nx].y;

                    t1 = (cell_ib[n].y_1 - cell_ib[n].y_m)/n_y;
                    t2 = (cell_ib[n].y_2 - cell_ib[n].y_m)/n_y;

                    cell_ib[n].x_1 = cell_ib[n].x_m + n_x*t1;
                    cell_ib[n].x_2 = cell_ib[n].x_m + n_x*t2;
                }
            }

            cell_ib[n].d1 = pow((pow(cell_ib[n].x_1 - cell_ib[n].x_m,2) + pow(cell_ib[n].y_1 - cell_ib[n].y_m,2)),0.5);
            cell_ib[n].d2 = pow((pow(cell_ib[n].x_2 - cell_ib[n].x_m,2) + pow(cell_ib[n].y_2 - cell_ib[n].y_m,2)),0.5);
        }   
    }

    FILE *fp = NULL;
    fp = fopen("Johansen_Colella_Mesh.dat", "w");
    fprintf(fp,"VARIABLES = \"x\", \"y\"\n");
    for(n=1; n<=N_ib; n++){
        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_m, cell_ib[n].y_m);
        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_1, cell_ib[n].y_1);
        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_2, cell_ib[n].y_2);
    }
    fclose(fp);

    for(n=1; n<=N_ib; n++){
        printf("%lf\t%lf\n", cell_ib[n].d1, cell_ib[n].d2);
    }

    // Finding surrounding cells to use for bilinear interpolation for cell 
    for (n=1; n<=N_ib; n++) {
        for (i=1; i<=cell_Nx; i++){
            n_ip = i;
            if (cell_ib[n].x_1<=cell[n_ip].x){
                I_ip = i;
                break;
            }
        }
        for (j=1; j<=cell_Ny; j++){
            n_ip = (j-1)*cell_Nx + 1;
            if (cell_ib[n].y_1<=cell[n_ip].y){
                J_ip = j;
                break;
            }
        }

        cell_ib[n].BI_cell_1[0] = gcell(I_ip-1,J_ip-1);
        cell_ib[n].BI_cell_1[1] = gcell(I_ip,J_ip-1);
        cell_ib[n].BI_cell_1[2] = gcell(I_ip,J_ip);
        cell_ib[n].BI_cell_1[3] = gcell(I_ip-1,J_ip);

        //Interpolation cell search ends
        cell_ib[n].dx1_1 = fabs(cell_ib[n].x_1 - cell[cell_ib[n].BI_cell_1[2]].x);
        cell_ib[n].dx2_1 = fabs(cell_ib[n].x_1 - cell[cell_ib[n].BI_cell_1[3]].x);
        cell_ib[n].dy1_1 = fabs(cell_ib[n].y_1 - cell[cell_ib[n].BI_cell_1[2]].y);
        cell_ib[n].dy2_1 = fabs(cell_ib[n].y_1 - cell[cell_ib[n].BI_cell_1[1]].y);
    }

    // Finding surrounding cells to use for bilinear interpolation for cell 
    for (n=1; n<=N_ib; n++) {
        for (i=1; i<=cell_Nx; i++){
            n_ip = i;
            if (cell_ib[n].x_2<=cell[n_ip].x){
                I_ip = i;
                break;
            }
        }
        for (j=1; j<=cell_Ny; j++){
            n_ip = (j-1)*cell_Nx + 1;
            if (cell_ib[n].y_2<=cell[n_ip].y){
                J_ip = j;
                break;
            }
        }

        cell_ib[n].BI_cell_2[0] = gcell(I_ip-1,J_ip-1);
        cell_ib[n].BI_cell_2[1] = gcell(I_ip,J_ip-1);
        cell_ib[n].BI_cell_2[2] = gcell(I_ip,J_ip);
        cell_ib[n].BI_cell_2[3] = gcell(I_ip-1,J_ip);
        //Interpolation cell search ends

        cell_ib[n].dx1_2 = fabs(cell_ib[n].x_2 - cell[cell_ib[n].BI_cell_2[2]].x);
        cell_ib[n].dx2_2 = fabs(cell_ib[n].x_2 - cell[cell_ib[n].BI_cell_2[3]].x);
        cell_ib[n].dy1_2 = fabs(cell_ib[n].y_2 - cell[cell_ib[n].BI_cell_2[2]].y);
        cell_ib[n].dy2_2 = fabs(cell_ib[n].y_2 - cell[cell_ib[n].BI_cell_2[1]].y);
    }
}

void johansen_colella_mesh_modified()
{
    int i, j, n;
    double dxp, dyp;
    double m_x, m_y;
    int n_ip, I_ip = 0, J_ip = 0;    

    for (n=1; n<=N_ib; n++) {

        cell_ib[n].d1 = dmin;
        cell_ib[n].d2 = 2.0*dmin;

        cell_ib[n].x_m = (cell_ib[n].x_cut[0] + cell_ib[n].x_cut[1])/2.0;
        cell_ib[n].y_m = (cell_ib[n].y_cut[0] + cell_ib[n].y_cut[1])/2.0;

        if(cell_ib[n].d_cut == 0){
            m_x = 0;
            m_y = 0;
        }
        else{
            m_x = (cell_ib[n].x_cut[1] - cell_ib[n].x_cut[0])/cell_ib[n].d_cut;
            m_y = (cell_ib[n].y_cut[1] - cell_ib[n].y_cut[0])/cell_ib[n].d_cut;
        }

        cell_ib[n].n_x = -m_y;
        cell_ib[n].n_y = m_x;

        dxp = -m_y*dmin;
        dyp = m_x*dmin;

        cell_ib[n].x_1 = cell_ib[n].x_m  + 1.0*dxp;
        cell_ib[n].y_1 = cell_ib[n].y_m  + 1.0*dyp;

        cell_ib[n].x_2 = cell_ib[n].x_m  + 2.0*dxp;
        cell_ib[n].y_2 = cell_ib[n].y_m  + 2.0*dyp;

        if(is_point_in_polygon(cell_ib[n].x_1, cell_ib[n].y_1) == 1){
            cell_ib[n].x_1 = cell_ib[n].x_m  - 1.0*dxp;
            cell_ib[n].y_1 = cell_ib[n].y_m  - 1.0*dyp;
            cell_ib[n].x_2 = cell_ib[n].x_m  - 2.0*dxp;
            cell_ib[n].y_2 = cell_ib[n].y_m  - 2.0*dyp;
            cell_ib[n].n_x = -cell_ib[n].n_x;
            cell_ib[n].n_y = -cell_ib[n].n_y;
        }
    }

    FILE *fp = NULL;
    fp = fopen("Johansen_Colella_Mesh.dat", "w");
    fprintf(fp,"VARIABLES = \"x\", \"y\"\n");
    for(n=1; n<=N_ib; n++){
        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_m, cell_ib[n].y_m);
        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_1, cell_ib[n].y_1);
        fprintf(fp,"%lf\t%lf\n", cell_ib[n].x_2, cell_ib[n].y_2);
    }
    fclose(fp);

    for (n=1; n<=N_ib; n++) {
        for (i=1; i<=cell_Nx; i++){
            n_ip = i;
            if (cell_ib[n].x_1<=cell[n_ip].x){
                I_ip = i;
                break;
            }
        }
        for (j=1; j<=cell_Ny; j++){
            n_ip = (j-1)*cell_Nx + 1;
            if (cell_ib[n].y_1<=cell[n_ip].y){
                J_ip = j;
                break;
            }
        }

        cell_ib[n].BI_cell_1[0] = gcell(I_ip-1,J_ip-1);
        cell_ib[n].BI_cell_1[1] = gcell(I_ip,J_ip-1);
        cell_ib[n].BI_cell_1[2] = gcell(I_ip,J_ip);
        cell_ib[n].BI_cell_1[3] = gcell(I_ip-1,J_ip);

        //Interpolation cell search ends
        cell_ib[n].dx1_1 = fabs(cell_ib[n].x_1 - cell[cell_ib[n].BI_cell_1[2]].x);
        cell_ib[n].dx2_1 = fabs(cell_ib[n].x_1 - cell[cell_ib[n].BI_cell_1[3]].x);
        cell_ib[n].dy1_1 = fabs(cell_ib[n].y_1 - cell[cell_ib[n].BI_cell_1[2]].y);
        cell_ib[n].dy2_1 = fabs(cell_ib[n].y_1 - cell[cell_ib[n].BI_cell_1[1]].y);
    }

    for (n=1; n<=N_ib; n++) {
        for (i=1; i<=cell_Nx; i++){
            n_ip = i;
            if (cell_ib[n].x_2<=cell[n_ip].x){
                I_ip = i;
                break;
            }
        }
        for (j=1; j<=cell_Ny; j++){
            n_ip = (j-1)*cell_Nx + 1;
            if (cell_ib[n].y_2<=cell[n_ip].y){
                J_ip = j;
                break;
            }
        }

        cell_ib[n].BI_cell_2[0] = gcell(I_ip-1,J_ip-1);
        cell_ib[n].BI_cell_2[1] = gcell(I_ip,J_ip-1);
        cell_ib[n].BI_cell_2[2] = gcell(I_ip,J_ip);
        cell_ib[n].BI_cell_2[3] = gcell(I_ip-1,J_ip);
        //Interpolation cell search ends

        cell_ib[n].dx1_2 = fabs(cell_ib[n].x_2 - cell[cell_ib[n].BI_cell_2[2]].x);
        cell_ib[n].dx2_2 = fabs(cell_ib[n].x_2 - cell[cell_ib[n].BI_cell_2[3]].x);
        cell_ib[n].dy1_2 = fabs(cell_ib[n].y_2 - cell[cell_ib[n].BI_cell_2[2]].y);
        cell_ib[n].dy2_2 = fabs(cell_ib[n].y_2 - cell[cell_ib[n].BI_cell_2[1]].y);
    }
}
