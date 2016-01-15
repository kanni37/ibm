#include "functions.h"
#include "ibm.h"

void ibm_mesh(int D_type)
{
    switch (D_type){
        case 1:
            ibm_FV_mesh();
            break;
        case 2:
            ibm_FD_mesh();
            break;
        default:
            printf("Error, bad input, quitting");
            break;
    }
}

void ibm_FD_mesh()
{
    int n;

    for(n=1;n<=Total_Grid_Points;n++){
        if(is_point_in_polygon(point[n].x, point[n].y) == 1){
            point[n].state = 1;
        }else{
            point[n].state = 0;
        }
    }

    // Identifying Point states
    N_fl = 0;           // Number of fluid points
    N_ib = 0;           // Number of immersed points
    N_sl = 0;           // Number of Solid points

    for(n=1; n<=Total_Grid_Points; n++){
        if (point[n].state == 1) {
            if( point[n+1].state == 0  || point[n-1].state == 0  ||
                    point[n+Nx].state == 0 || point[n-Nx].state == 0 ) {
                point[n].state = 2;
                N_ib++;
            }else{
                N_sl++;
            }
        }else{
            N_fl++;
        }
    }
}

void ibm_FV_mesh()
{
    int i, j, k, n;
    int count, intersection_count, count_face, cut_face_count;
    double m, t;
    double temp; 
    double Flux_center_x, Flux_center_y;
    double b[2];

    for(n=1;n<=Total_Grid_Points;n++){
        if(is_point_in_polygon(point[n].x, point[n].y) == 1){
            point[n].state = 1;
        }else{
            point[n].state = 0;
        }
    }

    // Creating cells and required information
    cell_Nx = Nx-1;
    cell_Ny = Ny-1;

    Total_cells = cell_Nx*cell_Ny;

    cell = malloc((Total_cells+1)*sizeof(struct cells));

    // Calculations for cell vertices
    for(j=1; j<=cell_Ny; j++){
        for(i=1; i<=cell_Nx; i++){
            n = (j-1)*cell_Nx + i;
            k = (j-1)*Nx + i;
            cell[n].vertex[0] = k;
            cell[n].vertex[1] = k+1;
            cell[n].vertex[2] = k+1+Nx;
            cell[n].vertex[3] = k+Nx;
            //printf("point[%d].x = %lf\n", n, point[n].x);
        }
    }
    // Calculations for cell centers
    for(n=1; n<=Total_cells; n++){
        cell[n].x = (point[cell[n].vertex[0]].x + point[cell[n].vertex[1]].x + point[cell[n].vertex[2]].x + point[cell[n].vertex[3]].x)/4;
        cell[n].y = (point[cell[n].vertex[0]].y + point[cell[n].vertex[1]].y + point[cell[n].vertex[2]].y + point[cell[n].vertex[3]].y)/4;
    }

    // Calculating dimensions of cells
    for(j=1; j<=cell_Ny; j++){
        for(i=1; i<=cell_Nx; i++){
            n = (j-1)*cell_Nx + i;
            cell[n].dx = fabs(point[cell[n].vertex[0]].x - point[cell[n].vertex[1]].x);
            cell[n].dy = fabs(point[cell[n].vertex[2]].y - point[cell[n].vertex[3]].y);
        }
    }
    
    // Defining faces of cells
    for(n=1; n<=Total_cells; n++){
        cell[n].face[0][0] = cell[n].vertex[0];
        cell[n].face[0][1] = cell[n].vertex[1];
        cell[n].face[1][0] = cell[n].vertex[1];
        cell[n].face[1][1] = cell[n].vertex[2];
        cell[n].face[2][0] = cell[n].vertex[3];
        cell[n].face[2][1] = cell[n].vertex[2];
        cell[n].face[3][0] = cell[n].vertex[0];
        cell[n].face[3][1] = cell[n].vertex[3];
    }

    // Identifying Cell states
    N_fl = 0;                   // Number of fluid cells
    N_ib = 0;                   // Number of immersed cells
    N_sl = 0;                   // Number of Solid Cells

    for(n=1;n<=Total_cells;n++){
        cell[n].state = 2;
        count = 0;

        for(k=0; k<4; k++){
            if(point[cell[n].vertex[k]].state == 1){
                count+= 1;
            }else if(point[cell[n].vertex[k]].state == 0){
                count+=-1;
            }
        }

        if(count==4){
            cell[n].state = 1;          //Solid
            N_sl++;
        }else if(count==-4){
            cell[n].state = 0;          //Fluid
            N_fl++;
        }
    }

    N_ib = Total_cells - N_sl - N_fl;

    printf("\nNo: of Fluid Cells = %d", N_fl);
    printf("\nNo: of Solid Cells = %d", N_sl);
    printf("\nNo: of Immersed Cells = %d\n", N_ib);

    cell_ib = malloc((N_ib+1)*sizeof(struct ib_cells));

    count = 1;
    for(n=1; n<=Total_cells; n++){
        if(cell[n].state == 2){
            cell_ib[count].cell_no = n;
            cell[n].ib_cell_no = count;

            cell_ib[count].x = cell[n].x;
            cell_ib[count].y = cell[n].y;

            for(k=0; k<=3; k++){
                cell_ib[count].vertex[k] = cell[n].vertex[k];
            }

            for(k=0; k<=3; k++){
                cell_ib[count].face[k][0] = cell[n].face[k][0];
                cell_ib[count].face[k][1] = cell[n].face[k][1];
            }
            count++;
        }
    }
    
    for(n=1;n<=N_ib;n++){
        if(is_point_in_polygon(cell_ib[n].x, cell_ib[n].y) == 1){
            cell_ib[n].position = 1;           // Solid
        }else{
            cell_ib[n].position = 0;           // Fluid
        }
    }
    
    /*
       for(n=1; n<=N_ib; n++){
       printf("\nIb_cell %d = %d, x = %lf, y = %lf", n, cell_ib[n].cell_no, cell_ib[n].x, cell_ib[n].y);
       }
       */

    //----------------------Finding ratios of length cut by the body on background cartesian grid (Finding alpha values)-----------------------

    for(n=1; n<=N_ib; n++){
        intersection_count = 0;
        count_face = 0;

        for(k=0; k<4; k++){
            cut_face_count = 0;

            if(point[cell_ib[n].face[k][0]].state == 0 && point[cell_ib[n].face[k][1]].state == 0){
                cell_ib[n].alpha[k] = 1;
            }else if(point[cell_ib[n].face[k][0]].state == 1 && point[cell_ib[n].face[k][1]].state == 1){
                cell_ib[n].alpha[k] = 0;
            }else{
                cell_ib[n].cut_face[count_face] = k;
                if(k==0 || k==2){
                    for(i=1; i<=Total_Body_Points;i++){
                        if(cut_face_count == 0){
                            if(max(point_Body[face_Body[i].node[0]].y,point_Body[face_Body[i].node[1]].y)<point[cell_ib[n].face[k][0]].y||
                                    min(point_Body[face_Body[i].node[0]].y,point_Body[face_Body[i].node[1]].y)>point[cell_ib[n].face[k][0]].y){
                                //Do Nothing   //Drop all the faces which do not intersect with the face
                            }else{
                                m = (point[cell_ib[n].face[k][0]].y-point_Body[face_Body[i].node[0]].y)
                                    /(point_Body[face_Body[i].node[1]].y-point_Body[face_Body[i].node[0]].y);
                                t = (point_Body[face_Body[i].node[0]].x - point[cell_ib[n].face[k][0]].x + 
                                        m*(point_Body[face_Body[i].node[1]].x-point_Body[face_Body[i].node[0]].x))
                                    /(point[cell_ib[n].face[k][1]].x - point[cell_ib[n].face[k][0]].x);
                                if((m>=0 && m<=1) && (t>=0 && t<=1)){
                                    if(point[cell_ib[n].face[k][0]].state == 1){
                                        cell_ib[n].alpha[k] = 1-t;
                                        cell_ib[n].x_cut[intersection_count] = point[cell_ib[n].face[k][0]].x + t*(point[cell_ib[n].face[k][1]].x - point[cell_ib[n].face[k][0]].x);
                                        cell_ib[n].y_cut[intersection_count] = point[cell_ib[n].face[k][0]].y + t*(point[cell_ib[n].face[k][1]].y - point[cell_ib[n].face[k][0]].y);
                                        cell_ib[n].slope_x_cut[intersection_count] = (point_Body[face_Body[i].node[1]].x - point_Body[face_Body[i].node[0]].x)/
                                            (pow(pow(point_Body[face_Body[i].node[1]].x - point_Body[face_Body[i].node[0]].x,2) + 
                                                 pow(point_Body[face_Body[i].node[1]].y - point_Body[face_Body[i].node[0]].y,2),0.5));
                                        cell_ib[n].slope_y_cut[intersection_count] = (point_Body[face_Body[i].node[1]].y - point_Body[face_Body[i].node[0]].y)/
                                            (pow(pow(point_Body[face_Body[i].node[1]].x - point_Body[face_Body[i].node[0]].x,2) + 
                                                 pow(point_Body[face_Body[i].node[1]].y - point_Body[face_Body[i].node[0]].y,2),0.5));
                                        //printf("\nIntersection count = %d   cell_ib[%d].face[%d]  xcut,ycut = %lf,%lf  x0,y0 = %lf,%lf   x1,y1 = %lf,%lf", 
                                        //                intersection_count, n, k, cell_ib[n].x_cut[intersection_count], cell_ib[n].y_cut[intersection_count], 
                                        //                point[cell_ib[n].face[k][0]].x, point[cell_ib[n].face[k][0]].y, point[cell_ib[n].face[k][1]].x, point[cell_ib[n].face[k][1]].y );
                                        intersection_count++;
                                        cut_face_count = 1;
                                        break;
                                    }else if(point[cell_ib[n].face[k][0]].state == 0){
                                        cell_ib[n].alpha[k] = t;
                                        cell_ib[n].x_cut[intersection_count] = point[cell_ib[n].face[k][0]].x + t*(point[cell_ib[n].face[k][1]].x - point[cell_ib[n].face[k][0]].x);
                                        cell_ib[n].y_cut[intersection_count] = point[cell_ib[n].face[k][0]].y + t*(point[cell_ib[n].face[k][1]].y - point[cell_ib[n].face[k][0]].y);
                                        cell_ib[n].slope_x_cut[intersection_count] = (point_Body[face_Body[i].node[1]].x - point_Body[face_Body[i].node[0]].x)/
                                            (pow(pow(point_Body[face_Body[i].node[1]].x - point_Body[face_Body[i].node[0]].x,2) + 
                                                 pow(point_Body[face_Body[i].node[1]].y - point_Body[face_Body[i].node[0]].y,2),0.5));
                                        cell_ib[n].slope_y_cut[intersection_count] = (point_Body[face_Body[i].node[1]].y - point_Body[face_Body[i].node[0]].y)/
                                            (pow(pow(point_Body[face_Body[i].node[1]].x - point_Body[face_Body[i].node[0]].x,2) + 
                                                 pow(point_Body[face_Body[i].node[1]].y - point_Body[face_Body[i].node[0]].y,2),0.5));
                                        //printf("\nIntersection count = %d   cell_ib[%d].face[%d]  xcut,ycut = %lf,%lf  x0,y0 = %lf,%lf   x1,y1 = %lf,%lf", 
                                        //                  intersection_count, n, k, cell_ib[n].x_cut[intersection_count], cell_ib[n].y_cut[intersection_count], 
                                        //                  point[cell_ib[n].face[k][0]].x, point[cell_ib[n].face[k][0]].y, point[cell_ib[n].face[k][1]].x, point[cell_ib[n].face[k][1]].y );
                                        intersection_count++;
                                        cut_face_count = 1;
                                        break;
                                    }
                                }
                            }
                        }
                    } 
                }
                else if(k==1 || k==3){
                    for(i=1; i<=Total_Body_Points;i++){
                        if(cut_face_count == 0){
                            if(max(point_Body[face_Body[i].node[0]].x,point_Body[face_Body[i].node[1]].x)<point[cell_ib[n].face[k][0]].x||
                                    min(point_Body[face_Body[i].node[0]].x,point_Body[face_Body[i].node[1]].x)>point[cell_ib[n].face[k][0]].x){
                                //Do Nothing   //Drop all the faces which do not intersect with the face
                            }else{
                                m = (point[cell_ib[n].face[k][0]].x-point_Body[face_Body[i].node[0]].x)/
                                    (point_Body[face_Body[i].node[1]].x-point_Body[face_Body[i].node[0]].x);
                                t = (point_Body[face_Body[i].node[0]].y - point[cell_ib[n].face[k][0]].y + 
                                        m*(point_Body[face_Body[i].node[1]].y-point_Body[face_Body[i].node[0]].y))/
                                    (point[cell_ib[n].face[k][1]].y - point[cell_ib[n].face[k][0]].y);
                                if((m>=0 && m<=1) && (t>=0 && t<=1)){
                                    if(point[cell_ib[n].face[k][0]].state == 1){
                                        cell_ib[n].alpha[k] = 1-t;
                                        cell_ib[n].x_cut[intersection_count] = point[cell_ib[n].face[k][0]].x + t*(point[cell_ib[n].face[k][1]].x - point[cell_ib[n].face[k][0]].x);
                                        cell_ib[n].y_cut[intersection_count] = point[cell_ib[n].face[k][0]].y + t*(point[cell_ib[n].face[k][1]].y - point[cell_ib[n].face[k][0]].y);
                                        cell_ib[n].slope_x_cut[intersection_count] = (point_Body[face_Body[i].node[1]].x - point_Body[face_Body[i].node[0]].x)/
                                            (pow(pow(point_Body[face_Body[i].node[1]].x - point_Body[face_Body[i].node[0]].x,2) + 
                                                 pow(point_Body[face_Body[i].node[1]].y - point_Body[face_Body[i].node[0]].y,2),0.5));
                                        cell_ib[n].slope_y_cut[intersection_count] = (point_Body[face_Body[i].node[1]].y - point_Body[face_Body[i].node[0]].y)/
                                            (pow(pow(point_Body[face_Body[i].node[1]].x - point_Body[face_Body[i].node[0]].x,2) + 
                                                 pow(point_Body[face_Body[i].node[1]].y - point_Body[face_Body[i].node[0]].y,2),0.5));
                                        //printf("\nIntersection count = %d   cell_ib[%d].face[%d]  xcut,ycut = %lf,%lf  x0,y0 = %lf,%lf   x1,y1 = %lf,%lf", 
                                        //                  intersection_count, n, k, cell_ib[n].x_cut[intersection_count], cell_ib[n].y_cut[intersection_count], 
                                        //                  point[cell_ib[n].face[k][0]].x, point[cell_ib[n].face[k][0]].y, point[cell_ib[n].face[k][1]].x, point[cell_ib[n].face[k][1]].y );
                                        intersection_count++;
                                        cut_face_count = 1;
                                        break;
                                    }else if(point[cell_ib[n].face[k][0]].state == 0){
                                        cell_ib[n].alpha[k] = t;
                                        cell_ib[n].x_cut[intersection_count] = point[cell_ib[n].face[k][0]].x + t*(point[cell_ib[n].face[k][1]].x - point[cell_ib[n].face[k][0]].x);
                                        cell_ib[n].y_cut[intersection_count] = point[cell_ib[n].face[k][0]].y + t*(point[cell_ib[n].face[k][1]].y - point[cell_ib[n].face[k][0]].y);
                                        cell_ib[n].slope_x_cut[intersection_count] = (point_Body[face_Body[i].node[1]].x - point_Body[face_Body[i].node[0]].x)/
                                            (pow(pow(point_Body[face_Body[i].node[1]].x - point_Body[face_Body[i].node[0]].x,2) + 
                                                 pow(point_Body[face_Body[i].node[1]].y - point_Body[face_Body[i].node[0]].y,2),0.5));
                                        cell_ib[n].slope_y_cut[intersection_count] = (point_Body[face_Body[i].node[1]].y - point_Body[face_Body[i].node[0]].y)/
                                            (pow(pow(point_Body[face_Body[i].node[1]].x - point_Body[face_Body[i].node[0]].x,2) + 
                                                 pow(point_Body[face_Body[i].node[1]].y - point_Body[face_Body[i].node[0]].y,2),0.5));
                                        //printf("\nIntersection count = %d   cell_ib[%d].face[%d]  xcut,ycut = %lf,%lf  x0,y0 = %lf,%lf   x1,y1 = %lf,%lf", 
                                        //                  intersection_count, n, k, cell_ib[n].x_cut[intersection_count], cell_ib[n].y_cut[intersection_count], 
                                        //                  point[cell_ib[n].face[k][0]].x, point[cell_ib[n].face[k][0]].y, point[cell_ib[n].face[k][1]].x, point[cell_ib[n].face[k][1]].y );
                                        intersection_count++;
                                        cut_face_count = 1;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                count_face++;
            }
        }
        cell_ib[n].d_cut = sqrt(pow(cell_ib[n].x_cut[1] - cell_ib[n].x_cut[0],2) + pow(cell_ib[n].y_cut[1] - cell_ib[n].y_cut[0],2));
//        if (cell_ib[n].d_cut < d_cut_tol){
//            cell_ib[n].d_cut = 0;
//        }

        if(cell_ib[n].x_cut[1] == cell_ib[n].x_cut[0] && cell_ib[n].y_cut[1] == cell_ib[n].y_cut[0]){
            cell_ib[n].slope_x_cut[0] = 0;
            cell_ib[n].slope_x_cut[1] = 0;
            cell_ib[n].slope_y_cut[0] = 0;
            cell_ib[n].slope_y_cut[1] = 0;
        }

        if((cell_ib[n].x_cut[1] - cell_ib[n].x_cut[0])*cell_ib[n].slope_x_cut[0] + 
                (cell_ib[n].y_cut[1] - cell_ib[n].y_cut[0])*cell_ib[n].slope_y_cut[0] < 0){    
            // Both direction of cut points and slopes should be in similar direction, 
            // if not then reverse the direction of slopes !

            //cell_ib[n].slope_x_cut[0] = -cell_ib[n].slope_x_cut[0];
            //cell_ib[n].slope_y_cut[0] = -cell_ib[n].slope_y_cut[0];
            //cell_ib[n].slope_x_cut[1] = -cell_ib[n].slope_x_cut[1];
            //cell_ib[n].slope_y_cut[1] = -cell_ib[n].slope_y_cut[1];
            swap(cell_ib[n].x_cut[1], cell_ib[n].x_cut[0], temp);
            swap(cell_ib[n].y_cut[1], cell_ib[n].y_cut[0], temp);
            swap(cell_ib[n].slope_x_cut[1], cell_ib[n].slope_x_cut[0], temp);
            swap(cell_ib[n].slope_y_cut[1], cell_ib[n].slope_y_cut[0], temp);
        }

        cell_ib[n].slope_push = 1*cell_ib[n].d_cut;
        //printf("\nslope xcut0,ycut0 = %lf,%lf  slope xcut1,ycut1 = %lf,%lf ", cell_ib[n].slope_x_cut[0],  cell_ib[n].slope_y_cut[0],  cell_ib[n].slope_x_cut[1],  cell_ib[n].slope_y_cut[1]);
        //printf("\nintersection_count = %d", intersection_count);
    }

    for(n=1; n<=N_ib; n++){
        for(k=0; k<=3; k++){
            //Finding the flux centers for each face of cell n
            if(k==0 || k==2){
                if(point[cell_ib[n].face[k][0]].state == 0){
                    Flux_center_x = point[cell_ib[n].face[k][0]].x + dx*cell_ib[n].alpha[k]/2;
                    Flux_center_y = point[cell_ib[n].face[k][0]].y;
                }else{
                    Flux_center_x = point[cell_ib[n].face[k][1]].x - dx*cell_ib[n].alpha[k]/2;
                    Flux_center_y = point[cell_ib[n].face[k][1]].y;
                }
            }else if(k==1 || k==3){
                if(point[cell_ib[n].face[k][0]].state == 0){
                    Flux_center_y = point[cell_ib[n].face[k][0]].y + dy*cell_ib[n].alpha[k]/2;
                    Flux_center_x = point[cell_ib[n].face[k][0]].x;
                }else{
                    Flux_center_y = point[cell_ib[n].face[k][1]].y - dy*cell_ib[n].alpha[k]/2;
                    Flux_center_x = point[cell_ib[n].face[k][1]].x;
                }
            }

            //Finding the cells from which to interpolate the properties at flux center of each face of cell n
            if(Flux_center_x > cell_ib[n].x && Flux_center_y > cell_ib[n].y){
                cell_ib[n].BI_flux_cell[k][0] = cell_ib[n].cell_no;
                cell_ib[n].BI_flux_cell[k][1] = cell_ib[n].cell_no+1;
                cell_ib[n].BI_flux_cell[k][2] = cell_ib[n].cell_no+1+cell_Nx;
                cell_ib[n].BI_flux_cell[k][3] = cell_ib[n].cell_no+cell_Nx;
            }else if(Flux_center_x > cell_ib[n].x && Flux_center_y <= cell_ib[n].y){
                cell_ib[n].BI_flux_cell[k][0] = cell_ib[n].cell_no;
                cell_ib[n].BI_flux_cell[k][1] = cell_ib[n].cell_no+1;
                cell_ib[n].BI_flux_cell[k][2] = cell_ib[n].cell_no+1-cell_Nx;
                cell_ib[n].BI_flux_cell[k][3] = cell_ib[n].cell_no-cell_Nx;
            }else if(Flux_center_x <= cell_ib[n].x && Flux_center_y > cell_ib[n].y){
                cell_ib[n].BI_flux_cell[k][0] = cell_ib[n].cell_no;
                cell_ib[n].BI_flux_cell[k][1] = cell_ib[n].cell_no-1;
                cell_ib[n].BI_flux_cell[k][2] = cell_ib[n].cell_no-1+cell_Nx;
                cell_ib[n].BI_flux_cell[k][3] = cell_ib[n].cell_no+cell_Nx;
            }else if(Flux_center_x <= cell_ib[n].x && Flux_center_y <= cell_ib[n].y){
                cell_ib[n].BI_flux_cell[k][0] = cell_ib[n].cell_no;
                cell_ib[n].BI_flux_cell[k][1] = cell_ib[n].cell_no-1;
                cell_ib[n].BI_flux_cell[k][2] = cell_ib[n].cell_no-1-cell_Nx;
                cell_ib[n].BI_flux_cell[k][3] = cell_ib[n].cell_no-cell_Nx;
            }

            //Find the dx1, dx2, dy1, dy2 for interpolation of fluxes correctly
            cell_ib[n].dx1[k] = fabs(Flux_center_x - cell[cell_ib[n].BI_flux_cell[k][2]].x);
            cell_ib[n].dx2[k] = fabs(Flux_center_x - cell[cell_ib[n].BI_flux_cell[k][3]].x);
            cell_ib[n].dy1[k] = fabs(Flux_center_y - cell[cell_ib[n].BI_flux_cell[k][2]].y);
            cell_ib[n].dy2[k] = fabs(Flux_center_y - cell[cell_ib[n].BI_flux_cell[k][1]].y);

           //  printf("\ndx1 = %lf, dx2 = %lf, dy1 = %lf, dy2 = %lf", cell_ib[n].dx1[k], cell_ib[n].dx2[k], cell_ib[n].dy1[k], cell_ib[n].dy2[k]);
        }
    }

    //----------------------Setting theta values for cut points on faces for each ib cell--------------------------

    
    for (n=1; n<=N_ib; n++) {
        for(k=0; k<=1; k++){
            if(cell_ib[n].y_cut[k] > yo){
                cell_ib[n].theta_cut[k] = acos((cell_ib[n].x_cut[k] - xo)/Radius);
            }else if(cell_ib[n].y_cut[k] == yo){
                if(cell_ib[n].y >= yo){
                    if((cell_ib[n].x_cut[k] - xo)/Radius >= 1.0){
                        cell_ib[n].theta_cut[k] = acos(1); 
                    }else if((cell_ib[n].x_cut[k] - xo)/Radius <= -1.0){ 
                        cell_ib[n].theta_cut[k] = acos(-1); 
                    }else{
                        cell_ib[n].theta_cut[k] = acos((cell_ib[n].x_cut[k] - xo)/Radius);
                    }
                }else{
                    if((cell_ib[n].x_cut[k] - xo)/Radius >= 1.0){
                        cell_ib[n].theta_cut[k] = 2*pi - acos(1); 
                    }else if((cell_ib[n].x_cut[k] - xo)/Radius <= -1.0){ 
                        cell_ib[n].theta_cut[k] = 2*pi - acos(-1); 
                    }else{
                    cell_ib[n].theta_cut[k] = 2*pi - acos((cell_ib[n].x_cut[k] - xo)/Radius);
                    }
                }
            }else{
                cell_ib[n].theta_cut[k] = 2*pi - acos((cell_ib[n].x_cut[k] - xo)/Radius);
            }
        }
    }

    cell_ib[0].theta = -1;
    for (n=1; n<=N_ib; n++) {
        cell_ib[n].N_local_division = cell_ib[n].d_cut/pow(pow(dx,2) + pow(dy,2),0.5)*N_division + 1;
        cell_ib[n].theta = (cell_ib[n].theta_cut[0] + cell_ib[n].theta_cut[1])/2;
//       printf("\ncell = %d, theta_0 = %lf, theta_1 = %lf, Local_division = %d", n, cell_ib[n].theta_cut[0], cell_ib[n].theta_cut[1], cell_ib[n].N_local_division);
    }
    
    //Defining the volume fractions of the cell
    //Defining cell types:  cell_type = 0:      Triangular
    //                      cell_type = 1:      Trapezoidal
    //                      cell_type = 2:      Pentagon type
    //                      cell_type = 4:      Regular
    
//    printf("\nPrinting immersed cell types:\n");
    for(n=1; n<=N_ib; n++){
        cell_ib[n].cell_vol_type = floor(cell_ib[n].alpha[0]) + floor(cell_ib[n].alpha[1]) +
                                   floor(cell_ib[n].alpha[2]) + floor(cell_ib[n].alpha[3]);
//        printf("IB cell = %d\tCell Type = %d\n", n, cell_ib[n].cell_vol_type);
    }

//    printf("\nPrinting immersed cell Volume Fractions:\n");
    for(n=1; n<=N_ib; n++){
        if(cell_ib[n].cell_vol_type == 0){
            count = 0;
            for(k=0; k<4; k++){
                b[count] = 0;
                if(cell_ib[n].alpha[k] != 0){
                    b[count] = cell_ib[n].alpha[k];
                    count++;     
                } 
            }
            cell_ib[n].vol_fr = 0.5*(b[0]*b[1]); 
        }
        else if(cell_ib[n].cell_vol_type == 1){
            b[0] = 0;
            b[1] = 1;
            for(k=0; k<4; k++){
                if(cell_ib[n].alpha[k] != 0 && cell_ib[n].alpha[k] != 1){
                    b[0] = b[0] + cell_ib[n].alpha[k];
                } 
            }
            cell_ib[n].vol_fr = 0.5*(b[0]*b[1]); 
        }
        else if(cell_ib[n].cell_vol_type == 2){
            count = 0;
            for(k=0; k<4; k++){
                b[count] = 0;
                if(cell_ib[n].alpha[k] != 1){
                    b[count] = 1 - cell_ib[n].alpha[k];
                    count++;
                } 
            }
            cell_ib[n].vol_fr = 1 - 0.5*(b[0]*b[1]); 
        }
        else if(cell_ib[n].cell_vol_type == 3 ||cell_ib[n].cell_vol_type == 4){
            cell_ib[n].vol_fr = 1;
        }
        else{
            printf("Something wrong with cell volume. Check ibm_mesh.c.");
        }
//        printf("IB cell = %d\tCell Type = %d\td_cut = %lf\tVol Fraction = %lf\n ", n, cell_ib[n].cell_vol_type, cell_ib[n].d_cut, cell_ib[n].vol_fr);
    }   
        
    for(n=1; n<=Total_cells; n++){
        if(cell[n].state == 2){
            k = cell[n].ib_cell_no;
            cell[n].vol_fr = cell_ib[k].vol_fr; 
        }else{
            cell[n].vol_fr = 1;
        }
    }
}
