#include "ibm.h"
#include "functions.h"

void place_body(int Immersed_body_type)
{
    switch (Immersed_body_type){
        case 1:
            place_nothing();
            break;
        case 2:
            place_circle();
            break;
        case 3:
            place_square();
            break;
        case 4:
            place_airfoil();
            break;
        default:
            printf("Error, bad Immersed_body_type input, quitting. Please check place_body.c");
            break;
    }
}

void place_circle()
{
    //-------------------------Defining Lagrangian Grid (Solid Body)-----------------------
    int n;
    
    xo = L/3;
    yo = H/2;
    
    printf("\nxo = %f,yo = %f", xo, yo);
    
    Total_Body_Points = 1023;
    //double Radius = 0.5
    
    printf("\n\nRadius = %lf", Radius);
    
    point_Body = malloc((Total_Body_Points+1)*sizeof(struct points));
    
    for(n=1; n<=Total_Body_Points; n++){
        point_Body[n].x = Radius*cos(n*2*pi/Total_Body_Points);
        point_Body[n].y = Radius*sin(n*2*pi/Total_Body_Points);
    }
    
    face_Body = malloc((Total_Body_Points+1)*sizeof(struct faces));
    
    for(n=1; n<=Total_Body_Points;n++){
        if(n==Total_Body_Points){
            face_Body[n].node[0] = n;
            face_Body[n].node[1] = 1;
        }else{
            face_Body[n].node[0] = n;
            face_Body[n].node[1] = n+1;
        }
    }
    
    move(xo, yo, point_Body, Total_Body_Points);
    
    xmax = point_Body[1].x;
    xmin = point_Body[1].x;
    ymax = point_Body[1].y;
    ymin = point_Body[1].y;
    
    for(n=1; n<=Total_Body_Points; n++){
        xmax = max(xmax,point_Body[n].x);
        xmin = min(xmin,point_Body[n].x);
        ymax = max(ymax,point_Body[n].y);
        ymin = min(ymin,point_Body[n].y);
    }
    
    xmax = xmax + 2*dx;
    xmin = xmin - 2*dx;
    ymax = ymax + 2*dy;
    ymin = ymin - 2*dy;
    
    printf("\n\nBounding Box Dimension:");
    printf("\nxmin = %f\txmax = %f", xmin, xmax);
    printf("\nymin = %f\tymax = %f\n", ymin, ymax);
}

void place_square()
{
    //-------------------------Defining Lagrangian Grid (Solid Body)-----------------------
    int n;
    
    xo = L/2 + dx/2;
    yo = H/2 + dy/2;
    
    printf("\nxo = %f,yo = %f", xo, yo);
    
    double side = 0.1;
    Total_Body_Points = 720;
    
    point_Body = malloc((Total_Body_Points+1)*sizeof(struct points));
    
    for(n=1; n<=Total_Body_Points; n++){
         if(n<=Total_Body_Points/4){
         point_Body[n].x = 0.5*side;
         point_Body[n].y = (-0.5+n*4/Total_Body_Points)*side;
         } else if(n>Total_Body_Points/4 && n<=Total_Body_Points/2){
         point_Body[n].x = (0.5 - (n-Total_Body_Points/4)*4/Total_Body_Points)*side;
         point_Body[n].y = 0.5*side;
         } else if(n>Total_Body_Points/2 && n<=3*Total_Body_Points/4){
         point_Body[n].x = -0.5*side;
         point_Body[n].y = (0.5 - (n-Total_Body_Points/2)*4/Total_Body_Points)*side;
         } else{
         point_Body[n].x = (-0.5 + (n-3*Total_Body_Points/4)*4/Total_Body_Points)*side;
         point_Body[n].y = -0.5*side;
         }
    }
    
    face_Body = malloc((Total_Body_Points+1)*sizeof(struct faces));
    
    for(n=1; n<=Total_Body_Points;n++){
        if(n==Total_Body_Points){
            face_Body[n].node[0] = n;
            face_Body[n].node[1] = 1;
        }else{
            face_Body[n].node[0] = n;
            face_Body[n].node[1] = n+1;
        }
    }
    
    move(xo, yo, point_Body, Total_Body_Points);
    
    xmax = point_Body[1].x;
    xmin = point_Body[1].x;
    ymax = point_Body[1].y;
    ymin = point_Body[1].y;
    
    for(n=1; n<=Total_Body_Points; n++){
        xmax = max(xmax,point_Body[n].x);
        xmin = min(xmin,point_Body[n].x);
        ymax = max(ymax,point_Body[n].y);
        ymin = min(ymin,point_Body[n].y);
    }
    
    xmax = xmax + 2*dx;
    xmin = xmin - 2*dx;
    ymax = ymax + 2*dy;
    ymin = ymin - 2*dy;
    
    printf("\nBounding Box Dimension:");
    printf("\nxmin = %f\txmax = %f", xmin, xmax);
    printf("\nymin = %f\tymax = %f\n", ymin, ymax);

}

void place_nothing()
{
    //-------------------------Defining Lagrangian Grid (Solid Body)-----------------------
    int n;
    
    xo = L/2;
    yo = H/2;
    
    printf("\nxo = %f,yo = %f", xo, yo);
    
    Total_Body_Points = 1;
    
    point_Body = malloc((Total_Body_Points+1)*sizeof(struct points));
    
    for(n=1; n<=Total_Body_Points; n++){
        point_Body[n].x = 0; 
        point_Body[n].y = 0; 
    }
    
    face_Body = malloc((Total_Body_Points+1)*sizeof(struct faces));
    
    for(n=1; n<=Total_Body_Points;n++){
        if(n==Total_Body_Points){
            face_Body[n].node[0] = n;
            face_Body[n].node[1] = 1;
        }else{
            face_Body[n].node[0] = n;
            face_Body[n].node[1] = n+1;
        }
    }
    
    move(4*xo,4*yo, point_Body, Total_Body_Points);
    
    xmax = point_Body[1].x;
    xmin = point_Body[1].x;
    ymax = point_Body[1].y;
    ymin = point_Body[1].y;
    
    for(n=1; n<=Total_Body_Points; n++){
        xmax = max(xmax,point_Body[n].x);
        xmin = min(xmin,point_Body[n].x);
        ymax = max(ymax,point_Body[n].y);
        ymin = min(ymin,point_Body[n].y);
    }
    
    xmax = xmax + 2*dx;
    xmin = xmin - 2*dx;
    ymax = ymax + 2*dy;
    ymin = ymin - 2*dy;
    
    printf("\n\nBounding Box Dimension:");
    printf("\nxmin = %f\txmax = %f", xmin, xmax);
    printf("\nymin = %f\tymax = %f\n", ymin, ymax);
}

void place_airfoil()
{
    //-------------------------Defining Lagrangian Grid (Solid Body)-----------------------
    int n;
    double a = 0.06;
    double c = 0.05;
    double n_lambda = 23.0/12.0;
    
    Total_Body_Points = 120;
   
    double theta[Total_Body_Points];  
    double complex zc[Total_Body_Points]; 
    double complex T[Total_Body_Points]; 
    double complex z[Total_Body_Points]; 

    xo = L/2;
    yo = H/2;
    
    printf("\nxo = %f,yo = %f", xo, yo);
    
    point_Body = malloc((Total_Body_Points+1)*sizeof(struct points));
    
    for(n=1; n<=Total_Body_Points; n++){
        theta[n] = n*2*pi/Total_Body_Points;
    }
   
    for(n=1; n<=Total_Body_Points; n++){
        zc[n] = a*cexp(theta[n]*I) + (c - a);
    }

    for(n=1; n<=Total_Body_Points; n++){
        T[n] = cpow(((zc[n] - c)/(zc[n] + c)), n_lambda);
    }
   
    for(n=1; n<=Total_Body_Points; n++){
        z[n] = n_lambda*c*(1 + T[n])/(1 - T[n]);
    }
    
    for(n=1; n<=Total_Body_Points; n++){
        point_Body[n].x = creal(z[n]);
        point_Body[n].y = cimag(z[n]); 
    }
    
    face_Body = malloc((Total_Body_Points+1)*sizeof(struct faces));
    
    for(n=1; n<=Total_Body_Points;n++){
        if(n==Total_Body_Points){
            face_Body[n].node[0] = n;
            face_Body[n].node[1] = 1;
        }else{
            face_Body[n].node[0] = n;
            face_Body[n].node[1] = n+1;
        }
    }
    
    move(xo, yo, point_Body, Total_Body_Points);
    
    xmax = point_Body[1].x;
    xmin = point_Body[1].x;
    ymax = point_Body[1].y;
    ymin = point_Body[1].y;
    
    for(n=1; n<=Total_Body_Points; n++){
        xmax = max(xmax,point_Body[n].x);
        xmin = min(xmin,point_Body[n].x);
        ymax = max(ymax,point_Body[n].y);
        ymin = min(ymin,point_Body[n].y);
    }
    
    xmax = xmax + 2*dx;
    xmin = xmin - 2*dx;
    ymax = ymax + 2*dy;
    ymin = ymin - 2*dy;
    
    printf("\n\nBounding Box Dimension:");
    printf("\nxmin = %f\txmax = %f", xmin, xmax);
    printf("\nymin = %f\tymax = %f\n", ymin, ymax);
    

    FILE *fp = NULL;
    fp = fopen("airfoil_data.dat", "w");
    fprintf(fp,"VARIABLES = \"x\", \"y\"\n");
    for(n=1; n<=Total_Body_Points; n++){
        //fprintf(fp,"%lf\t%lf\n", creal(T[n]), cimag(T[n]));
        //fprintf(fp,"%lf\t%lf\n", creal(zc[n]), cimag(zc[n]));
        fprintf(fp,"%lf\t%lf\n", point_Body[n].x, point_Body[n].y);
    }
    fclose(fp);
}
