#include "ibm.h"

// Defining functions for use in main file (ibm.c)

int move(double displacement_x,double displacement_y, struct points *point_OUT,int Total_Points_OUT)
{
    int i;
    for(i=1;i<=Total_Points_OUT;i++){                   // traversing each point one by one
        point_OUT[i].x=point_OUT[i].x+displacement_x;	// giving x displacement to body
		point_OUT[i].y=point_OUT[i].y+displacement_y;	// giving y displacement to body
	 }
    return 0;
}

int gcell(int i, int j){
    return ((j-1)*(Nx-1) + i);
}

int gpoint(int i, int j){
    return ((j-1)*Nx + i);
}

int comp_ib_cells(const void *a, const void *b)
{
    const struct ib_cells *sa = a;
    const struct ib_cells *sb = b;
    if (sa->theta > sb->theta) {
        return 1;
    } else if (sa->theta < sb->theta) {
        return -1;
    } else {
        return 0;
    }
}
