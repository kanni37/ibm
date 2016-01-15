#include "ibm.h"

int is_point_in_polygon(double point_x, double point_y)
{
    double m, t;
    int count, ray_above, ray_below;
    int i;
    
    //----------------------Identifying Solid and Fluid and Immersed Points------------------
    count=0;
    ray_above=0;
    ray_below=0;
        
    if(point_x<xmin||point_x>xmax||point_y<ymin||point_y>ymax){
        count=-2;
        //Do Nothing
    }else{
        for(i=1; i<=Total_Body_Points; i++){
            if(max(point_Body[face_Body[i].node[0]].x,point_Body[face_Body[i].node[1]].x)<point_x||
               max(point_Body[face_Body[i].node[0]].y,point_Body[face_Body[i].node[1]].y)<point_y||
               min(point_Body[face_Body[i].node[0]].y,point_Body[face_Body[i].node[1]].y)>point_y){
                //Do Nothing   //Trying to drop all the faces which do not intersect with the semi infinite ray
            }else{
                if(point_Body[face_Body[i].node[0]].y!=point_Body[face_Body[i].node[1]].y){
                    m = (point_y-point_Body[face_Body[i].node[0]].y)/(point_Body[face_Body[i].node[1]].y-point_Body[face_Body[i].node[0]].y);
                    t = point_Body[face_Body[i].node[0]].x + m*(point_Body[face_Body[i].node[1]].x-point_Body[face_Body[i].node[0]].x);
                    
                    if(t==point_x){
                        count=-1;
                        break;
                    }
                    if(t>point_x){
                        if(m>0 && m<1){
                            count++;
                        }else{
                            if(m*point_Body[face_Body[i].node[0]].y+(1-m)*point_Body[face_Body[i].node[1]].y-point_y > 0){
                                ray_above++;
                            }else{
                                ray_below++;
                            }
                        }
                    }
                }else if((point_Body[face_Body[i].node[1]].x-point_x)*(point_Body[face_Body[i].node[0]].x-point_x) <=0){
                    count = -1;
                    break;
                }
            }
        }
        if(count!=-1){
            count=count+(ray_above<=ray_below?ray_above:ray_below);
        }
    }
    if((count-2*(count/2))!=0){         // Count: odd = point inside : even = point outside 
        return 1;                       // Point inside
    }else{
        return 0;                       // Point outside
    }
}
