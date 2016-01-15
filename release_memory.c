#include "ibm.h"

void release_memory()
{
    free(cell);
    free(cell_ib);
    free(point);
    free(point_Body);
    free(face_Body);
}
