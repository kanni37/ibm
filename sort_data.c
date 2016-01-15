#include "ibm.h"
#include "functions.h"

void sort_data()
{
    int k, n;

    qsort(cell_ib, N_ib + 1, sizeof(cell_ib[0]), comp_ib_cells);

    for(n=1; n<=N_ib; n++){
        k = cell_ib[n].cell_no;
        cell[k].ib_cell_no = n;
    }
}
