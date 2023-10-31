#include <stdio.h>
#include <stdlib.h>

int main() {

    /* TEST east */

    int n = 4;
    int m = 4;
    int o = 4;
    
    int dim_chunk1 = 4;
    int dim_chunk2 = 4;
    int dim_chunk3 = 4;

    /*
    int *tab = malloc(sizeof(int)*n*m*o);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < o; ++k) {
                tab[i + j*n + k*n*m] = i + j*n + k*n*m;
            }
        }
    }

    int *east = malloc(sizeof(int)*m*o);

    for (int j = 0; j < dim_chunk2; ++j) {
        for (int k = 0; k < dim_chunk3; ++k) {
            east[j*dim_chunk3 + k] = tab[dim_chunk1-1 + j*dim_chunk1 + k*dim_chunk1*dim_chunk2];
            printf("%d ", east[j*dim_chunk3 + k]);
        }
        printf("\n");
    }
    printf("\n");
    */

    /* TEST xy_plan */

    /*
    int *tab = calloc(n*m*o, sizeof(int));

    int T_i = 0;
    int T_j = 0;
    int T_k = 1; // 0 cas particulier
    int inside = 1; // 0 pour tout, 1 pour l'intérieur

    
    if (T_k == 0) // we do not modify the z = 0 plane: it is maintained at constant temperature via water-cooling
        return 0;

    for (int j = inside; j < dim_chunk2-inside; ++j) {

        for (int i = inside; i < dim_chunk1-inside; ++i) {
            int u = T_k*dim_chunk2*dim_chunk1 + j * dim_chunk1 + i;
            tab[u] += 1;
            ++T_i;
        }

        T_i -= dim_chunk1 + 2*inside;
        ++T_j;
    }

    for (int j = 0; j < dim_chunk2; ++j) {
        for (int i = 0; i < dim_chunk1; ++i) {  
            printf("%d ", tab[i + j*dim_chunk1 + T_k*dim_chunk2*dim_chunk1]);
        }
        printf("\n");
    }
    printf("\n");
    */
    
    return 0;
}