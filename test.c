#include <stdio.h>
#include <stdlib.h>

int main() {
    
    int n = 2;
    int m = 2;
    int o = 2;
    
    int dim_chunk1 = 2;
    int dim_chunk2 = 2;
    int dim_chunk3 = 2;

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


    return 0;
}