#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {

    int n = 4;
    int m = 4;
    int o = 4;
    
    int dim_chunk0 = 4;
    int dim_chunk1 = 4;
    int dim_chunk2 = 4;

    /* TEST east */
    
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

    for (int j = 0; j < dim_chunk1; ++j) {
        for (int k = 0; k < dim_chunk2; ++k) {
            east[j*dim_chunk2 + k] = tab[dim_chunk0-1 + j*dim_chunk0 + k*dim_chunk0*dim_chunk1];
            printf("%d ", east[j*dim_chunk2 + k]);
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

    for (int j = inside; j < dim_chunk1-inside; ++j) {

        for (int i = inside; i < dim_chunk0-inside; ++i) {
            int u = T_k*dim_chunk1*dim_chunk0 + j * dim_chunk0 + i;
            tab[u] += 1;
            ++T_i;
        }

        T_i -= dim_chunk0 + 2*inside;
        ++T_j;
    }

    for (int j = 0; j < dim_chunk1; ++j) {
        for (int i = 0; i < dim_chunk0; ++i) {  
            printf("%d ", tab[i + j*dim_chunk0 + T_k*dim_chunk1*dim_chunk0]);
        }
        printf("\n");
    }
    printf("\n");
    */

    /* TEST halo */

    /*
    int *halo = calloc(sizeof(int),(dim_chunk0+2)*(dim_chunk1+2)*(dim_chunk2+2));
    int *west = malloc(sizeof(int)*dim_chunk1*dim_chunk2);

    for (int j = 0; j < dim_chunk1; ++j) {
        for (int k = 0; k < dim_chunk2; ++k) {
            west[j*dim_chunk2 + k] = j*dim_chunk2 + k;
        }
    }

    for (int j = 1; j < dim_chunk1+1; ++j) {
        for (int k = 1; k < dim_chunk2+1; ++k) {
            halo[j*(dim_chunk0+2) + k*(dim_chunk0+2)*(dim_chunk1+2)] = west[(j-1)*dim_chunk2 + k-1];
            printf("%d ",west[(j-1)*dim_chunk2 + k]);
        }
        printf("\n");
    }

    for (int j = 0; j < dim_chunk1+2; ++j) {
        for (int k = 0; k < dim_chunk2+2; ++k) {
            printf("%d ",halo[j*(dim_chunk0+2) + k*(dim_chunk0+2)*(dim_chunk1+2)]);
        }
        printf("\n");
    }
    printf("\n");
    */

    /*
    n = 168;
    
    int i = floor(cbrt(n));
    while (n%i != 0){
        i = i-1;
    }
    int j = floor(sqrt(n/i));
    while ((n/i)%j != 0){
        j = j-1;
    }
    int k = floor(n/i/j);

    int *tab = malloc(3*sizeof(int));
    if (tab == NULL){
	perror("malloc");
	exit(1);	
    }
    printf(" k = %d ",k);
    printf("j= %d ",j);
    printf("i=%d\n ",i);
    */


    int *tab = malloc(sizeof(int)*10);
    free(tab);
    free(tab);
    

    return 0;
}