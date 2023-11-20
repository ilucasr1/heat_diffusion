#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <mpi.h>

/* AUTHOR : Charles Bouillaguet <charles.bouillaguet@lip6.fr>
   USAGE  : compile with -lm (and why not -O3)
            redirect the standard output to a text file
            gcc heatsink.c -O3 -lm -o heatsink
            ./heatsink > steady_state.txt
            then run the indicated python script for graphical rendering

   DISCLAIMER : this code does not claim to an absolute realism.
                this code could be obviously improved, but has been written so as
				to make as clear as possible the physics principle of the simulation.
*/

/* one can change the matter of the heatsink, its size, the power of the CPU, etc. */
#define ALUMINIUM
#define FAST         /* MEDIUM is faster, and FAST is even faster (for debugging) */
#define DUMP_STEADY_STATE

const double L = 0.15;      /* length (x) of the heatsink (m) */
const double l = 0.12;      /* height (z) of the heatsink (m) */
const double E = 0.008;     /* width (y) of the heatsink (m) */
const double watercooling_T = 20;   /* temperature of the fluid for water-cooling, (°C) */
const double CPU_TDP = 280; /* power dissipated by the CPU (W) */

/* dl: "spatial step" for simulation (m) */
/* dt: "time step" for simulation (s) */
#ifdef FAST
double dl = 0.004;
double dt = 0.004;
#endif

#ifdef MEDIUM
double dl = 0.002;
double dt = 0.002;
#endif

#ifdef NORMAL
double dl = 0.001;
double dt = 0.001;
#endif

#ifdef CHALLENGE
double dl = 0.0001;
double dt = 0.00001;
#endif

/* sink_heat_capacity: specific heat capacity of the heatsink (J / kg / K) */
/* sink_density: density of the heatsink (kg / m^3) */
/* sink_conductivity: thermal conductivity of the heatsink (W / m / K) */
/* euros_per_kg: price of the matter by kilogram */
#ifdef ALUMINIUM
double sink_heat_capacity = 897;
double sink_density = 2710;
double sink_conductivity = 237;
double euros_per_kg = 1.594;
#endif

#ifdef COPPER
double sink_heat_capacity = 385;
double sink_density = 8960;
double sink_conductivity = 390;
double euros_per_kg = 5.469;
#endif

#ifdef GOLD
double sink_heat_capacity = 128;
double sink_density = 19300;
double sink_conductivity = 317;
double euros_per_kg = 47000;
#endif

#ifdef IRON
double sink_heat_capacity = 444;
double sink_density = 7860;
double sink_conductivity = 80;
double euros_per_kg = 0.083;
#endif

const double Stefan_Boltzmann = 5.6703e-8;  /* (W / m^2 / K^4), radiation of black body */
const double heat_transfer_coefficient = 10;    /* coefficient of thermal convection (W / m^2 / K) */
double CPU_surface;

int r = 0; // pour des test
int it = 0; //pour des test


/* 
 * Return True if the CPU is in contact with the heatsink at the point (x,y).
 * This describes an AMD EPYC "Rome".
 */
static inline bool CPU_shape(double x, double y)
{
    x -= (L - 0.0754) / 2;
    y -= (l - 0.0585) / 2;
    bool small_y_ok = (y > 0.015 && y < 0.025) || (y > 0.0337 && y < 0.0437);
    bool small_x_ok = (x > 0.0113 && x < 0.0186) || (x > 0.0193 && x < 0.0266)
        || (x > 0.0485 && x < 0.0558) || (x > 0.0566 && x < 0.0639);
    bool big_ok = (x > 0.03 && x < 0.045 && y > 0.0155 && y < 0.0435);
    return big_ok || (small_x_ok && small_y_ok);
}

/* returns the total area of the surface of contact between CPU and heatsink (in m^2) */
double CPU_contact_surface()
{
    double S = 0;
    for (double x = dl / 2; x < L; x += dl)
        for (double y = dl / 2; y < l; y += dl)
            if (CPU_shape(x, y))
                S += dl * dl;
    return S;
}

/* Returns the new temperature of the cell (i, j, k). For this, there is an access to neighbor
 * cells (left, right, top, bottom, front, back), except if (i, j, k) is on the external surface. */
static inline double update_temperature(const double c, const double f, const double b, const double s, 
                                        const double no, const double w, const double e, 
                                        int n, int m, int o, int i, int j, int k)
{
		/* quantity of thermal energy that must be brought to a cell to make it heat up by 1°C */
    // if (k == 15 && j == 0 && i == 37) {
    //     fprintf(stderr, "it %d my_rank %d c %f f %f b %f s %f n %f w %f e %f\n", 
    //                     it, r, c, f, b, s, no, w, e);
    // }
    const double cell_heat_capacity = sink_heat_capacity * sink_density * dl * dl * dl; /* J.K */
    const double dl2 = dl * dl;
    double thermal_flux = 0;

    if (i > 0)
        thermal_flux += (w - c) * sink_conductivity * dl; /* neighbor x-1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(c, 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (c - watercooling_T);
    }

    if (i < n - 1)
        thermal_flux += (e - c) * sink_conductivity * dl; /* neighbor x+1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(c, 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (c - watercooling_T);
    }

    if (j > 0)
        thermal_flux += (s - c) * sink_conductivity * dl; /* neighbor y-1 */
    else {
        /* Bottom cell: does it receive it from the CPU ? */
	if (CPU_shape(i * dl, k * dl))
            thermal_flux += CPU_TDP / CPU_surface * dl2;
        else {
            thermal_flux -= Stefan_Boltzmann * dl2 * pow(c, 4);
            thermal_flux -= heat_transfer_coefficient * dl2 * (c - watercooling_T);
        }
    }

    if (j < m - 1)
        thermal_flux += (no - c) * sink_conductivity * dl; /* neighbor y+1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(c, 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (c - watercooling_T);
    }

    if (k > 0)
        thermal_flux += (f - c) * sink_conductivity * dl; /* neighbor z-1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(c, 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (c - watercooling_T);
    }

    if (k < o - 1)
        thermal_flux += (b - c) * sink_conductivity * dl; /* neighbor z+1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(c, 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (c - watercooling_T);
    }

    /* adjust temperature depending on the heat flux */
    return c + thermal_flux * dt / cell_heat_capacity;
}

/* Run the simulation on the k-th xy plane.
 * v is the index of the start of the k-th xy plane in the arrays T and R. */
static inline void do_xy_plane(const double *T, double *R, 
                               int n, int m, int o, int k, int *dim_chunk, int my_rank, int inside)
{
    int T_i = (my_rank*dim_chunk[0])%n;
    int T_j = ((my_rank*dim_chunk[0]/n)*dim_chunk[1])%m;
    int T_k = ( ( (my_rank*dim_chunk[0]/n)*dim_chunk[1] )/m )*dim_chunk[2] + k;

    if (T_k == 0) // we do not modify the z = 0 plane: it is maintained at constant temperature via water-cooling
        return;

    int u;
    for (int j = inside; j < dim_chunk[1]-inside; ++j) {   // y
        for (int i = inside; i < dim_chunk[0]-inside; ++i) {   // x
            u = k*dim_chunk[1]*dim_chunk[0] + j * dim_chunk[0] + i;
            R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], T[u + dim_chunk[0]*dim_chunk[1]],
                                        T[u - dim_chunk[0]], T[u + dim_chunk[0]], T[u - 1], T[u + 1], n, m, o, T_i, T_j, T_k);
            ++T_i;
        }
        T_i -= dim_chunk[0] + 2*inside;
        ++T_j;
    }
}

// static inline void do_xz_plane(const double *T, double *R, int *edges, 
//                                int n, int m, int o, int j, int *dim_chunk, int p, int inside) // A FAIRE
// {
//     int T_i = (p*dim_chunk[0])%n;
//     int T_j = ((p*dim_chunk[0]/n)*dim_chunk[1])%m + j;
//     int T_k = ( ( (p*dim_chunk[0]/n)*dim_chunk[1] )/m )*dim_chunk[2];

//     for (int k = inside; k < dim_chunk[2]-inside; ++k) {   // z
//         if (T_k == 0){
//             ++T_k;
//             continue;
//         }
//         for (int i = inside; i < dim_chunk[0]-inside; ++i) {   // x
//             int u = k*dim_chunk[1]*dim_chunk[0] + j * dim_chunk[0] + i;
//             R[u] = update_temperature(T, u, n, m, o, T_i, T_j, T_k);
//             ++T_i;
//         }
//         T_i -= dim_chunk[0] + 2*inside;
//         ++T_k;
//     }
// }

// static inline void do_yz_plane(const double *T, double *R,  int *edges, 
//                                int n, int m, int o, int i, int *dim_chunk, int p, int inside) // A FAIRE
// {
//     int T_i = (p*dim_chunk[0])%n + i;
//     int T_j = ((p*dim_chunk[0]/n)*dim_chunk[1])%m;
//     int T_k = ( ( (p*dim_chunk[0]/n)*dim_chunk[1] )/m )*dim_chunk[2];

//     for (int k = inside; k < dim_chunk[2]-inside; ++k) {   // z
//         if (T_k == 0){
//             ++T_k;
//             continue;
//         }
//         for (int j = inside; j < dim_chunk[1]-inside; ++j) {   // y
//             int u = k*dim_chunk[1]*dim_chunk[0] + j * dim_chunk[0] + i;
//             R[u] = update_temperature(T, u, n, m, o, T_i, T_j, T_k);
//             ++T_j;
//         }
//         T_j -= dim_chunk[1] + 2*inside;
//         ++T_k;
//     }
// }


static inline void do_halo(const double *T, double *R,  double *front_r, double *back_r, 
                            double *south_r, double *north_r, double *west_r, double *east_r, 
                            int n, int m, int o, int *dim_chunk, int my_rank) {
/**
 * @param T chunk at time t
 * @param R chunk at time t+1
 * @param n dim of the whole data on x
 * @param m dim of the whole data on y
 * @param o dim of the whole data on z
 * @param my_rank current processor's rank
 */

    int T_i = (my_rank*dim_chunk[0])%n;
    int T_j = ((my_rank*dim_chunk[0]/n)*dim_chunk[1])%m;
    int T_k = ( ( (my_rank*dim_chunk[0]/n)*dim_chunk[1] )/m )*dim_chunk[2];

    int u;
    double f;
    double b;
    double s;
    double no;
    double w;
    double e;

    // FACE FRONT
    // 3 cas : centre / arete / coins
    // cas centre :
    if (T_k == 0) goto end_front;

    for (int j = 1; j < dim_chunk[1] - 1; ++j) {
        for (int i = 1; i < dim_chunk[0] - 1; ++i) {
            u = j*dim_chunk[0] + i;
            f = (front_r == NULL) ? 0 : front_r[u];
            b = (dim_chunk[2] > 1) ? T[u + dim_chunk[0]*dim_chunk[1]] :
                ((back_r == NULL) ? 0 : back_r[u]);

            R[u] = update_temperature(T[u], f, b, T[u - dim_chunk[0]],
            T[u + dim_chunk[0]], T[u - 1], T[u + 1], n, m, o, T_i + i, T_j + j, T_k);
        }
    }

    // cas aretes :
    // arete bas et haut
    for (int i = 1; i < dim_chunk[0]-1; ++i) {

        //j = 0
        u = i;
        f = (front_r == NULL) ? 0 : front_r[u];
        b = (dim_chunk[2] > 1) ? T[u + dim_chunk[0]*dim_chunk[1]] :
            ((back_r == NULL) ? 0 : back_r[u]);
        s = (south_r == NULL) ? 0 : south_r[i];
        no = (dim_chunk[1] > 1) ? T[u + dim_chunk[0]] :
            ((north_r == NULL) ? 0 : north_r[i]);
        
        R[u] = update_temperature(T[u], f, b, s, no, 
                                  T[u - 1], T[u + 1], n, m, o, T_i + i, T_j, T_k);
                                    
        //j = dim_chunk[1]-1
        u = dim_chunk[0]*(dim_chunk[1]-1) + i;
        f = (front_r == NULL) ? 0 : front_r[u];
        b = (dim_chunk[2] > 1) ? T[u + dim_chunk[0]*dim_chunk[1]] :
            ((back_r == NULL) ? 0 : back_r[u]);
        s = (dim_chunk[1] > 1) ? T[u - dim_chunk[0]] :
            ((south_r == NULL) ? 0 : south_r[i]);
        no = (north_r == NULL) ? 0 : north_r[i];
        
        R[u] = update_temperature(T[u], f, b, s, no, T[u - 1], T[u + 1], 
                                n, m, o, T_i + i, T_j + dim_chunk[1] - 1, T_k);
    }
    
    //arete gauche et droite
    for (int j = 1; j < dim_chunk[1]-1; ++j) {

        //i = 0
        u = j*dim_chunk[0];
        f = (front_r == NULL) ? 0 : front_r[u];
        b = (dim_chunk[2] > 1) ? T[u + dim_chunk[0]*dim_chunk[1]] :
            ((back_r == NULL) ? 0 : back_r[u]);
        w = (west_r == NULL)  ? 0 : west_r[j];
        e = (dim_chunk[0] > 1) ? T[u + 1] :
            ((east_r == NULL) ? 0 : east_r[j]);
        
        R[u] = update_temperature(T[u], f, b, T[u - dim_chunk[0]],
                                  T[u + dim_chunk[0]], w, e, n, m, o, T_i, T_j + j, T_k);

        //i = dim_chunk[0]-1
        u = j*dim_chunk[0] + dim_chunk[0] - 1;
        f = (front_r == NULL) ? 0 : front_r[u];        
        b = (dim_chunk[2] > 1) ? T[u + dim_chunk[0]*dim_chunk[1]] :
            ((back_r == NULL) ? 0 : back_r[u]);
        w = (dim_chunk[0] > 1) ? T[u - 1] :
            ((west_r == NULL) ? 0 : west_r[j]);
        e = (east_r == NULL) ? 0 : east_r[j];

        R[u] = update_temperature(T[u], f, b, T[u - dim_chunk[0]], T[u + dim_chunk[0]], 
                                    w, e, n, m, o, T_i + dim_chunk[0] - 1, T_j + j, T_k);
    
    }

    // cas coins : 
    //coin bas gauche
    u = 0;
    f = (front_r == NULL) ? 0 : front_r[u];
    b = (dim_chunk[2] > 1) ? T[u + dim_chunk[0]*dim_chunk[1]] :
        ((back_r == NULL) ? 0 : back_r[u]);
    s = (south_r == NULL) ? 0 : south_r[0];
    no = (dim_chunk[1] > 1) ? T[u + dim_chunk[0]] :
        ((north_r == NULL) ? 0 : north_r[0]);
    w = (west_r == NULL)  ? 0 : west_r[0];
    e = (dim_chunk[0] > 1) ? T[u + 1] :
        ((east_r == NULL) ? 0 : east_r[0]);

    R[u] = update_temperature(T[u], f, b, s, no, w, e, n, m, o, T_i, T_j, T_k);
    // fprintf(stderr, "my_rank %d i %d j  %d k %d c %f f %f b %f s %f n %f w %f e %f\n", r, T_i, T_j, T_k, T[u], f, b, s, T[u + dim_chunk[0]], w, T[u + 1]);
    //coin bas droite 
    u = dim_chunk[0] - 1;
    f = (front_r == NULL) ? 0 : front_r[u];
    b = (dim_chunk[2] > 1) ? T[u + dim_chunk[0]*dim_chunk[1]] :
        ((back_r == NULL) ? 0 : back_r[u]);
    s = (south_r == NULL) ? 0 : south_r[u];
    no = (dim_chunk[1] > 1) ? T[u + dim_chunk[0]] : 
        ((north_r == NULL) ? 0 : north_r[u]);
    w = (dim_chunk[0] > 1) ? T[u -1] : 
        ((west_r == NULL) ? 0 : west_r[0]);
    e = (east_r == NULL)  ? 0 : east_r[0];

    R[u] = update_temperature(T[u], f, b, s, no, w, e, n, m, o,
                              T_i + dim_chunk[0] - 1, T_j, T_k);

    //coin haut gauche
    u = dim_chunk[0]*(dim_chunk[1] - 1);
    f = (front_r == NULL) ? 0 : front_r[u];
    b = (dim_chunk[2] > 1) ? T[u + dim_chunk[0]*dim_chunk[1]] :
        ((back_r == NULL) ? 0 : back_r[u]);
    s = (dim_chunk[1] > 1) ? T[u - dim_chunk[0]] :
        (south_r == NULL) ? 0 : south_r[0];
    no = (north_r == NULL) ? 0 : north_r[0];
    w = (west_r == NULL)  ? 0 : west_r[dim_chunk[1] - 1];
    e = (dim_chunk[0] > 1) ? T[u + 1] :
        ((east_r == NULL) ? 0 : east_r[dim_chunk[1] - 1]);
    
    R[u] = update_temperature(T[u], f, b, s, no, w, e, n, m, o, 
                              T_i, T_j + dim_chunk[1] - 1, T_k);
                                
    // coin haut droite
    u = dim_chunk[0]*dim_chunk[1] - 1;
    f  = (front_r == NULL) ? 0 : front_r[u];
    b = (dim_chunk[2] > 1) ? T[u + dim_chunk[0]*dim_chunk[1]] :
        ((back_r == NULL) ? 0 : back_r[u]);
    s = (dim_chunk[1] > 1) ? T[u - dim_chunk[0]] :
        ((south_r == NULL) ? 0 : south_r[dim_chunk[0] - 1]);
    no = (north_r == NULL) ? 0 : north_r[dim_chunk[0] - 1];
    w = (dim_chunk[0] > 1) ?  T[u - 1] :
        ((west_r == NULL) ? 0 : west_r[dim_chunk[1] - 1]);
    e = (east_r == NULL)  ? 0 : east_r[dim_chunk[1] - 1];
    
    R[u] = update_temperature(T[u], f, b, s,
                              no, w, e, n, m, o, T_i + dim_chunk[0] - 1, T_j + dim_chunk[1] - 1, T_k);

    end_front :

    if (dim_chunk[2] == 1) goto end_back;
    // FACE BACK
    // 3 cas : centre / arete / coins

    // cas centre :
    for (int j = 1; j < dim_chunk[1] - 1; ++j) {
        for (int i = 1; i < dim_chunk[0] - 1; ++i) {
            u = (dim_chunk[2]-1)*dim_chunk[1]*dim_chunk[0] + j*dim_chunk[0] + i;
            b = (back_r == NULL) ? 0 : back_r[j*dim_chunk[0] + i];
            
            R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], b, T[u - dim_chunk[0]],
                                      T[u + dim_chunk[0]], T[u - 1], T[u + 1], n, m, o, T_i + i, T_j + j, T_k + dim_chunk[2]-1);
        }
    }

    // cas aretes :
    // arete bas et haut
    for (int i = 1; i < dim_chunk[0]-1; ++i) {
        //j = 0
        u = (dim_chunk[2]-1)*dim_chunk[1]*dim_chunk[0] + i;
        b = (back_r == NULL)  ? 0 : back_r[i];
        s = (south_r == NULL) ? 0 : south_r[(dim_chunk[2]-1)*dim_chunk[0] + i];
        no = (dim_chunk[1] > 1) ? T[u + dim_chunk[0]] :
            ((north_r == NULL) ? 0 : north_r[(dim_chunk[2]-1)*dim_chunk[0] + i]);
        
        R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], b, s,
                                  no, T[u - 1], T[u + 1], n, m, o, T_i + i, T_j, T_k + dim_chunk[2]-1);
                                       
        //j = dim_chunk[1]-1
        u = (dim_chunk[2]-1)*dim_chunk[1]*dim_chunk[0] + dim_chunk[0]*(dim_chunk[1]-1) + i;
        b = (back_r == NULL)  ? 0 : back_r[dim_chunk[0]*(dim_chunk[1]-1) + i];
        s = (dim_chunk[1] > 1) ? T[u - dim_chunk[0]] :
            ((south_r == NULL) ? 0 : south_r[(dim_chunk[2]-1)*dim_chunk[0] + i]);
        no = (north_r == NULL) ? 0 : north_r[(dim_chunk[2]-1)*dim_chunk[0] + i];
        
        R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], b, s,
                                    no, T[u - 1], T[u + 1], n, m, o, T_i + i, T_j + dim_chunk[1] - 1, T_k + dim_chunk[2]-1);
    }
    
    //arete gauche et droite
    for (int j = 1; j < dim_chunk[1]-1; ++j) {
        //i = 0
        u = (dim_chunk[2]-1)*dim_chunk[1]*dim_chunk[0] + j*dim_chunk[0];
        b = (back_r == NULL) ? 0 : back_r[j*dim_chunk[0]];
        w = (west_r == NULL) ? 0 : west_r[(dim_chunk[2]-1)*dim_chunk[1] + j];
        e = (dim_chunk[0] > 1) ? T[u + 1] :
            ((east_r == NULL) ? 0 : east_r[(dim_chunk[2]-1)*dim_chunk[1] + j]);
        
        R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], b, T[u - dim_chunk[0]],
                                    T[u + dim_chunk[0]], w, e, n, m, o, T_i, T_j + j, T_k + dim_chunk[2]-1);

        //i = dim_chunk[0]-1
        u = (dim_chunk[2]-1)*dim_chunk[1]*dim_chunk[0] + j*dim_chunk[0] + dim_chunk[0] - 1;
        b = (back_r == NULL) ? 0 : back_r[j*dim_chunk[0] + dim_chunk[0] - 1];
        w = (dim_chunk[0] > 1) ?  T[u - 1]: 
            ((west_r == NULL) ? 0 : west_r[(dim_chunk[2]-1)*dim_chunk[1] + j]);
        e = (east_r == NULL) ? 0 : east_r[(dim_chunk[2]-1)*dim_chunk[1] + j];

        R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], b, T[u - dim_chunk[0]],
                                    T[u + dim_chunk[0]], w, e, n, m, o, T_i + dim_chunk[0] - 1, T_j + j, T_k + dim_chunk[2]-1);
    
    }

    // cas coins : 
    //coin bas gauche
    u = (dim_chunk[2]-1)*dim_chunk[1]*dim_chunk[0];
    b = (back_r == NULL)  ? 0 : back_r[0];
    s = (south_r == NULL) ? 0 : south_r[(dim_chunk[2]-1)*dim_chunk[0]];
    no = (dim_chunk[1] > 1) ? T[u + dim_chunk[0]] :
         ((north_r == NULL) ? 0 : north_r[(dim_chunk[2]-1)*dim_chunk[0]]);
    w = (west_r == NULL)  ? 0 : west_r[(dim_chunk[2]-1)*dim_chunk[1]];
    e = (dim_chunk[0] > 1) ? T[u + 1] :
        ((east_r == NULL) ? 0 : east_r[(dim_chunk[2]-1)*dim_chunk[1]]);

    R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], b, s,
                                no, w, e, n, m, o, T_i, T_j, T_k + dim_chunk[2]-1);

    //coin bas droite
    u = (dim_chunk[2]-1)*dim_chunk[1]*dim_chunk[0] + dim_chunk[0] - 1;
    b = (back_r == NULL)  ? 0 : back_r[dim_chunk[0] - 1];
    s = (south_r == NULL) ? 0 : south_r[dim_chunk[2]*dim_chunk[0] - 1];
    no = (dim_chunk[1] > 1) ? T[u + dim_chunk[0]] :
        ((north_r == NULL) ? 0 : north_r[dim_chunk[2]*dim_chunk[0] - 1]);
    w = (dim_chunk[0] > 1) ? T[u - 1] : 
        ((west_r == NULL) ? 0 : west_r[(dim_chunk[2]-1)*dim_chunk[1]]);
    e = (east_r == NULL)  ? 0 : east_r[(dim_chunk[2]-1)*dim_chunk[1]];

    R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], b, s,
                                no, w, e, n, m, o, T_i + dim_chunk[0] - 1, T_j, T_k + dim_chunk[2]-1);

    //coin haut gauche
    u = (dim_chunk[2]-1)*dim_chunk[1]*dim_chunk[0] + dim_chunk[0]*(dim_chunk[1] - 1);
    b  = (back_r == NULL)  ? 0 : back_r[(dim_chunk[1] - 1)*dim_chunk[0]];
    s = (dim_chunk[1] >  1) ? T[u - dim_chunk[0]] :
        ((south_r == NULL) ? 0 : south_r[(dim_chunk[2]-1)*dim_chunk[0]]);
    no = (north_r == NULL) ? 0 : north_r[(dim_chunk[2]-1)*dim_chunk[0]];
    w = (west_r == NULL)  ? 0 : west_r[dim_chunk[1]*dim_chunk[2] - 1];
    e = (dim_chunk[0] > 1) ? T[u + 1]:
        ((east_r == NULL) ? 0 : east_r[dim_chunk[1]*dim_chunk[2] - 1]);

    R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], b, s,
                                no, w, e, n, m, o, T_i, T_j + dim_chunk[1] - 1, T_k + dim_chunk[2]-1);
                                
    // coin haut droite
    u = dim_chunk[2]*dim_chunk[1]*dim_chunk[0] - 1;
    b = (back_r == NULL)  ? 0 : back_r[dim_chunk[1]*dim_chunk[0] - 1];
    s = (dim_chunk[1] >  1) ? T[u - dim_chunk[0]] :
        ((south_r == NULL) ? 0 : south_r[dim_chunk[0]*dim_chunk[2] - 1]);
    no = (north_r == NULL) ? 0 : north_r[dim_chunk[0]*dim_chunk[2] - 1];
    w = (dim_chunk[0] > 1) ? T[u - 1] :
        ((west_r == NULL) ? 0 : west_r[dim_chunk[1]*dim_chunk[2] - 1]);
    e = (east_r == NULL)  ? 0 : east_r[dim_chunk[1]*dim_chunk[2] - 1];
    
    R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], b, s,
                                no, w, e, n, m, o, T_i + dim_chunk[0] - 1, T_j + dim_chunk[1] - 1, T_k + dim_chunk[2]-1);

    end_back :

    // FACE SOUTH
    // 2 cas : centre / arete

    // cas centre :
    for (int k = 1; k < dim_chunk[2] - 1; ++k) {
        for (int i = 1; i < dim_chunk[0] - 1; ++i) {
            u = k*dim_chunk[1]*dim_chunk[0] + i;
            s = (south_r == NULL) ? 0 : south_r[k*dim_chunk[0] + i];
            no = (dim_chunk[1] > 1) ? T[u + dim_chunk[0]] : 
                 ((north_r == NULL) ? 0 : north_r[k*dim_chunk[0] + i]);
            
            R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], T[u + dim_chunk[0]*dim_chunk[1]], s,
                                      no, T[u - 1], T[u + 1], n, m, o, T_i + i, T_j, T_k + k);
        }
    }
    
    //arete bas-gauche et bas-droite
    for (int k = 1; k < dim_chunk[2] - 1; ++k) {
        //i = 0
        u = k*dim_chunk[1]*dim_chunk[0];
        s = (south_r == NULL) ? 0 : south_r[k*dim_chunk[0]];
        no = (dim_chunk[1] > 1) ? T[u + dim_chunk[0]] : 
             ((north_r == NULL) ? 0 : north_r[k*dim_chunk[0]]);
        w = (west_r == NULL)  ? 0 : west_r[k*dim_chunk[1] + dim_chunk[1] - 1];
        e = (dim_chunk[0] > 1) ? T[u + 1] :
            ((east_r == NULL) ? 0 : east_r[k*dim_chunk[1] + dim_chunk[1] - 1]);
        
        R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], T[u + dim_chunk[0]*dim_chunk[1]], s,
                                  no, w, e, n, m, o, T_i, T_j, T_k + k);

        //i = dim_chunk[0]-1
        u = k*dim_chunk[1]*dim_chunk[0] + dim_chunk[0] - 1;
        s = (south_r == NULL) ? 0 : south_r[k*dim_chunk[0] + dim_chunk[0] - 1];
        no = (dim_chunk[1] > 1) ? T[u + dim_chunk[0]] : 
             ((north_r == NULL) ? 0 : north_r[k*dim_chunk[0] + dim_chunk[0] - 1]);
        w = (dim_chunk[0] > 1) ? T[u - 1] :
            ((west_r == NULL) ? 0 : west_r[k*dim_chunk[1] + dim_chunk[1] - 1]);
        e = (east_r == NULL)  ? 0 : east_r[k*dim_chunk[1] + dim_chunk[1] - 1];

        //fprintf(stderr, "no_s = %lf no = %lf\n", no, north_r[k*dim_chunk[0] + dim_chunk[0] - 1]);
        R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], T[u + dim_chunk[0]*dim_chunk[1]], s,
                                   no, T[u - 1], e, n, m, o, T_i + dim_chunk[0] - 1, T_j, T_k + k);
    
    }

    

    // FACE NORTH
    // 2 cas : centre / arete
    if (dim_chunk[1] == 1) goto end_north;

    // cas centre :
    for (int k = 1; k < dim_chunk[2] - 1; ++k) {
        for (int i = 1; i < dim_chunk[0] - 1; ++i) {
            u = k*dim_chunk[1]*dim_chunk[0] + dim_chunk[0]*(dim_chunk[1]-1) + i;
            no = (north_r == NULL) ? 0 : north_r[k*dim_chunk[0] + i];
            
            R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], T[u + dim_chunk[0]*dim_chunk[1]], T[u - dim_chunk[0]],
                                      no, T[u - 1], T[u + 1], n, m, o, T_i + i, T_j + dim_chunk[1] - 1, T_k + k);
        }
    }
    
    //arete haut-gauche et haut-droite
    for (int k = 1; k < dim_chunk[2] - 1; ++k) {
        //i = 0
        u = k*dim_chunk[1]*dim_chunk[0] + (dim_chunk[1]-1)*dim_chunk[0];
        no = (north_r == NULL) ? 0 : north_r[k*dim_chunk[0]];
        w = (west_r == NULL)  ? 0 : west_r[k*dim_chunk[1] + dim_chunk[1] - 1];
        e = (dim_chunk[0] > 1) ? T[u + 1] :
            ((east_r == NULL) ? 0 : east_r[k*dim_chunk[1] + dim_chunk[1] - 1]);
        
        R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], T[u + dim_chunk[0]*dim_chunk[1]], T[u - dim_chunk[0]],
                                  no, w, e, n, m, o, T_i, T_j + dim_chunk[1] - 1, T_k + k);

        //i = dim_chunk[0]-1
        u = k*dim_chunk[1]*dim_chunk[0] + (dim_chunk[1]-1)*dim_chunk[0] + dim_chunk[0] - 1;
        no = (north_r == NULL) ? 0 : north_r[k*dim_chunk[0] + dim_chunk[0] - 1];
        w = (dim_chunk[0] > 1) ? T[u - 1] :
            ((west_r == NULL) ? 0 : west_r[k*dim_chunk[1] + dim_chunk[1] - 1]);
        e =  (east_r == NULL) ? 0 : east_r[k*dim_chunk[1] + dim_chunk[1] - 1];

        R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], T[u + dim_chunk[0]*dim_chunk[1]], T[u - dim_chunk[0]],
                                  no, w, e, n, m, o, T_i + dim_chunk[0] - 1, T_j + dim_chunk[1] - 1, T_k + k);
    
    }

    end_north :

    // FACE WEST
    // 2 cas : centre / arete

    // cas centre :
    for (int k = 1; k < dim_chunk[2] - 1; ++k) {
        for (int j = 1; j < dim_chunk[1] - 1; ++j) {
            u = k*dim_chunk[1]*dim_chunk[0] + j*dim_chunk[0];
            w = (west_r == NULL) ? 0 : west_r[k*dim_chunk[1] + j];
            e = (dim_chunk[0] > 1) ? T[u + 1] : 
                ((east_r == NULL) ? 0 : east_r[k*dim_chunk[1] + j]);
            
            R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], T[u + dim_chunk[0]*dim_chunk[1]], T[u - dim_chunk[0]],
                                      T[u + dim_chunk[0]], w, e, n, m, o, T_i, T_j + j, T_k + k);
        }
    }

    if (dim_chunk[0] == 1) return;
    // FACE EAST
    // 2 cas : centre / arete

    // cas centre :
    for (int k = 1; k < dim_chunk[2] - 1; ++k) {
        for (int j = 1; j < dim_chunk[1] - 1; ++j) {
            u = k*dim_chunk[1]*dim_chunk[0] + j*dim_chunk[0] + dim_chunk[0] - 1;
            e = (east_r == NULL) ? 0 : east_r[k*dim_chunk[1] + j];
            
            R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], T[u + dim_chunk[0]*dim_chunk[1]], T[u - dim_chunk[0]],
                                      T[u + dim_chunk[0]], T[u - 1], e, n, m, o, T_i + dim_chunk[0] - 1, T_j + j, T_k + k);
        }
    }
}

int *get_3_dividers(int n)
{
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

    if (k<j) {
        if (i<j) {
            tab[2] = j;
            if (k<i) {
                tab[1] = i;
                tab[0] = k;
            }
            else {
                tab[1] = k;
                tab[0] = i;
            }
        }
        else {
            tab[2] = i;
            tab[1] = j;
            tab[0] = k;
        }
        
    }

    else {
        if (i<k) {
            tab[2] = k;
            if (j<i) {
                tab[1] = i;
                tab[0] = j;
            }
            else {
                tab[1] = j;
                tab[0] = i;
            }
        }
        else {
            tab[2] = i;
            tab[1] = k;
            tab[0] = j;
        }
        
    }

    return tab;
}
    
 

// void fill_halo(double *tab, int *dim_chunk, double *halo, int index) { // A MODIFIER
//     /* index = -1 (center) / 0 (front) / 1 (back) / 2 (west) / 3 (east) / 4 (south) /5 (north) */
    
//     switch(index) {
        
//         case -1 : // center
//             for (int i = 1; i  < dim_chunk[0]+1; ++i) {

//                 for (int j = 1; j < dim_chunk[1]+1; ++j) {

//                     for (int k = 1; k < dim_chunk[2]+1; ++k) {
//                         halo[i + j*(dim_chunk[0]+2) + k*(dim_chunk[0]+2)*(dim_chunk[1]+2)] 
//                         = tab[(i-1) + (j-1)*dim_chunk[0] + (k-1)*dim_chunk[0]*dim_chunk[1]];
//                     }
//                 }
//             }
//             break;
        
//         case 0 : // east
//             for (int j = 1; j < dim_chunk[1]+1; ++j) {

//                 for (int k = 1; k < dim_chunk[2]+1; ++k) { // i = dim_chunk[0]-1
//                     halo[dim_chunk[0]-1 + j*(dim_chunk[0]+2) + k*(dim_chunk[0]+2)*(dim_chunk[1]+2)]
//                         = tab[j*dim_chunk[2] + k];
//                 }
//             }
//             break;

//         case 1 : // south
//             for (int i = 1; i < dim_chunk[0]+1; ++i) {

//                 for (int k = 1 ; k < dim_chunk[2]+1; ++k) { // j = 0
//                     halo[i + k*(dim_chunk[0]+2)*(dim_chunk[1]+2)] 
//                         = tab[i*dim_chunk[2] + k];
//                 }

//             }
//             break;

//         case 2 : // west
//             for (int j = 1; j < dim_chunk[1]+1; ++j) {

//                 for (int k = 1; k < dim_chunk[2]+1; ++k) { // i = 0
//                     halo[j*(dim_chunk[0]+2) + k*(dim_chunk[0]+2)*(dim_chunk[1]+2)]
//                         = tab[j*dim_chunk[2] + k];
//                 }
//             }
//             break;
        
//         case 3 : // north
//             for (int i = 1; i < dim_chunk[0]+1; ++i) {

//                 for (int k = 1 ; k < dim_chunk[2]+1; ++k) { // j = dim_chunk[1]-1
//                     halo[i + (dim_chunk[1]-1)*dim_chunk[0] + k*(dim_chunk[0]+2)*(dim_chunk[1]+2)] 
//                         = tab[i*dim_chunk[2] + k];
//                 }

//             }
//             break;
            
//         default :
//             return;
//     }
// }

int main(int argc, char **argv)
{
   	
    int my_rank; /* rank of the process */
    int p; /* number of processes */
    //int dest; /*rank of the dest */
    // int source; /* rank of the source */
    // int tag = 0; /* tag of the message */
    //char message[100];
    //MPI_Status status;
    //char hostname[SIZE_H_N] ;
    //gethostname(hostname, SIZE_H_N);	
       
    /* Initialisation */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    r = my_rank; // pour le test

    CPU_surface = CPU_contact_surface();
    double V = L * l * E;
    int n = ceil(L / dl);
    int m = ceil(E / dl);
    int o = ceil(l / dl);
    int *dividers = get_3_dividers(p);
    int dim = 3;
    for (int i = 0; i < 3; i++){
    	if (dividers[i] == 1) dim--;
    }
    int *dim_chunk = malloc(3*sizeof(int));
    if (dim_chunk == NULL){
    	perror("malloc");
	exit(1);
    }
    dim_chunk[0] = n/dividers[0];
    dim_chunk[1] = m/dividers[1];
    dim_chunk[2] = o/dividers[2];
    fprintf(stderr, "n %d, dim_chunk[0] %d\n", n, dim_chunk[0]);
    fprintf(stderr, "m %d, dim_chunk[1] %d\n", m, dim_chunk[1]);
    fprintf(stderr, "o %d, dim_chunk[2] %d\n", o, dim_chunk[2]);

    fprintf(stderr, "dimchunk : %d %d %d\n",dim_chunk[0], dim_chunk[1], dim_chunk[2]);
    
    fprintf(stderr, "HEATSINK\n");
    fprintf(stderr, "\tDimension (cm) [x,y,z] = %.1f x %.1f x %.1f\n", 100 * L, 100 * E, 100 * l);
    fprintf(stderr, "\tVolume = %.1f cm^3\n", V * 1e6);
    fprintf(stderr, "\tWeight = %.2f kg\n", V * sink_density);
    fprintf(stderr, "\tPrice = %.2f €\n", V * sink_density * euros_per_kg);
    fprintf(stderr, "SIMULATION\n");
    fprintf(stderr, "\tGrid (x,y,z) = %d x %d x %d (%.1fMo)\n", n, m, o, 7.6293e-06 * n * m * o);
    fprintf(stderr, "\tdt = %gs\n", dt);
    fprintf(stderr, "CPU\n");
    fprintf(stderr, "\tPower = %.0fW\n", CPU_TDP);
    fprintf(stderr, "\tArea = %.1f cm^2\n", CPU_surface * 10000);

    //temperature of each cell, in degree Kelvin
    /*
    double *T = malloc(n * m * o * sizeof(*T));
    double *R = malloc(n * m * o * sizeof(*R));
    if (T == NULL || R == NULL) {
        perror("T or R could not be allocated");
        exit(1);
    }
    */

    int size_chunk = dim_chunk[0]*dim_chunk[1]*dim_chunk[2];
    double *chunk = NULL;   // chunk(i,j,k) = chunk[i + j*chunk0 + k*chunk0*chunk1]
    double *new_chunk = NULL;

    double *front_r = NULL; // id = 0 // front(i,j) = front[j*chunk0 + i]
    double *back_r = NULL;  // id = 1 // back(i,j)  = back[j*chunk0 + i]
    double *south_r = NULL; // id = 2 // south(i,k) = south[k*chunk0 + i] = chunk(i,-1,k)
    double *north_r = NULL; // id = 3 // north(i,k) = north[k*chunk0 + i] = chunk(i,dim_chunk1,k)
    double *west_r = NULL;  // id = 4 // west(j,k)  = west[k*chunk1 + j]  = chunk(-1,j,k)
    double *east_r = NULL;  // id = 5 // east(j,k)  = east[k*chunk1 + j]  = chunk(dim_chunk0,j,k)

    double *front_s = NULL; // front(i,j) = front[j*chunk0 + i]
    double *back_s = NULL;  // back(i,j)  = back[j*chunk0 + i]
    double *south_s = NULL; // south(i,k) = south[k*chunk0 + i] = chunk(i,0,k)
    double *north_s = NULL; // north(i,k) = north[k*chunk0 + i] = chunk(i,0,k)
    double *west_s = NULL;  // west(j,k)  = west[k*chunk1 + j]  = chunk(0,j,k)    
    double *east_s = NULL;  // east(j,k)  = east[k*chunk1 + j]  = chunk(0,j,k)

    chunk = malloc(size_chunk*sizeof(double));
    if (chunk == NULL){
        perror("malloc_blocs : malloc");
	    exit(1);
    }
    new_chunk = malloc(size_chunk*sizeof(double));
    if (new_chunk == NULL){
        perror("malloc_blocs : malloc");
	    exit(1);
    }
    
    if (dividers[2] > 1) {
        front_r = malloc((dim_chunk[0])*(dim_chunk[1])*sizeof(double));
        if (front_r == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        back_r = malloc((dim_chunk[0])*(dim_chunk[1])*sizeof(double));
        if (back_r == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        front_s = malloc((dim_chunk[0])*(dim_chunk[1])*sizeof(double));
        if (front_s == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        back_s = malloc((dim_chunk[0])*(dim_chunk[1])*sizeof(double));
        if (back_s == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }

    }

    if (dividers[1] > 1) {
        south_r = malloc((dim_chunk[0])*(dim_chunk[2])*sizeof(double));
        if (south_r == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        north_r = malloc((dim_chunk[0])*(dim_chunk[2])*sizeof(double));
        if (north_r == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        south_s = malloc((dim_chunk[0])*(dim_chunk[2])*sizeof(double));
        if (south_s == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        north_s = malloc((dim_chunk[0])*(dim_chunk[2])*sizeof(double));
        if (north_s == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
    }

    if (dividers[0] > 1) {
        west_r = malloc((dim_chunk[1])*(dim_chunk[2])*sizeof(double));
        if (west_r == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        east_r = malloc((dim_chunk[1])*(dim_chunk[2])*sizeof(double));
        if (east_r == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        west_s = malloc((dim_chunk[1])*(dim_chunk[2])*sizeof(double));
        if (west_s == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        east_s = malloc((dim_chunk[1])*(dim_chunk[2])*sizeof(double));
        if (east_s == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
    }


    /* initially the heatsink is at the temperature of the water-cooling fluid */
    for (int u = 0; u < size_chunk; u++)
       new_chunk[u] = chunk[u] = watercooling_T + 273.15;

    /* let's go! we switch the CPU on and launch the simulation until it reaches a stationary state. */
    double t = 0;
    int n_steps = 0;
    int convergence = 0;


    /* simulating time steps */
    while (convergence == 0) {
        it = n_steps;
        /* Update all cells. xy planes are processed, for increasing values of z. */
	/* Ask for all the Irecv (check the border conditions)
	 * prepare the data and Isend them
	 * Start the calculations for the inside of the block
	 * wait for the Irecv and calculate each as soon as it's sent
	 * maybe wait for the Isend before going to t+1 because the values Isent must not change before they're effectively sent
	 * */

        /* receive the sides from others */
        // int test_f = my_rank >= dividers[2]*dividers[1];
        // int test_b = my_rank < dividers[2]*dividers[1]*(dividers[0]-1);



        // int test_s = my_rank%(dividers[2]*dividers[1]) >= dividers[2]*(dividers[1]-1);
        // int test_n = my_rank%(dividers[2]*dividers[1]) < dividers[2];



        // int test_w = my_rank%dividers[2] != 0;
        // int test_e = my_rank%dividers[2] != dividers[2]-1;

        int test_f = my_rank >= dividers[0]*dividers[1];
        int test_b = my_rank < dividers[0]*dividers[1]*(dividers[2]-1);
        int test_s = my_rank%(dividers[0]*dividers[1])  - dividers[0] >= 0 ;
        int test_n = my_rank%(dividers[0]*dividers[1])  + dividers[0] < dividers[0]*dividers[1];
        int test_w = my_rank%dividers[0] != 0;
        int test_e = my_rank%dividers[0] != dividers[0]-1;

        int nb_neigh = test_f + test_b + test_s + test_n + test_w + test_e;
        int cpt = 0;
        
        MPI_Request *rec = calloc(nb_neigh, sizeof(MPI_Request));
        if (test_f) {
            // printf("0\n");
            MPI_Irecv(front_r, (dim_chunk[0])*(dim_chunk[1]),
                    MPI_DOUBLE, my_rank-dividers[0]*dividers[1], 0, MPI_COMM_WORLD,
                    &rec[cpt++]);
        }
        if (test_b) {
            // printf("1\n");
            MPI_Irecv(back_r, (dim_chunk[0])*(dim_chunk[1]),
                    MPI_DOUBLE, my_rank+dividers[0]*dividers[1], 0, MPI_COMM_WORLD,
                    &rec[cpt++]);
        }
        if (test_s) {
            // printf("2\n");
            MPI_Irecv(south_r, (dim_chunk[0])*(dim_chunk[2]),
                    MPI_DOUBLE, my_rank-dividers[0], 0, MPI_COMM_WORLD,
                    &rec[cpt++]);
        } 
        if (test_n) {
            // printf("3\n");
            MPI_Irecv(north_r, (dim_chunk[0])*(dim_chunk[2]),
                    MPI_DOUBLE, my_rank+dividers[0], 0, MPI_COMM_WORLD,
                    &rec[cpt++]);
        }
        if (test_w) {
            // printf("4\n");
            MPI_Irecv(west_r, (dim_chunk[1])*(dim_chunk[2]),
                    MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD,
                    &rec[cpt++]);
        }
        if (test_e) {
            // printf("5\n");
            MPI_Irecv(east_r, (dim_chunk[1])*(dim_chunk[2]),
                    MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD,
                    &rec[cpt]);
        }


        /* create the sides for neighbours */
        if (test_f && test_b) {
            for (int i = 0; i < dim_chunk[0]; ++i) {
                for (int j = 0; j < dim_chunk[1]; ++j) { /* (x,y) pour les faces front back*/
                    front_s[j*dim_chunk[0] + i] = chunk[j*dim_chunk[0] + i];
                    back_s[j*dim_chunk[0] + i]  = chunk[j*dim_chunk[0] + i + dim_chunk[0]*dim_chunk[1]*(dim_chunk[2]-1)];
                }
            }
        }
        else if (test_f) {
            for (int i = 0; i < dim_chunk[0]; ++i) {
                for (int j = 0; j < dim_chunk[1]; ++j) { /* (x,y) pour les faces front back*/
                    front_s[j*dim_chunk[0] + i] = chunk[j*dim_chunk[0] + i];
                }
            }
        }
        else if (test_b) {
            for (int i = 0; i < dim_chunk[0]; ++i) {
                for (int j = 0; j < dim_chunk[1]; ++j) { /* (x,y) pour les faces front back*/
                    back_s[j*dim_chunk[0] + i] = chunk[j*dim_chunk[0] + i + dim_chunk[0]*dim_chunk[1]*(dim_chunk[2]-1)];
                }
            }
        } 

        if (test_n && test_s) {
            for (int i = 0; i < dim_chunk[0]; ++i) {
                for (int k = 0; k < dim_chunk[2]; ++k) {
                    south_s[k*dim_chunk[0] + i] = chunk[k*dim_chunk[0]*dim_chunk[1] + i];
                    north_s[k*dim_chunk[0] + i] = chunk[k*dim_chunk[0]*dim_chunk[1] + i + (dim_chunk[1]-1)*dim_chunk[0]];
                }
            }
        }
        else if (test_n) {
            for (int i = 0; i < dim_chunk[0]; ++i) {
                for (int k = 0; k < dim_chunk[2]; ++k) {
                    north_s[k*dim_chunk[0] + i] = chunk[k*dim_chunk[0]*dim_chunk[1] + i + (dim_chunk[1]-1)*dim_chunk[0]];
                }
            }
        }
        else if (test_s) {
            for (int i = 0; i < dim_chunk[0]; ++i) {
                for (int k = 0; k < dim_chunk[2]; ++k) {
                    south_s[k*dim_chunk[0] + i] = chunk[k*dim_chunk[0]*dim_chunk[1] + i];
                }
            }
        }
        
        if (test_w && test_e) {
            for (int j = 0; j < dim_chunk[1]; ++j) { /* (y,z) pour les faces west east */
                for (int k = 0; k < dim_chunk[2]; ++k) {
                    west_s[k*dim_chunk[1] + j] = chunk[k*dim_chunk[0]*dim_chunk[1] + j*dim_chunk[0]];
                    east_s[k*dim_chunk[1] + j] = chunk[k*dim_chunk[0]*dim_chunk[1] + j*dim_chunk[0] + dim_chunk[0]-1];
                }
            } 
        }
        else if (test_w) {
            for (int j = 0; j < dim_chunk[1]; ++j) { /* (y,z) pour les faces west east */
                for (int k = 0; k < dim_chunk[2]; ++k) {
                    west_s[k*dim_chunk[1] + j] = chunk[k*dim_chunk[0]*dim_chunk[1] + j*dim_chunk[0]];
                }
            } 
        }
        else if (test_e) {
            for (int j = 0; j < dim_chunk[1]; ++j) { /* (y,z) pour les faces west east */
                for (int k = 0; k < dim_chunk[2]; ++k) {
                    east_s[k*dim_chunk[1] + j] = chunk[k*dim_chunk[0]*dim_chunk[1] + j*dim_chunk[0] + dim_chunk[0]-1];
                }
            } 
        }

        
        cpt = 0;
        /* sends the sides to neighbours */
        MPI_Request *sen = calloc(nb_neigh, sizeof(MPI_Request));
        if (test_f) {
            MPI_Isend(front_s, (dim_chunk[0])*(dim_chunk[1]),
                    MPI_DOUBLE, my_rank-dividers[0]*dividers[1], 0, MPI_COMM_WORLD,
                    &sen[cpt++]);
        }
        if (test_b) {
            MPI_Isend(back_s, (dim_chunk[0])*(dim_chunk[1]),
                    MPI_DOUBLE, my_rank+dividers[0]*dividers[1], 0, MPI_COMM_WORLD,
                    &sen[cpt++]);
        }
        if (test_s) {
            MPI_Isend(south_s, (dim_chunk[0])*(dim_chunk[2]),
                    MPI_DOUBLE, my_rank-dividers[0], 0, MPI_COMM_WORLD,
                    &sen[cpt++]);
        } 
        if (test_n) {
            MPI_Isend(north_s, (dim_chunk[0])*(dim_chunk[2]),
                    MPI_DOUBLE, my_rank+dividers[0], 0, MPI_COMM_WORLD,
                    &sen[cpt++]);
        }
        if (test_w) {
            MPI_Isend(west_s, (dim_chunk[1])*(dim_chunk[2]),
                    MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD,
                    &sen[cpt++]);
        }
        if (test_e) {
            MPI_Isend(east_s, (dim_chunk[1])*(dim_chunk[2]),
                    MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD,
                    &sen[cpt]);
        }
        

        // int *edges = calloc(sizeof(int), 4);
        // if (edges == NULL) {
        //     perror("malloc main");
        //     return EXIT_FAILURE;
        // }
        
        /* compute the inside of the chunk */
        for (int k = 1; k < dim_chunk[0]-1; ++k) {
            do_xy_plane(chunk, new_chunk, n, m, o, k, dim_chunk, my_rank, 1);
        }

        // MPI_Barrier(MPI_COMM_WORLD);
        int indx;

        for (int i = 0; i < nb_neigh; ++i) {
            MPI_Waitany(nb_neigh, rec, &indx, MPI_STATUSES_IGNORE);
            // fprintf(stderr, "test13 1/2\n");


            // edges[indx] = 1;

            // switch(indx) { // A REVOIR
            //     case 0 : //east
            //         fill_halo(east_r, dim_chunk, halo, 0);
            //         do_yz_plane(halo, new_chunk, edges, n, m, o, dim_chunk[0]-1, dim_chunk, p, 0);
            //         break;

            //     case 1 : //south
            //         fill_halo(east_r, dim_chunk, halo, 0);
            //         do_xz_plane(halo, new_chunk, edges, n, m, o, 0, dim_chunk, p, 0);
            //         break;

            //     case 2 : //west
            //         fill_halo(east_r, dim_chunk, halo, 0);
            //         do_yz_plane(halo, new_chunk, edges, n, m, o, 0, dim_chunk, p, 0);
            //         break;

            //     case 3 : //north
            //         fill_halo(east_r, dim_chunk, halo, 0);
            //         do_xz_plane(halo, new_chunk, edges, n, m, o, dim_chunk[1]-1, dim_chunk, p, 0);
            //         break;
                    
            //     default :
            //         return EXIT_FAILURE;
            // }
        }
        // if (test_n) {
        //     for (int i = 0; i < dim_chunk[0]; ++i) {
        //         for (int k = 0; k < dim_chunk[2]; ++k) {
        //             fprintf(stderr, "north i %d k %d | val = %f\n", i, k, north_r[k*dim_chunk[0] + i]-273.15);
        //         }
        //     }
        // }
        // if (test_s) {
        //     for (int i = 0; i < dim_chunk[0]; ++i) {
        //         for (int k = 0; k < dim_chunk[2]; ++k) {
        //             fprintf(stderr, "south i %d k %d | val = %f\n", i, k, south_r[k*dim_chunk[0] + i]-273.15);
        //         }
        //     }
        // }
        // fprintf(stderr, "test14\n");

        do_halo(chunk, new_chunk, front_r, back_r, south_r, north_r, west_r, east_r, n, m, o, dim_chunk, my_rank);
        // fprintf(stderr, "test15\n");

        /* each second, we test the convergence, and print a short progress report */
        //test it locally, and sum all the partial values with MPI_allReduceand MPI_SUM
        //Be careful : the sqrt has to be done on the total sum 

        //each processor compute the sum
        if (n_steps % ((int)(1 / dt)) == 0) {
            double delta_T = 0;
            double max = -INFINITY;

            for (int u = 0; u < dim_chunk[0]*dim_chunk[1]*dim_chunk[2]; u++) {
                // if (new_chunk[u] < chunk[u])
                // {fprintf(stderr,"new_c %f c %f \n", new_chunk[u], chunk[u]);}
                delta_T += (new_chunk[u] - chunk[u]) * (new_chunk[u] - chunk[u]);
                if (new_chunk[u] > max)
                    max = new_chunk[u];
            }

            /*send all sums to 0*/
            double sum = 0;
            double max_all = 0;
            MPI_Reduce(&delta_T, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&max, &max_all, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


            //only 0
            if (my_rank == 0) {

                delta_T = sqrt(sum) / dt;
                fprintf(stderr, "t = %.1fs ; T_max = %.1f°C ; convergence = %g\n", t, max_all - 273.15, delta_T);
                
                if (delta_T < 0.1) {
                    convergence = 1;
                }
                

            }
            MPI_Bcast(&convergence, 1, MPI_INT, 0, MPI_COMM_WORLD);
	    }

        /* the new temperatures are in R */
        double *tmp = new_chunk;
        new_chunk = chunk; 
        chunk = tmp;

        t += dt;
        n_steps += 1;

        free(rec);
        free(sen);
    }

    // fprintf(stderr, "fin while %d dim %d\n", my_rank, dim);


    if (dim == 3) {
        MPI_Comm row_comm;
        MPI_Comm_split(MPI_COMM_WORLD, my_rank/dividers[0], my_rank%dividers[0], &row_comm);
        double *row_gather = NULL;
        if (my_rank%dividers[0] == 0) {
            row_gather = malloc(dividers[0]*size_chunk*sizeof(double));
            if (row_gather == NULL) {
                perror("main : malloc row_gather");
                return EXIT_FAILURE;
            }
        }
	      // fprintf(stderr, "entry gather dim 3 my_rank %d \n", my_rank);
        MPI_Gather(chunk, size_chunk, MPI_DOUBLE, row_gather,
                   size_chunk, MPI_DOUBLE, 0, row_comm);
	      // fprintf(stderr, "exit gather dim 3 my_rank %d \n", my_rank);
        free(chunk);
        chunk = NULL;

        size_chunk = dividers[0]*size_chunk;

        if (my_rank%dividers[0] == 0) {
            chunk = malloc(size_chunk * sizeof(double));
            
            /* reorganisation of the chunk so it is in the write order */
            for (int k = 0; k < dim_chunk[2]; ++k) {
                for (int j = 0; j < dim_chunk[1]; ++j) {
                    for (int i = 0; i < dim_chunk[0]; ++i) {
                        for (int l = 0; l < dividers[0] ; ++l) {

                            chunk[ dividers[0]*k*dim_chunk[0]*dim_chunk[1] + 
                                    dividers[0]*j*dim_chunk[0] + l*dim_chunk[0] + i ] =
                                
                                    row_gather[ l*(dim_chunk[0]*dim_chunk[1]*dim_chunk[2]) +
                                                k*dim_chunk[0]*dim_chunk[1] + j*dim_chunk[0] + i ];
                        }
                    }
                }
            }
        }
        free(row_gather);

        dim_chunk[0] = dividers[0]*dim_chunk[0];
    }

    if (dim >= 2) {
        MPI_Comm col_comm;
        int color = (my_rank/(dividers[0]*dividers[1]))*dividers[0] + 
                    my_rank%dividers[0];
        int key = (my_rank%(dividers[0]*dividers[1]))/dividers[0];
        MPI_Comm_split(MPI_COMM_WORLD, color, key, &col_comm);
        double *col_gather = NULL;
        if (key == 0) {
            col_gather = malloc(dividers[1]*size_chunk*sizeof(double));
            if (col_gather == NULL) {
                perror("main : malloc col_gather");
                return EXIT_FAILURE;
            }
        }

        // fprintf(stderr, "entry gather dim 2 my_rank %d color %d key %d\n", my_rank, color, key);
        if (my_rank%dividers[0] == 0) {
            // fprintf(stderr, "inside gather dim 2 my_rank %d\n", my_rank);

            MPI_Gather(chunk, size_chunk, MPI_DOUBLE, col_gather,
                       size_chunk, MPI_DOUBLE, 0, col_comm);
        }
        // fprintf(stderr, "exit gather dim2 my_rank %d\n", my_rank);

        
        free(chunk);
        chunk = NULL;
        size_chunk = dividers[1]*size_chunk;

        if (key == 0) {
            // fprintf(stderr, "inside loop my_rank %d\n", my_rank);
            chunk = malloc(size_chunk * sizeof(double));

            for (int k = 0; k < dim_chunk[2]; ++k) {
                for (int j = 0; j < dim_chunk[1]; ++j) {
                    for (int i = 0; i < dim_chunk[0]; ++i) {
                        for (int l = 0; l < dividers[1] ; ++l) {
                            
                            chunk[ dividers[1]*k*dim_chunk[0]*dim_chunk[1] + 
                                    l*dim_chunk[0]*dim_chunk[1] + 
                                    dividers[1]*j*dim_chunk[0] + i] = 
                                
                                col_gather[l*(dim_chunk[0]*dim_chunk[1]*dim_chunk[2]) + 
                                    k*dim_chunk[0]*dim_chunk[1] + j*dim_chunk[0] + i];
                        }
                    }
                }
            }
        }
        free(col_gather);

        dim_chunk[1] = dim_chunk[1]*dividers[1];
    }
    // fprintf(stderr, "fin zonedim2 %d dim %d\n", my_rank, dim);

    double *T = NULL;
    if (my_rank == 0) {
        T = malloc(dividers[2]*size_chunk*sizeof(double));
        if (T == NULL) {
            perror("main: malloc T");
            return EXIT_FAILURE;
        }
    }
    
    MPI_Comm final_comm;
    MPI_Comm_split(MPI_COMM_WORLD, (my_rank%(dividers[0]*dividers[1])), 
            my_rank/(dividers[0]*dividers[1]), &final_comm);
        

     //fprintf(stderr, "n1 %d, n2 %d, n3 %d\n", dim_chunk[0]*p,dim_chunk[1]*p,dim_chunk[2]*p);
    // fprintf(stderr, "entry gather %d\n", my_rank);
    if (my_rank%(dividers[0]*dividers[1]) == 0) {
        //  fprintf(stderr, "inside gather my_rank %d\n", my_rank);
        MPI_Gather(chunk, size_chunk, MPI_DOUBLE, T,
                    size_chunk, MPI_DOUBLE, 0, 
                    final_comm);
    }
    //  fprintf(stderr, "exit gather %d\n", my_rank);


    if (my_rank == 0) {
        #ifdef DUMP_STEADY_STATE
        fprintf(stderr,"###### STEADY STATE; t = %.1f\n", t);
        for (int k = 0; k < o; ++k) {   // z
            fprintf(stderr, "# z = %g\n", k * dl);
            for (int j = 0; j < m; ++j) {   // y
                for (int i = 0; i < n; ++i) {   // x
                    
                    fprintf(stderr,"%.1f ", T[k * n * m + j * n + i] - 273.15);
                }
                fprintf(stderr,"\n");
            }
        }

        // for (int k = 0; k < dim_chunk[2]; ++k) {   // z
        //     fprintf(stderr, "# z = %g\n", k * dl);
        //     for (int j = 0; j < dim_chunk[1]; ++j) {   // y
        //         for (int i = 0; i < dim_chunk[0]; ++i) {   // x
        //             fprintf(stderr, "%.1f ", chunk[k*dim_chunk[0]*dim_chunk[1] + j*dim_chunk[0] + i] - 273.15);
                    
        //         }
        //         fprintf(stderr,"\n");
        //     }
        // }
        fprintf(stderr,"\n");
        fprintf(stderr, "For graphical rendering: python3 rendu_picture_steady.py [filename.txt] %d %d %d\n", n, m, o);
        #endif
    }

     fprintf(stderr, "end prog %d\n", my_rank);

    MPI_Finalize();


    free(chunk);

    free(new_chunk);
    free(dim_chunk);
    free(dividers);
    free(T);
    free(front_r);
    free(back_r);
    free(south_r);
    free(north_r);
    free(west_r);
    free(east_r);
    free(front_s);
    free(back_s);
    free(south_s);
    free(north_s);
    free(west_s);
    free(east_s);

    exit(EXIT_SUCCESS);
}
