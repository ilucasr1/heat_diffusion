#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
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
#define CHALLENGE         /* MEDIUM is faster, and FAST is even faster (for debugging) */
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
double dt = 0.0001;
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
                               int n, int m, int o, int k, int *dim_chunk, int my_rank, int* excess)
{
    int T_i = (my_rank*dim_chunk[0])%(n+excess[0]);
    int T_j = ((my_rank*dim_chunk[0]/(n+excess[0]))*dim_chunk[1])%(m+excess[1]);
    int T_k = ( ( (my_rank*dim_chunk[0]/(n+excess[0]))*dim_chunk[1] )/(m+excess[1]) )*dim_chunk[2] + k;
    
    if (T_k == 0) // we do not modify the z = 0 plane: it is maintained at constant temperature via water-cooling
        return;

    int u;
    for (int j = 1; j < dim_chunk[1]-1; ++j) {   // y
        for (int i = 1; i < dim_chunk[0]-1; ++i) {   // x
            u = k*dim_chunk[1]*dim_chunk[0] + j * dim_chunk[0] + i;
            R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], T[u + dim_chunk[0]*dim_chunk[1]],
                                        T[u - dim_chunk[0]], T[u + dim_chunk[0]], T[u - 1], T[u + 1], n, m, o, T_i+i, T_j+j, T_k);
        }
    }
}


static inline void do_halo(const double *T, double *R,  double *front_r, double *back_r, 
                            double *south_r, double *north_r, double *west_r, double *east_r, 
                            int n, int m, int o, int *dim_chunk, int my_rank, int*excess) {
/**
 * @param T chunk at time t
 * @param R chunk at time t+1
 * @param n dim of the whole data on x
 * @param m dim of the whole data on y
 * @param o dim of the whole data on z
 * @param my_rank current processor's rank
 */

    int T_i = (my_rank*dim_chunk[0])%(n+excess[0]);
    int T_j = ((my_rank*dim_chunk[0]/(n+excess[0]))*dim_chunk[1])%(m+excess[1]);
    int T_k = (((my_rank*dim_chunk[0]/(n+excess[0]))*dim_chunk[1])/(m+excess[1]))*dim_chunk[2];

    int     u;
    double  f;
    double  b;
    double  s;
    double no;
    double  w;
    double  e;

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
                                
    // coin haut droit
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
        w = (west_r == NULL)  ? 0 : west_r[k*dim_chunk[1]];
        e = (dim_chunk[0] > 1) ? T[u + 1] :
            ((east_r == NULL) ? 0 : east_r[k*dim_chunk[1]]);
        
        R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], T[u + dim_chunk[0]*dim_chunk[1]], s,
                                  no, w, e, n, m, o, T_i, T_j, T_k + k);

        //i = dim_chunk[0]-1
        u = k*dim_chunk[1]*dim_chunk[0] + dim_chunk[0] - 1;
        s = (south_r == NULL) ? 0 : south_r[k*dim_chunk[0] + dim_chunk[0] - 1];
        no = (dim_chunk[1] > 1) ? T[u + dim_chunk[0]] : 
             ((north_r == NULL) ? 0 : north_r[k*dim_chunk[0] + dim_chunk[0] - 1]);
        w = (dim_chunk[0] > 1) ? T[u - 1] :
            ((west_r == NULL) ? 0 : west_r[k*dim_chunk[1]]);
        e = (east_r == NULL)  ? 0 : east_r[k*dim_chunk[1]];

        R[u] = update_temperature(T[u], T[u - dim_chunk[0]*dim_chunk[1]], T[u + dim_chunk[0]*dim_chunk[1]], s,
                                   no, w, e, n, m, o, T_i + dim_chunk[0] - 1, T_j, T_k + k);
    
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


int load_data(double *new_chunk, double *chunk, int *dim_chunk, int my_rank, int p, int *n_steps, double *t) 
{
    FILE *f = fopen("save.txt", "r");
    if (f == NULL) {
        perror("Erreur lors de l'ouverture du fichier");
        return EXIT_FAILURE;
    }

    int save = 0;
    fscanf(f, "%d", &save);

    /* case with no saving */
    if (save == 0 || save != p) {
        for (int u = 0; u < dim_chunk[0]*dim_chunk[1]*dim_chunk[2]; u++) {
            new_chunk[u] = watercooling_T + 273.15;
            chunk[u] = watercooling_T + 273.15;
        }
        return EXIT_SUCCESS;

    }
    
    /* case with a previous save run on same nbr of processus */
    int rank = -1;
    double val = 0.;
    int nbr_char = 2048;

    #ifdef CHALLENGE
    nbr_char *= 2;
    #endif

    char *line = malloc(nbr_char*sizeof(char));
    if (line == NULL) {
        fprintf(stderr, "load malloc failed\n");
        return EXIT_FAILURE;
    }

    fscanf(f, "%lf %d\n", t, n_steps);
    
    while (fgets(line, nbr_char, f) != NULL) {
        if (sscanf(line, "%d", &rank) == EOF) {
            fprintf(stderr, "load scanf failure\n");
            fclose(f);
            return EXIT_FAILURE;
        }

        if (rank != my_rank) {
            for (int k = 0; k < dim_chunk[2]; ++k) {
                // read z = ...
                if (fgets(line, nbr_char*sizeof(char), f) == NULL) {
                    fprintf(stderr, "Impossible to jump %d lines. The file is too short.\n", dim_chunk[2]);
                    fclose(f);
                    return EXIT_FAILURE;
                }
                // read the data
                for (int j = 0; j < dim_chunk[1]; ++j) {
                    if (fgets(line, nbr_char*sizeof(char), f) == NULL) {
                        fprintf(stderr, "Impossible to jump %d lines. The file is too short.\n", dim_chunk[2]);
                        fclose(f);
                        return EXIT_FAILURE;
                    }
                }
            }
        }

        
        else {
            for (int k = 0; k < dim_chunk[2]; ++k) {
                // read z = ...
                if (fgets(line, nbr_char*sizeof(char), f) == NULL) {
                    fprintf(stderr, "Impossible to jump %d lines. The file is too short.\n", dim_chunk[2]);
                    fclose(f);
                    return EXIT_FAILURE;
                }

                // read the data
                for (int j = 0; j < dim_chunk[1]; ++j) {
                    if (fgets(line, nbr_char*sizeof(char), f) == NULL) {
                        fprintf(stderr, "Impossible to jump %d lines. The file is too short.\n", dim_chunk[2]);
                        fclose(f);
                        return EXIT_FAILURE;
                    }
                    char *tmp = malloc(nbr_char);
                    strcpy(tmp, line);
                    for (int i = 0; i < dim_chunk[0]; ++i) {

                        if (i == 0)
                            val = atof(strtok(tmp, " \n"));
                        else
                            val = atof(strtok(NULL, " \n"));

                        chunk[k*dim_chunk[1]*dim_chunk[0] + j*dim_chunk[0] + i] = val;
                        new_chunk[k*dim_chunk[1]*dim_chunk[0] + j*dim_chunk[0] + i] = val;
                    }
                    free(tmp);
                    sscanf(line, "\n");
                }
                sscanf(line, "\n");
            }
            return EXIT_SUCCESS;
        }
    }
    
    return 0;
}


int store_data(int my_rank, int p, int *dim_chunk, double *chunk, double t, int n_steps) {
    
    FILE *f;
    
    if (my_rank == 0) {
        f = fopen("save.txt", "w");
        if (f == NULL) {
            perror("Erreur lors de l'ouverture du fichier");
            return EXIT_FAILURE;
        }
        fprintf(f, "%d\n%.3lf %d\n", p, t, n_steps);
    }
    else {
        f = fopen("save.txt", "a");
        if (f == NULL) {
            perror("Erreur lors de l'ouverture du fichier");
            return EXIT_FAILURE;
        }
    }
    
    int msg = 0;
    for (int i = 0; i < my_rank; ++i) {
        MPI_Recv(&msg, 1, MPI_INT, i , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    fprintf(f,"%d\n", my_rank);
    for (int k = 0; k < dim_chunk[2]; ++k) {
        fprintf(f, "# z = %g\n", k * dl);
        for (int j = 0; j < dim_chunk[1]; ++j) {
            for (int i = 0; i < dim_chunk[0]; ++i) {
                fprintf(f, "%lf ", chunk[k*dim_chunk[0]*dim_chunk[1] + j*dim_chunk[0] + i]);
            }
            fprintf(f, "\n");
        }
    }
    
    for (int i = my_rank+1; i < p; ++i) {
        MPI_Send(&msg, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }

    fclose(f);
    return 0;
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
            tab[0] = j;
            if (k<i) {
                tab[2] = i;
                tab[1] = k;
            }
            else {
                tab[2] = k;
                tab[1] = i;
            }
        }
        else {
            tab[0] = i;
            tab[2] = j;
            tab[1] = k;
        }
        
    }

    else {
        if (i<k) {
            tab[0] = k;
            if (j<i) {
                tab[2] = i;
                tab[1] = j;
            }
            else {
                tab[2] = j;
                tab[1] = i;
            }
        }
        else {
            tab[0] = i;
            tab[2] = k;
            tab[1] = j;
        }
        
    }

    return tab;
}


int main(int argc, char **argv)
{
    time_t begin = time(NULL);
    int my_rank; /* rank of the process */
    int p; /* number of processes */
    int dim = 3;
    if (argc > 1)
        dim = atoi(argv[1]);

    /* Initialisation */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    CPU_surface = CPU_contact_surface();
    double V = L * l * E;
    int n = ceil(L / dl);
    int m = ceil(E / dl);
    int o = ceil(l / dl);
    int *dividers = get_3_dividers(p);

    int *dim_chunk = malloc(3*sizeof(int));
    if (dim_chunk == NULL){
    	perror("malloc");
	exit(1);
    }
    if (dim <= 2) {
        if (dim == 1) {
            dividers[0] *= dividers[2];
            dividers[2] = 1;
        }
        dividers[0] *= dividers[1];
        dividers[1] = 1;
    }    
    
    dividers[0]=15;
    dividers[1]=4;
    dividers[2]=12;
    dim_chunk[0] = ceil(n/(float)dividers[0]);
    dim_chunk[1] = ceil(m/(float)dividers[1]);
    dim_chunk[2] = ceil(o/(float)dividers[2]);
    
    int *excess = calloc(3,sizeof(int));
    if (excess == NULL) {
        perror("main malloc : excess");
        return EXIT_FAILURE;
    }
    excess[0] = dim_chunk[0]*dividers[0] - n;
    excess[1] = dim_chunk[1]*dividers[1] - m;
    excess[2] = dim_chunk[2]*dividers[2] - o;
    
    if (my_rank == 0) {
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
        fprintf(stderr, "MODELISATION\n");
        fprintf(stderr, "\texcess %d %d %d\n", excess[0], excess[1], excess[2]);
        fprintf(stderr, "\tn %d, dim_chunk[0] %d\n", n, dim_chunk[0]);
        fprintf(stderr, "\tm %d, dim_chunk[1] %d\n", m, dim_chunk[1]);
        fprintf(stderr, "\to %d, dim_chunk[2] %d\n", o, dim_chunk[2]);
    }

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

    chunk = calloc(size_chunk, sizeof(double));
    if (chunk == NULL){
        perror("malloc_blocs : malloc");
	    exit(1);
    }
    new_chunk = calloc(size_chunk, sizeof(double));
    if (new_chunk == NULL){
        perror("malloc_blocs : malloc");
	    exit(1);
    }
    
    if (dividers[2] > 1) {
        front_r = calloc((dim_chunk[0])*(dim_chunk[1]), sizeof(double));
        if (front_r == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        back_r = calloc((dim_chunk[0])*(dim_chunk[1]), sizeof(double));
        if (back_r == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        front_s = calloc((dim_chunk[0])*(dim_chunk[1]), sizeof(double));
        if (front_s == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        back_s = calloc((dim_chunk[0])*(dim_chunk[1]), sizeof(double));
        if (back_s == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }

    }

    if (dividers[1] > 1) {
        south_r = calloc((dim_chunk[0])*(dim_chunk[2]), sizeof(double));
        if (south_r == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        north_r = calloc((dim_chunk[0])*(dim_chunk[2]), sizeof(double));
        if (north_r == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        south_s = calloc((dim_chunk[0])*(dim_chunk[2]), sizeof(double));
        if (south_s == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        north_s = calloc((dim_chunk[0])*(dim_chunk[2]), sizeof(double));
        if (north_s == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
    }

    if (dividers[0] > 1) {
        west_r = calloc((dim_chunk[1])*(dim_chunk[2]), sizeof(double));
        if (west_r == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        east_r = calloc((dim_chunk[1])*(dim_chunk[2]), sizeof(double));
        if (east_r == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        west_s = calloc((dim_chunk[1])*(dim_chunk[2]), sizeof(double));
        if (west_s == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
        east_s = calloc((dim_chunk[1])*(dim_chunk[2]), sizeof(double));
        if (east_s == NULL){
            perror("malloc_blocs : malloc");
            exit(1);
        }
    }
    
    /* let's go! we switch the CPU on and launch the simulation until it reaches a stationary state. */
    double t = 0;
    int n_steps = 0;
    int convergence = 0;

    /* initially the heatsink is at the temperature of the water-cooling fluid */
    load_data(new_chunk, chunk, dim_chunk, my_rank, p, &n_steps, &t);

    /* simulating time steps */
    while (convergence == 0) {

        /* Update all cells. xy planes are processed, for increasing values of z. */
	/* Ask for all the Irecv (check the border conditions)
	 * prepare the data and Isend them
	 * Start the calculations for the inside of the block
	 * wait for the Irecv and calculate each as soon as it's sent
	 * maybe wait for the Isend before going to t+1 because the values Isent must not change before they're effectively sent
	 * */

        /* receive the sides from others */
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
            MPI_Irecv(front_r, (dim_chunk[0])*(dim_chunk[1]),
                    MPI_DOUBLE, my_rank-dividers[0]*dividers[1], 0, MPI_COMM_WORLD,
                    &rec[cpt++]);
        }
        if (test_b) {
            MPI_Irecv(back_r, (dim_chunk[0])*(dim_chunk[1]),
                    MPI_DOUBLE, my_rank+dividers[0]*dividers[1], 0, MPI_COMM_WORLD,
                    &rec[cpt++]);
        }
        if (test_s) {
            MPI_Irecv(south_r, (dim_chunk[0])*(dim_chunk[2]),
                    MPI_DOUBLE, my_rank-dividers[0], 0, MPI_COMM_WORLD,
                    &rec[cpt++]);
        } 
        if (test_n) {
            MPI_Irecv(north_r, (dim_chunk[0])*(dim_chunk[2]),
                    MPI_DOUBLE, my_rank+dividers[0], 0, MPI_COMM_WORLD,
                    &rec[cpt++]);
        }
        if (test_w) {
            MPI_Irecv(west_r, (dim_chunk[1])*(dim_chunk[2]),
                    MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD,
                    &rec[cpt++]);
        }
        if (test_e) {
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
        
        /* compute the inside of the chunk */
        for (int k = 1; k < dim_chunk[2]-1; ++k) {
            do_xy_plane(chunk, new_chunk, n, m, o, k, dim_chunk, my_rank, excess);
        }
        
        MPI_Waitall(nb_neigh, rec, MPI_STATUSES_IGNORE);
        
        do_halo(chunk, new_chunk, front_r, back_r, south_r, north_r, west_r, east_r, n, m, o, dim_chunk, my_rank, excess);
       
       
        /* each second, we test the convergence, and print a short progress report */
        //test it locally, and sum all the partial values with MPI_allReduceand MPI_SUM
        //Be careful : the sqrt has to be done on the total sum 

        //each processor compute the sum
        if (n_steps % ((int)(1 / dt)) == 0) {
            double delta_T = 0;
            double max = -INFINITY;
            int cpt = 0;
            for (int u = 0; u < dim_chunk[0]*dim_chunk[1]*dim_chunk[2]; ++u) {
                int cond1 = (n%dividers[0] != 0 && my_rank%dividers[0] == dividers[0]-1 && 
                         u%dim_chunk[0] >= dim_chunk[0]-excess[0]);

                int cond2 = (m % dividers[1] != 0 && (my_rank/dividers[0])%dividers[1] == dividers[1]-1 && 
                        (u/dim_chunk[0])%dim_chunk[1] >= dim_chunk[1]-excess[1]);
                        
                int cond3 = (o % dividers[2] != 0 && my_rank/(dividers[0]*dividers[1]) == dividers[2]-1 && 
                         u/(dim_chunk[0]*dim_chunk[1]) >= dim_chunk[2]-excess[2]);
                         
                if (cond1 || cond2 || cond3) continue;
                cpt++;
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

        if (n_steps % (int)(1/dt)*10 == 0) {
            store_data(my_rank, p, dim_chunk, chunk, t, n_steps);
        }
    }

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
    
    double *T = NULL;
    double *R = NULL;
    int* sizes_chunk = NULL;
    int *positions = NULL;
    if (my_rank == 0) {
        T = calloc(n*m*o, sizeof(double));
        if (T == NULL) {
            perror("main: malloc T");
            return EXIT_FAILURE;
        }

        R = calloc(n*m*o, sizeof(double));
        if (R == NULL) {
            perror("main: malloc R");
            return EXIT_FAILURE;
        }

        sizes_chunk = malloc(p*sizeof(int));
        if (sizes_chunk == NULL) {
            perror("main : malloc sizes");
            return EXIT_FAILURE;
        }
        positions = calloc(p,sizeof(int));
        if (positions == NULL) {
            perror("main : calloc positions");
            return EXIT_FAILURE;
        }
    }
    
    int tmp_size_chunk = (dim_chunk[0]- ((my_rank%dividers[0] == dividers[0]-1 && n % dividers[0] != 0)?excess[0]:0))*
                         (dim_chunk[1]- (((my_rank/dividers[0])%dividers[1] == dividers[1]-1 && m % dividers[1] != 0)?excess[1]:0))*
                         (dim_chunk[2]- ((my_rank/(dividers[0]*dividers[1]) == dividers[2]-1 && o % dividers[2] != 0)?excess[2]:0));

    if ( (my_rank%dividers[0])*dim_chunk[0] >= n ||
        ((my_rank/dividers[0])%dividers[1])*dim_chunk[1] >= m ||
         (my_rank/(dividers[0]*dividers[1]))*dim_chunk[2] >= o)
         tmp_size_chunk = 0;

    MPI_Gather(&tmp_size_chunk, 1, MPI_INT, sizes_chunk, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (my_rank == 0) {
        for (int i = 1; i < p; ++i) {
            positions[i] = positions[i-1] + sizes_chunk[i-1];
        }
    }

    int face_right = (n%dividers[0]   != 0 && my_rank%dividers[0] == dividers[0]-1);
    int face_north = (m % dividers[1] != 0 && (my_rank/dividers[0])%dividers[1] == dividers[1]-1);
    int face_back  = (o % dividers[2] != 0 && my_rank/(dividers[0]*dividers[1]) == dividers[2]-1);
    if ((face_right || face_north || face_back) && tmp_size_chunk > 0) {
        double *tmp_chunk = calloc( (dim_chunk[0]-face_right*excess[0])*
                                    (dim_chunk[1]-face_north*excess[1])*
                                    (dim_chunk[2]-face_back*excess[2]),sizeof(double));
        if (tmp_chunk == NULL) {
            perror("main : malloc tmp_chunk");
            return EXIT_FAILURE;
        }
        for (int k = 0; k < dim_chunk[2]-face_back*excess[2]; ++k) {
            for (int j = 0; j < dim_chunk[1]-face_north*excess[1]; ++j) {
                for (int i = 0; i < dim_chunk[0]-face_right*excess[0]; ++i) {
                    tmp_chunk[i + j*(dim_chunk[0]-face_right*excess[0])
                                + k*(dim_chunk[0]-face_right*excess[0])*(dim_chunk[1]-face_north*excess[1])] 
                    = chunk[i + j*(dim_chunk[0]) 
                            + k*(dim_chunk[0])*(dim_chunk[1])];
                }
            }
        }
        free(chunk);
        chunk = tmp_chunk;
    }
    
    MPI_Gatherv(chunk, tmp_size_chunk, MPI_DOUBLE,T, sizes_chunk, positions,MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    int size_line = dim_chunk[0]*dividers[0] - excess[0];
    int size_face = (dim_chunk[0]*dividers[0] - excess[0])*(dim_chunk[1]*dividers[1] - excess[1]);
    if (my_rank == 0) {
        for (int d2 = 0; d2 < dividers[2]; ++d2) {
            for (int d1 = 0; d1 < dividers[1]; ++d1) {
                for (int d0 = 0; d0 < dividers[0] ; ++d0) {// parcourt 0,1,2,...,p-1
                    int rank  = d0 + d1*dividers[0] + d2*dividers[0]*dividers[1];
                    if (sizes_chunk[rank] == 0) continue;
                    int cond0 = (n%dividers[0]   != 0 && rank%dividers[0] == dividers[0]-1);//east
                    int cond1 = (m % dividers[1] != 0 && (rank/dividers[0])%dividers[1] == dividers[1]-1);//north
                    int cond2 = (o % dividers[2] != 0 && rank/(dividers[0]*dividers[1]) == dividers[2]-1);//back
                    for (int k = 0; k < ((cond2)?dim_chunk[2]-excess[2]:dim_chunk[2]); ++k) {
                        for (int j = 0; j < ((cond1)?dim_chunk[1]-excess[1]:dim_chunk[1]); ++j) {
                            for (int i = 0; i < ((cond0)?dim_chunk[0]-excess[0]:dim_chunk[0]); ++i) {//parcourt T[i,j,k]

                                int u = d2*dim_chunk[2]*size_face
                                        + d1*dim_chunk[1]*size_line
                                        + d0*dim_chunk[0]
                                        + k*(size_face)
                                        + j*(size_line) + i;

                                int v = positions[d0 +d1*dividers[0] + d2*dividers[0]*dividers[1]] 
                                        + k*(dim_chunk[0]-excess[0]*cond0)*(dim_chunk[1]-excess[1]*cond1)
                                        + j*(dim_chunk[0]-excess[0]*cond0) 
                                        + i;
                                
                                R[u] = T[v];


                            }
                        }
                    }
                }
            }
        }

        #ifdef DUMP_STEADY_STATE
        FILE *f = fopen("data.txt", "w");
        fprintf(stderr,"###### STEADY STATE; t = %.1f\n", t);
        for (int k = 0; k < o; ++k) {   // z
            fprintf(f, "# z = %g\n", k * dl);
            for (int j = 0; j < m; ++j) {   // y
                for (int i = 0; i < n; ++i) {   // x
                    fprintf(f,"%.3f ", R[k * n * m + j * n + i] - 273.15);
                }
                fprintf(f,"\n");
            }
        }
        fclose(f);

        fprintf(stderr,"\n");
        fprintf(stderr, "For graphical rendering: python3 rendu_picture_steady.py data.txt %d %d %d\n", n, m, o);
        #endif
    }
    
    FILE *f = fopen("save.txt", "w");
    fprintf(f, "0");
    fclose(f);

    MPI_Finalize();

    free(chunk);
    free(excess);
    free(new_chunk);
    free(dim_chunk);
    free(dividers);
    free(positions);
    free(sizes_chunk);
    free(T);
    free(R);

    time_t end = time(NULL);
    if (my_rank == 0) {
        FILE *f = fopen("data_test.txt", "a");
        fprintf(f, "CHALLENGE %d processus : Temps d'execution : %ld secondes\n",p, (end-begin));
        fclose(f);
    }
    exit(EXIT_SUCCESS);
}
