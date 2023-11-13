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
static inline double update_temperature(const double *T, int u, int n, int m, int o, int i, int j, int k)
{
		/* quantity of thermal energy that must be brought to a cell to make it heat up by 1°C */
    const double cell_heat_capacity = sink_heat_capacity * sink_density * dl * dl * dl; /* J.K */
    const double dl2 = dl * dl;
    double thermal_flux = 0;

    if (i > 0)
        thermal_flux += (T[u - 1] - T[u]) * sink_conductivity * dl; /* neighbor x-1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (i < n - 1)
        thermal_flux += (T[u + 1] - T[u]) * sink_conductivity * dl; /* neighbor x+1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (j > 0)
        thermal_flux += (T[u - n] - T[u]) * sink_conductivity * dl; /* neighbor y-1 */
    else {
        /* Bottom cell: does it receive it from the CPU ? */
	if (CPU_shape(i * dl, k * dl))
            thermal_flux += CPU_TDP / CPU_surface * dl2;
        else {
            thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
            thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
        }
    }

    if (j < m - 1)
        thermal_flux += (T[u + n] - T[u]) * sink_conductivity * dl; /* neighbor y+1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (k > 0)
        thermal_flux += (T[u - n * m] - T[u]) * sink_conductivity * dl; /* neighbor z-1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (k < o - 1)
        thermal_flux += (T[u + n * m] - T[u]) * sink_conductivity * dl; /* neighbor z+1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    /* adjust temperature depending on the heat flux */
    return T[u] + thermal_flux * dt / cell_heat_capacity;
}

/* Run the simulation on the k-th xy plane.
 * v is the index of the start of the k-th xy plane in the arrays T and R. */
static inline void do_xy_plane(const double *T, double *R, int *edges, 
                               int n, int m, int o, int k, int *dim_chunk, int p, int inside) // A FAIRE
{
    int T_i = (p*dim_chunk[0])%n;
    int T_j = ((p*dim_chunk[0]/n)*dim_chunk[1])%m;
    int T_k = ( ( (p*dim_chunk[0]/n)*dim_chunk[1] )/m )*dim_chunk[2] + k;

    if (T_k == 0) // we do not modify the z = 0 plane: it is maintained at constant temperature via water-cooling
        return;

    for (int j = inside; j < dim_chunk[1]-inside; ++j) {   // y
        for (int i = inside; i < dim_chunk[0]-inside; ++i) {   // x
            int u = k*dim_chunk[1]*dim_chunk[0] + j * dim_chunk[0] + i;
            R[u] = update_temperature(T, u, n, m, o, T_i, T_j, T_k);
            ++T_i;
        }
        T_i -= dim_chunk[0] + 2*inside;
        ++T_j;
    }
}

static inline void do_xz_plane(const double *T, double *R, int *edges, 
                               int n, int m, int o, int j, int *dim_chunk, int p, int inside) // A FAIRE
{
    int T_i = (p*dim_chunk[0])%n;
    int T_j = ((p*dim_chunk[0]/n)*dim_chunk[1])%m + j;
    int T_k = ( ( (p*dim_chunk[0]/n)*dim_chunk[1] )/m )*dim_chunk[2];

    for (int k = inside; k < dim_chunk[2]-inside; ++k) {   // z
        if (T_k == 0){
            ++T_k;
            continue;
        }
        for (int i = inside; i < dim_chunk[0]-inside; ++i) {   // x
            int u = k*dim_chunk[1]*dim_chunk[0] + j * dim_chunk[0] + i;
            R[u] = update_temperature(T, u, n, m, o, T_i, T_j, T_k);
            ++T_i;
        }
        T_i -= dim_chunk[0] + 2*inside;
        ++T_k;
    }
}

static inline void do_yz_plane(const double *T, double *R,  int *edges, 
                                int n, int m, int o, int i, int *dim_chunk, int p, int inside) // A FAIRE
{
    int T_i = (p*dim_chunk[0])%n + i;
    int T_j = ((p*dim_chunk[0]/n)*dim_chunk[1])%m;
    int T_k = ( ( (p*dim_chunk[0]/n)*dim_chunk[1] )/m )*dim_chunk[2];

    for (int k = inside; k < dim_chunk[2]-inside; ++k) {   // z
        if (T_k == 0){
            ++T_k;
            continue;
        }
        for (int j = inside; j < dim_chunk[1]-inside; ++j) {   // y
            int u = k*dim_chunk[1]*dim_chunk[0] + j * dim_chunk[0] + i;
            R[u] = update_temperature(T, u, n, m, o, T_i, T_j, T_k);
            ++T_j;
        }
        T_j -= dim_chunk[1] + 2*inside;
        ++T_k;
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
            tab[0] = j;
            if (k<i) {
                tab[1] = i;
                tab[2] = k;
            }
            else {
                tab[1] = k;
                tab[2] = i;
            }
        }
        else {
            tab[0] = i;
            tab[1] = j;
            tab[2] = k;
        }
        
    }

    else if (j<k) {
        if (i<k) {
            tab[0] = k;
            if (j<i) {
                tab[1] = i;
                tab[2] = j;
            }
            else {
                tab[1] = j;
                tab[2] = i;
            }
        }
        else {
            tab[0] = i;
            tab[1] = k;
            tab[2] = j;
        }
        
    }

    return tab;
}


void malloc_blocs(int size_chunk, int *dim_chunk,
		double *chunk, double *new_chunk, double *halo,
		double *east_r, double *west_r, double *north_r, double *south_r,
        double *east_s, double *west_s, double *north_s, double *south_s) 
        // TODO 1: malloc que si nécessaire (pas les bords) -> modifier TODO 2
{
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
    halo = calloc(sizeof(double),(dim_chunk[0]+2)*(dim_chunk[1]+2)*(dim_chunk[2]+2));
    if (chunk == NULL){
        perror("malloc_blocs : malloc");
        exit(1);
    }
    east_r = malloc((dim_chunk[1])*(dim_chunk[2])*sizeof(double));
    if (east_r == NULL){
        perror("malloc_blocs : malloc");
	    exit(1);
    }
    west_r = malloc((dim_chunk[1])*(dim_chunk[2])*sizeof(double));
    if (west_r == NULL){
        perror("malloc_blocs : malloc");
	    exit(1);
    }
    north_r = malloc((dim_chunk[0])*(dim_chunk[2])*sizeof(double));
    if (north_r == NULL){
	    perror("malloc_blocs : malloc");
	    exit(1);
    }
    south_r = malloc((dim_chunk[0])*(dim_chunk[2])*sizeof(double));
   	if (south_r == NULL){
        perror("malloc_blocs : malloc");
	    exit(1);
    }
    east_s = malloc((dim_chunk[1])*(dim_chunk[2])*sizeof(double));
    if (east_s == NULL){
        perror("malloc_blocs : malloc");
	    exit(1);
    }
    west_s = malloc((dim_chunk[1])*(dim_chunk[2])*sizeof(double));
    if (west_s == NULL){
        perror("malloc_blocs : malloc");
	    exit(1);
    }
    north_s = malloc((dim_chunk[0])*(dim_chunk[2])*sizeof(double));
    if (north_s == NULL){
        perror("malloc_blocs : malloc");
	    exit(1);
    }
    south_s = malloc((dim_chunk[0])*(dim_chunk[2])*sizeof(double));
   	if (south_s == NULL){
        perror("malloc_blocs : malloc");
	    exit(1);
    }
}

void fill_halo(double *tab, int *dim_chunk, double *halo, int index) { 
    /* index = -1 (center) / 0 (east) / 1 (south) / 2 (west) / 3 (north) */
    
    switch(index) {
        
        case -1 : // center
            for (int i = 1; i  < dim_chunk[0]+1; ++i) {

                for (int j = 1; j < dim_chunk[1]+1; ++j) {

                    for (int k = 1; k < dim_chunk[2]+1; ++k) {
                        halo[i + j*(dim_chunk[0]+2) + k*(dim_chunk[0]+2)*(dim_chunk[1]+2)] 
                        = tab[(i-1) + (j-1)*dim_chunk[0] + (k-1)*dim_chunk[0]*dim_chunk[1]];
                    }
                }
            }
            break;
        
        case 0 : // east
            for (int j = 1; j < dim_chunk[1]+1; ++j) {

                for (int k = 1; k < dim_chunk[2]+1; ++k) { // i = dim_chunk[0]-1
                    halo[dim_chunk[0]-1 + j*(dim_chunk[0]+2) + k*(dim_chunk[0]+2)*(dim_chunk[1]+2)]
                        = tab[j*dim_chunk[2] + k];
                }
            }
            break;

        case 1 : // south
            for (int i = 1; i < dim_chunk[0]+1; ++i) {

                for (int k = 1 ; k < dim_chunk[2]+1; ++k) { // j = 0
                    halo[i + k*(dim_chunk[0]+2)*(dim_chunk[1]+2)] 
                        = tab[i*dim_chunk[2] + k];
                }

            }
            break;

        case 2 : // west
            for (int j = 1; j < dim_chunk[1]+1; ++j) {

                for (int k = 1; k < dim_chunk[2]+1; ++k) { // i = 0
                    halo[j*(dim_chunk[0]+2) + k*(dim_chunk[0]+2)*(dim_chunk[1]+2)]
                        = tab[j*dim_chunk[2] + k];
                }
            }
            break;
        
        case 3 : // north
            for (int i = 1; i < dim_chunk[0]+1; ++i) {

                for (int k = 1 ; k < dim_chunk[2]+1; ++k) { // j = dim_chunk[1]-1
                    halo[i + (dim_chunk[1]-1)*dim_chunk[0] + k*(dim_chunk[0]+2)*(dim_chunk[1]+2)] 
                        = tab[i*dim_chunk[2] + k];
                }

            }
            break;
            
        default :
            return;
    }
}

int main(int argc, char **argv)
{
   	
    int my_rank; /* rank of the process */
    int p; /* number of processes */
    int dest; /*rank of the dest */
    int source; /* rank of the source */
    int tag = 0; /* tag of the message */
    char message[100];
    MPI_Status status;
    //char hostname[SIZE_H_N] ;
    //gethostname(hostname, SIZE_H_N);	
       
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
    double *chunk;   // chunk(i,j,k) = chunk[i + j*chunk0 + k*chunk0*chunk1]
    double *new_chunk;

    double *east_r;  // id = 0 // east(j,k) = east[j*chunk2 + k] = chunk(dim_chunk0,j,k)
    double *south_r; // id = 1 // south(i,k) = south[i*chunk2 + k] = chunk(i,-1,k)
    double *west_r;  // id = 2 // west(j,k) = west[j*chunk2 + k] = chunk(-1,j,k)
    double *north_r; // id = 3 // north(i,k) = north[i*chunk2 + k] = chunk(i,dim_chunk1,k)

    double *east_s;  // east(j,k) = east[j*chunk2 + k] = chunk(0,j,k)
    double *west_s;  // west(j,k) = west[j*chunk2 + k] = chunk(0,j,k)
    double *north_s; // north(i,k) = north[i*chunk2 + k] = chunk(i,0,k)
    double *south_s; // south(i,k) = south[i*chunk2 + k] = chunk(i,0,k)

    double *halo;  // halo(i,j,k) = halo[i +j*chunk0+2 + k*chunk0+2*chunk1+2] 

    malloc_blocs(size_chunk, dim_chunk, chunk, new_chunk,
		    halo, east_r, west_r, north_r, south_r, 
            east_s, west_s, north_s, south_s);

    /* initially the heatsink is at the temperature of the water-cooling fluid */
    for (int u = 0; u < size_chunk; u++)
       new_chunk[u] = chunk[u] = watercooling_T + 273.15;

    /* let's go! we switch the CPU on and launch the simulation until it reaches a stationary state. */
    double t = 0;
    int n_steps = 0;
    int convergence = 0;


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
        MPI_Request rec[6];
        if (p%dividers[0] != dividers[0]-1) {
            MPI_Irecv(east_r, (dim_chunk[1])*(dim_chunk[2]),
                    MPI_DOUBLE, p+1, 0, MPI_COMM_WORLD,
                    &rec[0]);
        }
        if ((p/dividers[0])%dividers[1] != 0) {
            MPI_Irecv(south_r, (dim_chunk[0])*(dim_chunk[2]),
                MPI_DOUBLE, p-n, 1, MPI_COMM_WORLD,
                &rec[1]);
        } 
        if (p%dividers[0] != 0) {
            MPI_Irecv(west_r, (dim_chunk[1])*(dim_chunk[2]),
                    MPI_DOUBLE, p-1, 2, MPI_COMM_WORLD,
                    &rec[2]);
        }
        if ((p/dividers[0])%dividers[1] != dividers[1]-1) {
            MPI_Irecv(north_r, (dim_chunk[0])*(dim_chunk[2]),
                    MPI_DOUBLE, p+n, 3, MPI_COMM_WORLD,
                    &rec[3]);
        }

        /* create the sides for neighbours */
        // TODO 2: remplir que les bords
        for (int j = 0; j < dim_chunk[1]; ++j) { /* (y,z) pour les faces east west */
            for (int k = 0; k < dim_chunk[2]; ++k) {
                east_s[j*dim_chunk[2] + k] = chunk[dim_chunk[0]-1 + j*dim_chunk[0] + k*dim_chunk[0]*dim_chunk[1]];
                west_s[j*dim_chunk[2] + k] = chunk[j*dim_chunk[0] + k*dim_chunk[0]*dim_chunk[1]];
            }
        } 

        for (int i = 0; i < dim_chunk[0]; ++i) {
            for (int k = 0; k < dim_chunk[2]; ++k) {
                south_s[i*dim_chunk[2] + k] = chunk[i + k*dim_chunk[0]*dim_chunk[1]];
                north_s[i*dim_chunk[2] + k] = chunk[(dim_chunk[1]-1)*dim_chunk[0] + i + k*dim_chunk[0]*dim_chunk[1]];
            }
        } 

        /* sends the sides to neighbours */
        MPI_Request sen[6];
        if (p%dividers[0] != dividers[0]-1) {
            MPI_Isend(east_s, (dim_chunk[1])*(dim_chunk[2]),
                    MPI_DOUBLE, p+1, 0, MPI_COMM_WORLD,
                    &sen[0]);
        }
        if ((p/dividers[0])%dividers[1] != 0) {
            MPI_Isend(south_s, (dim_chunk[0])*(dim_chunk[2]),
                MPI_DOUBLE, p-n, 1, MPI_COMM_WORLD,
                &sen[1]);
        } 
        if (p%dividers[0] != 0) {
            MPI_Isend(west_s, (dim_chunk[1])*(dim_chunk[2]),
                    MPI_DOUBLE, p-1, 2, MPI_COMM_WORLD,
                    &sen[2]);
        }
        if ((p/dividers[0])%dividers[1] != dividers[1]-1) {
            MPI_Isend(north_s, (dim_chunk[0])*(dim_chunk[2]),
                    MPI_DOUBLE, p+n, 3, MPI_COMM_WORLD,
                    &sen[3]);
        }
        
        int *edges = calloc(sizeof(int), 4);
        if (edges == NULL) {
            perror("malloc main");
            return EXIT_FAILURE;
        }
        
        /* compute the inside of the chunk */
        for (int k = 1; k < dim_chunk[0]-1; ++k) {
            do_xy_plane(chunk, new_chunk, edges, n, m, o, k, dim_chunk, p, 1);
        }

        int indx;
        fill_halo(chunk, dim_chunk, halo, -1);

        for (int i = 0; i < 4; i++) {
            MPI_Waitany(4, rec, &indx, MPI_STATUS_IGNORE);

            edges[indx] = 1;

            switch(indx) {
                case 0 : //east
                    fill_halo(east_r, dim_chunk, halo, 0);
                    do_yz_plane(halo, new_chunk, edges, n, m, o, dim_chunk[0]-1, dim_chunk, p, 0);
                    break;

                case 1 : //south
                    fill_halo(east_r, dim_chunk, halo, 0);
                    do_xz_plane(halo, new_chunk, edges, n, m, o, 0, dim_chunk, p, 0);
                    break;

                case 2 : //west
                    fill_halo(east_r, dim_chunk, halo, 0);
                    do_yz_plane(halo, new_chunk, edges, n, m, o, 0, dim_chunk, p, 0);
                    break;

                case 3 : //north
                    fill_halo(east_r, dim_chunk, halo, 0);
                    do_xz_plane(halo, new_chunk, edges, n, m, o, dim_chunk[1]-1, dim_chunk, p, 0);
                    break;
                    
                default :
                    return EXIT_FAILURE;
            }
        }

            /* each second, we test the convergence, and print a short progress report */
        //test it locally, and sum all the partial values with MPI_allReduceand MPI_SUM
        //Be careful : the sqrt has to be done on the total sum 
        if (n_steps % ((int)(1 / dt)) == 0) {
            double delta_T = 0;
            double max = -INFINITY;

            for (int u = 0; u < n * m * o; u++) {
                delta_T += (R[u] - T[u]) * (R[u] - T[u]);
                if (R[u] > max)
                    max = R[u];
            }

            delta_T = sqrt(delta_T) / dt;
            fprintf(stderr, "t = %.1fs ; T_max = %.1f°C ; convergence = %g\n", t, max - 273.15, delta_T);
            
            if (delta_T < 0.1)
                convergence = 1;
	    }

        /* the new temperatures are in R */
        double *tmp = R;
        R = T;
        T = tmp;

        t += dt;
        n_steps += 1;
    }

#ifdef DUMP_STEADY_STATE
    printf("###### STEADY STATE; t = %.1f\n", t);
    for (int k = 0; k < o; k++) {   // z
        printf("# z = %g\n", k * dl);
        for (int j = 0; j < m; j++) {   // y
            for (int i = 0; i < n; i++) {   // x
                printf("%.1f ", T[k * n * m + j * n + i] - 273.15);
            }
            printf("\n");
        }
    }
    printf("\n");
    fprintf(stderr, "For graphical rendering: python3 rendu_picture_steady.py [filename.txt] %d %d %d\n", n, m, o);
#endif
    
    MPI_Finalize();

    exit(EXIT_SUCCESS);
}
