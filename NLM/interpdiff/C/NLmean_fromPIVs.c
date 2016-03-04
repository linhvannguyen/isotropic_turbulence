/*
 * compile: gcc -O3 -std=c99 -o out -fopenmp NLM_interpdiff.c -lm -I/usr/include/netcdf-3/ -L/usr/lib64/ -lnetcdf -lnetcdf_c++
 * in the terminal: export OMP_NUM_THREADS=3
*/

#include<stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include <omp.h>

/* This is the name of the data file we will read. */
#define FILENAME_RD "/data/ISOTROPIC/NLM/interpdiff/Udiff_spacespacing3.nc"
#define FILENAME_WR "/data/ISOTROPIC/NLM/interpdiff/NLM_interpdiff_simpatch0_accumpatch0_searchbox0_tau0100.nc"

/* all constants */
#define N_HR 96

#define SCALE_FACTOR_SPACE 3
#define SCALE_FACTOR_TIME 4

#define ESTIMATE_RANGE_Y_MIN 0
#define ESTIMATE_RANGE_Y_MAX (N_HR - 1)
#define ESTIMATE_RANGE_Z_MIN 0
#define ESTIMATE_RANGE_Z_MAX (N_HR - 1)

#define SIM_HAFTSIZE 12 // dont go further than 12
#define ACC_HAFTSIZE 12
#define NEIGHBOR_HAFTSIZE 8

#define SIM_FULLSIZE (2 * SIM_HAFTSIZE + 1)
#define ACC_FULLSIZE (2 * ACC_HAFTSIZE + 1)
#define NEIGHBOR_FULLSIZE (2 * NEIGHBOR_HAFTSIZE + 1)

#define TAU 0.1

#define NUM_VARS 1
#define NUM_SCALES 2


#define NUM_3DSNAPS 1  /* #3D snapshots */
#define NUM_BLOCKS 2 /* #(1:SCALE_FACTOR_TIME:N_HR) - 1*/
#define NUM_2DSNAPS (SCALE_FACTOR_TIME * NUM_BLOCKS + 1) /* #2D snapshots in each 3D block */
#define NDIMS 4

/* Handle errors by printing an error message and exiting with a non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

/* **********************************************************************************/
/* ****************************** USEFUL FUNCTIONS **********************************/
/* **********************************************************************************/

/*
 * get_onesnap: take part of a big array(arr1) and put to small one (arr2): arr2 = arr1[id_start:id_end]
 */
void get_onesnap(double *arr1,double *arr2, int id_start, int id_end)
{
    for (int i = id_start; i < id_end + 1; i++)
        arr2[i - id_start] = arr1[i];
}


/*
 * put_onesnap: assign small array (arr2) into biger one (arr1): arr1[id_start:id_end] = arr2
 */
void put_onesnap(double *arr1,double *arr2, int id_start, int id_end)
{
    for (int i = id_start; i < id_end + 1; i++)
        arr1[i] = arr2[i - id_start];
}


/*
 * norm_by_weight: normalize x[dim] by weight W[dim]
 */

void norm_by_weight(int dim, double *x, double *W)
{
    for (int k = 0; k < dim; k++)
        x[k] = x[k]/W[k];
}


/* **********************************************************************************/
/* ****************************** NETCDF UTILS **************************************/
/* **********************************************************************************/

/*
 * creat_netcdf: create the netcdf file [filename] contain [num_vars] variables
 * variable names are [varname]
*/

void create_netcdf(char *filename, int num_vars, char *varname[num_vars])
{
    int ncid_wr, retval_wr;
    int vel_varid_wr;
    int Nt, Nx, Ny, Nz;
    int dimids[NDIMS];

    /* Create the file. */
    if ((retval_wr = nc_create(filename, NC_CLOBBER, &ncid_wr)))
       ERR(retval_wr);

    /* Define the dimensions. The record dimension is defined to have
     * unlimited length - it can grow as needed.*/
    if ((retval_wr = nc_def_dim(ncid_wr, "Ny", N_HR, &Ny)))
        ERR(retval_wr);
    if ((retval_wr = nc_def_dim(ncid_wr, "Nz", N_HR, &Nz)))
        ERR(retval_wr);
    if ((retval_wr = nc_def_dim(ncid_wr, "Nt", NC_UNLIMITED, &Nt)))
        ERR(retval_wr);

    /* Define the netCDF variables for the data. */
    dimids[0] = Nt;
    dimids[1] = Nx;
    dimids[2] = Ny;
    dimids[3] = Nz;

    for (int i = 0; i<num_vars; i++)
    {
        if ((retval_wr = nc_def_var(ncid_wr, varname[i], NC_FLOAT, NDIMS, dimids, &vel_varid_wr)))
            ERR(retval_wr);
    }

    /* End define mode (SHOULD NOT FORGET THIS!). */
    if ((retval_wr = nc_enddef(ncid_wr)))
        ERR(retval_wr);

    /* Close the file. */
    if ((retval_wr = nc_close(ncid_wr)))
        ERR(retval_wr);
    printf("\n *** SUCCESS creating file: %s!\n", filename);
}


/*
 * write_netcdf:
 * write into [filename], variable [varname] [snap_end - snap_start + 1 ] snapshots [snaps] started at [snap_start]
*/

void write_netcdf(char *filename, char *varname, size_t *start, size_t *count, double *snaps)
{
    int ncid_wr, retval_wr;
    int vel_varid_wr;

    /* Open the file. NC_WRITE tells netCDF we want read-only access to the file.*/
    if ((retval_wr = nc_open(filename, NC_WRITE, &ncid_wr)))
        ERR(retval_wr);

    /* Get variable*/
    if ((retval_wr = nc_inq_varid(ncid_wr, varname, &vel_varid_wr)))
        ERR(retval_wr);;

    /* Put variable*/
    if ((retval_wr = nc_put_vara_double(ncid_wr, vel_varid_wr, start, count, &snaps[0])))
        ERR(retval_wr);

    /* Close the file. */
    if ((retval_wr = nc_close(ncid_wr)))
        ERR(retval_wr);

    printf("\n *** SUCCESS writing variables \"%s\" to \"%s\"!\n", varname, filename);
}


/*
 * read_netcdf: read from [filename], variable [varname] [snap_end - snap_start + 1 ] snapshots [snaps]
 * started at [snap_start]
*/

void read_netcdf(char *filename, char *varname, size_t *start, size_t *count, double *snaps)
{
    int ncid_rd, retval_rd;
    int vel_varid_rd;

    /*  ******** PREPARE TO READ ************* */
    /* Open the file. NC_NOWRITE tells netCDF we want read-only access to the file.*/
    if ((retval_rd = nc_open(filename, NC_NOWRITE, &ncid_rd)))
        ERR(retval_rd);

    /* Get the varids of the velocity in netCDF */
    if ((retval_rd = nc_inq_varid(ncid_rd, varname, &vel_varid_rd)))
        ERR(retval_rd);

    if ((retval_rd = nc_get_vara_double(ncid_rd, vel_varid_rd, start, count, &snaps[0])))
        ERR(retval_rd);

    /* Close the file, freeing all resources. */
    if ((retval_rd = nc_close(ncid_rd)))
        ERR(retval_rd);

    printf("\n *** SUCCESS reading variables \"%s\" from \"%s\" \n", varname, filename);
}


/* **********************************************************************************/
/* ****************************** ESTIMATE_DISTANCE *********************************/
/* **********************************************************************************/
/*
 * estimate_distance: estimate the distances between ref patch and moving patches (prev and after)
 * patches are of fixed size (2*SIM_HAFTSIZE+1) x (2*SIM_HAFTSIZE+1)
 * reference patch are centered at [center_ref_idy, center_ref_idz]
 * moving patches are centered at [center_moving_idy, center_moving_idz]
 * dist_all contain 2 elements: distances to moving patches in the prev and after plane
 * x_ref: reference plane
 * x_prev: previous plane
 * x_after: plane after
 * ref_ids_y(z): indices of points in reference patch
 * moving_ids_y(z): indices of points in moving patch
 */
void generate_grids(int *gridpatches_y, int *gridpatches_z, int * acc_ids)
{
    int neighbor_id, sim_id;
    int gridyoffset_neighbor[NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE], gridzoffset_neighbor[NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE];

    for (int m = 0; m < NEIGHBOR_FULLSIZE; m++)
    {
        for (int n = 0; n < NEIGHBOR_FULLSIZE; n++)
        {
            gridyoffset_neighbor[m * NEIGHBOR_FULLSIZE + n] = m - NEIGHBOR_HAFTSIZE;
            gridzoffset_neighbor[m * NEIGHBOR_FULLSIZE + n] = n - NEIGHBOR_HAFTSIZE;
        }
    }

    int gridyoffset_sim[SIM_FULLSIZE * SIM_FULLSIZE], gridzoffset_sim[SIM_FULLSIZE * SIM_FULLSIZE];
    for (int p = 0; p < SIM_FULLSIZE; p++)
    {
        for (int q = 0; q < SIM_FULLSIZE; q++)
        {
            gridyoffset_sim[p * SIM_FULLSIZE + q] = p - SIM_HAFTSIZE;
            gridzoffset_sim[p * SIM_FULLSIZE + q] = q - SIM_HAFTSIZE;
        }
    }

    int grid_sim[SIM_FULLSIZE][SIM_FULLSIZE];
    for (int p = 0; p < SIM_FULLSIZE; p++)
        for (int q = 0; q < SIM_FULLSIZE; q++)
            grid_sim[p][q] = p * SIM_FULLSIZE + q;
    for (int p = 0; p < ACC_FULLSIZE; p++)
        for (int q = 0; q < ACC_FULLSIZE; q++)
            acc_ids[p * ACC_FULLSIZE + q] = grid_sim[SIM_HAFTSIZE - ACC_HAFTSIZE + p][SIM_HAFTSIZE - ACC_HAFTSIZE + q];

    int valy, valz;
    long int grid_id;
    for (int i = 0; i < N_HR; i++)
    {
        for (int j = 0; j < N_HR; j++)
        {
            for (int neighbor_id = 0; neighbor_id < NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE; neighbor_id++)
            {
                for (int sim_id = 0; sim_id < SIM_FULLSIZE * SIM_FULLSIZE; sim_id++)
                {
                    grid_id = i * N_HR * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE
                            + j * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE
                            + neighbor_id * SIM_FULLSIZE * SIM_FULLSIZE + sim_id;

                    valy = i + gridyoffset_neighbor[neighbor_id] + gridyoffset_sim[sim_id];
                    valz = j + gridzoffset_neighbor[neighbor_id] + gridzoffset_sim[sim_id];

                    if (valy < 0)
                        gridpatches_y[grid_id] = (N_HR - 1) + valy;
                    else if (valy > (N_HR - 1))
                        gridpatches_y[grid_id] = valy - (N_HR - 1);
                    else
                        gridpatches_y[grid_id] = valy;

                    if (valz < 0)
                        gridpatches_z[grid_id] = (N_HR - 1) + valz;
                    else if (valz > (N_HR - 1))
                        gridpatches_z[grid_id] = valz - (N_HR - 1);
                    else
                        gridpatches_z[grid_id] = valz;
                }
            }
        }
    }
    //printf("\n  gridpatches_z: %i \n", gridpatches_y[0]);
}


/* **********************************************************************************/
/* ****************************** ESTIMATE_DISTANCE *********************************/
/* **********************************************************************************/
/*
 * estimate_distance: estimate the distances between ref patch and moving patches (prev and after)
 * patches are of fixed size (2*SIM_HAFTSIZE+1) x (2*SIM_HAFTSIZE+1)
 * reference patch are centered at [center_ref_idy, center_ref_idz]
 * moving patches are centered at [center_moving_idy, center_moving_idz]
 * dist_all contain 2 elements: distances to moving patches in the prev and after plane
 * x_ref: reference plane
 * x_prev: previous plane
 * x_after: plane after
 * ref_ids_y(z): indices of points in reference patch
 * moving_ids_y(z): indices of points in moving patch
 */
/*void fusion(double *x_NLM, double *weight_NLM, double *x_ref, double *x_moving, double *x_fusion,
              int gridpatches_y[N_HR][N_HR][NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE][SIM_FULLSIZE * SIM_FULLSIZE],
              int gridpatches_z[N_HR][N_HR][NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE][SIM_FULLSIZE * SIM_FULLSIZE],
              int acc_ids[NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE], int est_idy, int est_idz)*/
void NLmean(double *x_NLM, double *weight_NLM, double *x_ref, double *x_moving, double *x_fusion, int *gridpatches_y, int *gridpatches_z, int *acc_ids)
{
    double norm_fact = 1.0/((double) (SIM_FULLSIZE * SIM_FULLSIZE));
    int ri = NEIGHBOR_HAFTSIZE * NEIGHBOR_FULLSIZE + NEIGHBOR_HAFTSIZE;

    int est_idy;
    #pragma omp parallel for private (est_idy)
    for (est_idy = 0; est_idy < N_HR; est_idy++)
        for (int est_idz = 0; est_idz < N_HR; est_idz++)
            for (int ni = 0; ni < NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE; ni++)
            {
                int ref_idy, ref_idz, moving_idy, moving_idz;
                double du;
                double d = 0.0;
                long int grid_rid, grid_nid;

                for (int si = 0; si < SIM_FULLSIZE * SIM_FULLSIZE; si++)
                {
                    grid_rid = est_idy * N_HR * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE
                             + est_idz * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE + ri * SIM_FULLSIZE * SIM_FULLSIZE + si ;
                    grid_nid = est_idy * N_HR * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE
                             + est_idz * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE + ni * SIM_FULLSIZE * SIM_FULLSIZE + si;

                    ref_idy = gridpatches_y[grid_rid];
                    moving_idy = gridpatches_y[grid_nid];
                    ref_idz = gridpatches_z[grid_rid];
                    moving_idz = gridpatches_z[grid_nid];

                    //compute distance btw reference patch and fusion patch
                    du = x_ref[ref_idy * N_HR + ref_idz] - x_moving[moving_idy * N_HR + moving_idz];
                    d = d + norm_fact*du*du;
                }

                double w = exp(-d/(2.0*TAU*TAU));
                for(int k = 0; k < ACC_FULLSIZE * ACC_FULLSIZE; k++)
                {
                    int ai =  acc_ids[k];
                    grid_rid = est_idy * N_HR * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE
                             + est_idz * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE + ri * SIM_FULLSIZE * SIM_FULLSIZE + ai ;
                    grid_nid = est_idy * N_HR * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE
                             + est_idz * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE + ni * SIM_FULLSIZE * SIM_FULLSIZE + ai;

                    ref_idy = gridpatches_y[grid_rid];
                    moving_idy = gridpatches_y[grid_nid];
                    ref_idz = gridpatches_z[grid_rid];
                    moving_idz = gridpatches_z[grid_nid];

                    x_NLM[ref_idy * N_HR + ref_idz] = x_NLM[ref_idy * N_HR + ref_idz] + w*x_fusion[moving_idy * N_HR + moving_idz];
                    weight_NLM[ref_idy * N_HR + ref_idz] = weight_NLM[ref_idy * N_HR + ref_idz] + w;
                }
            }
}


/* **********************************************************************************/
/* ********************************** MAIN FUNCTION *********************************/
/* **********************************************************************************/

int main()
{
    /* Creat the file to save results */
    char *varnames[NUM_VARS] = {"x_rec_all"};
    create_netcdf(FILENAME_WR, NUM_VARS, varnames);

    /* Allocate memory */
    double *x_HR_fusion_largescales_all = (double*)malloc(NUM_3DSNAPS * NUM_2DSNAPS * N_HR * N_HR * sizeof(double));
    double *x_HR_fusion_smallscales_all = (double*)malloc(NUM_3DSNAPS * NUM_2DSNAPS * N_HR * N_HR * sizeof(double));
    double *x_NLM_smallscales_all = (double*)malloc(NUM_3DSNAPS * NUM_2DSNAPS * N_HR * N_HR * sizeof(double));

    /* read all snapshots */
    size_t start_ids[4] = {0, 0, 0, 0};
    size_t count_ids[4] = {NUM_3DSNAPS, NUM_2DSNAPS, N_HR, N_HR };
    read_netcdf(FILENAME_RD, "Uinterp_all", start_ids, count_ids, x_HR_fusion_largescales_all);
    read_netcdf(FILENAME_RD, "Udiff_all", start_ids, count_ids, x_HR_fusion_smallscales_all);

    double time_all_start = omp_get_wtime();

    double *x_prev_largescales = (double*)malloc(N_HR * N_HR * sizeof(double));
    double *x_after_largescales = (double*)malloc(N_HR * N_HR * sizeof(double));
    double *x_prev_smallscales = (double*)malloc(N_HR * N_HR * sizeof(double));
    double *x_after_smallscales = (double*)malloc(N_HR * N_HR * sizeof(double));
    double *x_current_largescales = (double*)malloc(N_HR * N_HR * sizeof(double));
    double *x_first_smallscales = (double*)malloc(N_HR * N_HR * sizeof(double));
    double *x_NLM_smallscales, *W; // pre-declare required for calloc

    long int grid_size = N_HR * N_HR * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE;
    int *gridpatches_y = (int*)malloc(grid_size * sizeof(int));
    int *gridpatches_z = (int*)malloc(grid_size * sizeof(int));
    int *acc_ids = (int*)malloc(ACC_FULLSIZE * ACC_FULLSIZE * sizeof(int));
    generate_grids(gridpatches_y, gridpatches_z, acc_ids);


    int snap3d_id;
    for(snap3d_id = 0; snap3d_id<NUM_3DSNAPS; snap3d_id++)
    {
        int t_offset = snap3d_id * NUM_2DSNAPS * N_HR*N_HR;

        // put first PIV
        get_onesnap(x_HR_fusion_smallscales_all, x_first_smallscales, t_offset + 0 * N_HR * N_HR, t_offset + 1 * N_HR * N_HR - 1);
        put_onesnap(x_NLM_smallscales_all, x_first_smallscales, t_offset + 0 * N_HR * N_HR, t_offset + 1 * N_HR * N_HR - 1);
        printf("\n Putting PIV snapshot at t= %i \n", 0);

        int block_id;
        //#pragma omp parallel for private(block_id)
        for(block_id = 0; block_id < NUM_BLOCKS; block_id++)
        {
            int t_prev = SCALE_FACTOR_TIME*block_id;
            int t_after = SCALE_FACTOR_TIME*(block_id+1);
            int t_est;

            get_onesnap(x_HR_fusion_largescales_all, x_prev_largescales, t_offset + t_prev * N_HR * N_HR, t_offset + (t_prev + 1) * N_HR * N_HR - 1);
            get_onesnap(x_HR_fusion_largescales_all, x_after_largescales, t_offset + t_after * N_HR * N_HR, t_offset + (t_after + 1) * N_HR * N_HR - 1);
            get_onesnap(x_HR_fusion_smallscales_all, x_prev_smallscales, t_offset + t_prev * N_HR * N_HR, t_offset + (t_prev + 1) * N_HR * N_HR - 1);
            get_onesnap(x_HR_fusion_smallscales_all, x_after_smallscales, t_offset + t_after * N_HR * N_HR, t_offset + (t_after + 1) * N_HR * N_HR - 1);

            // #pragma omp parallel for private(t, t_est, x_NLM_smallscales, W)
            for (int t = 1; t < SCALE_FACTOR_TIME; t++)
            {
                t_est = t_prev + t;
                //printf("\n Start estimating snapshot id %i bounded by 2 PIV planes at [%i, %i] \n", t_est, t_prev, t_after);
                double start_time = omp_get_wtime();
                get_onesnap(x_HR_fusion_largescales_all, x_current_largescales, t_offset + t_est * N_HR * N_HR, t_offset + (t_est + 1) * N_HR * N_HR - 1);

                //Initialize with zeros
                x_NLM_smallscales = (double*)calloc(N_HR * N_HR, sizeof(double)); // initialize with zeros using calloc
                W = (double*)calloc(N_HR * N_HR, sizeof(double)); // initialize with zeros using calloc

                // Propagation from 2 PIV planes
                NLmean(x_NLM_smallscales, W, x_current_largescales, x_prev_largescales, x_prev_smallscales, gridpatches_y, gridpatches_z, acc_ids);
                NLmean(x_NLM_smallscales, W, x_current_largescales, x_after_largescales, x_after_smallscales, gridpatches_y, gridpatches_z, acc_ids);

                // Normalize and put back
                norm_by_weight(N_HR*N_HR, x_NLM_smallscales, W);
                put_onesnap(x_NLM_smallscales_all, x_NLM_smallscales, t_offset + t_est * N_HR * N_HR, t_offset + (t_est + 1) * N_HR * N_HR - 1);
                double elapped_time = omp_get_wtime() - start_time;
                printf("\n Complete estimating snapshot id %i in 3D snapshot %i bounded by 2 PIV planes at [%i, %i] in %f seconds \n", t_est, snap3d_id, t_prev, t_after, elapped_time);
            }

            // Put last PIV of the block
            put_onesnap(x_NLM_smallscales_all, x_after_smallscales, t_offset + t_after * N_HR * N_HR, t_offset + (t_after + 1) * N_HR * N_HR - 1);
            printf("\n Putting PIV snapshot at t= %i \n", t_after);
        }
    }

    // Write to file
    write_netcdf(FILENAME_WR, "x_rec_all", start_ids, count_ids, x_NLM_smallscales_all);

    /* free memory */
    free(x_NLM_smallscales); free(W);
    free(x_current_largescales);
    free(x_first_smallscales);
    free(x_prev_largescales); free(x_after_largescales);
    free(x_prev_smallscales); free(x_after_smallscales);
    free(x_NLM_smallscales_all); free(x_HR_fusion_largescales_all); free(x_HR_fusion_smallscales_all);

    free(gridpatches_y); free(gridpatches_z); free(acc_ids);
    printf("\n FINISH ALL COMPUTATION IN %f SECONDS \n", (double)omp_get_wtime() - time_all_start);

    return 1;
}
