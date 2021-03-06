/* ------------
 * The code is adapted from the XSEDE online course Applications of Parallel Computing. 
 * The copyright belongs to all the XSEDE and the University of California Berkeley staff
 * that moderate this online course as well as the University of Toronto CSC367 staff.
 * This code is provided solely for the use of students taking the CSC367 course at 
 * the University of Toronto.
 * Copying for purposes other than this use is expressly prohibited. 
 * All forms of distribution of this code, whether as given or with 
 * any changes, are expressly prohibited. 
 * -------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"

#include <iostream>
#include <vector>

using namespace std;

//
//  tuned constants
//
#define density 0.0005
#define cutoff  0.01

//
//  benchmarking program
//
int main( int argc, char **argv )
{   
    int navg,nabsavg=0,numthreads; 
    double dmin, absmin=1.0,davg,absavg=0.0;
	
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" ); 
        printf( "-no turns off all correctness checks and particle output\n");   
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;      

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    //
    // Initialize grid and bins
    //
    double size = sqrt(n * density);
    double bin_size = cutoff;
    int num_bins = ceil(size / bin_size);

    vector<particle_t*> grid[num_bins][num_bins];

    //
    // put particles in the bins where they belong
    //
    for (int p = 0; p < n; p ++) {
        int row = particles[p].x / bin_size;
        int col = particles[p].y / bin_size;
        grid[row][col].push_back(&particles[p]);
    }

    //
    // initialize locks for each bin
    //
    omp_lock_t locks[num_bins][num_bins] ;
    for(int r = 0; r < num_bins; r ++){
        for(int c = 0; c < num_bins; c ++){
            omp_init_lock(&locks[r][c]);
        }
    }

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    #pragma omp parallel private(dmin) 
    {
    numthreads = omp_get_num_threads();
    for( int step = 0; step < 1000; step++ )
    {
        navg = 0;
        davg = 0.0;
	dmin = 1.0;
        //
        //  compute all forces
        //
        #pragma omp for reduction (+:navg) reduction(+:davg)
        for(int r = 0; r < num_bins; r ++) {
            for (int c = 0; c < num_bins; c ++) {
                for (int p = 0; p < grid[r][c].size(); p ++) { // for each particle in each bin
                    grid[r][c][p]->ax = 0;
                    grid[r][c][p]->ay = 0;
                    for (int new_r = max(0, r - 1); new_r <= min(r + 1, num_bins - 1); new_r ++) {
                        for (int new_c = max(0, c - 1); new_c <= min(c + 1, num_bins - 1); new_c ++) {
                            for (int new_p = 0; new_p < grid[new_r][new_c].size(); new_p ++) {
                                apply_force(*grid[r][c][p], *grid[new_r][new_c][new_p], &dmin, &davg, &navg);
                            }
                        }
                    }
                }
            }
        }
        
		
        //
        //  move particles
        //
        #pragma omp for
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );

        //
        // update bins
        //
        #pragma omp for
        // for each bin
        for (int r = 0; r < num_bins; r ++) {
            for (int c = 0; c < num_bins; c ++) {
                // for each particle in the bin
                for (int p = 0; p < grid[r][c].size(); p ++) {
                    omp_set_lock(&locks[r][c]);
                    int new_r = grid[r][c][p]->x / bin_size;
                    int new_c = grid[r][c][p]->y / bin_size;
                    // move the particle to a new bin if necessary
                    if (r != new_r || c != new_c) {
                        grid[r][c].erase(grid[r][c].begin() + p);
                    }
                    omp_unset_lock(&locks[r][c]);
                    
                    if (r != new_r || c != new_c) {
                        omp_set_lock(&locks[new_r][new_c]);
                        grid[new_r][new_c].push_back(grid[r][c][p]);
                        omp_unset_lock(&locks[new_r][new_c]);
                    }
                }
            }
        }
  
        if( find_option( argc, argv, "-no" ) == -1 ) 
        {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) { 
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
	  if (dmin < absmin) absmin = dmin; 
		
          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
}
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");
    
    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
