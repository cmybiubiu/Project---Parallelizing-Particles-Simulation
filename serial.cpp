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

#include <iostream>
#include <vector>
#include <algorithm>
#include <unistd.h>


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
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
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

    vector<int> grid[num_bins][num_bins];

    //
    // put particles in the bins they belong
    //
    for (int p = 0; p < n; p ++) {
        int row = particles[p].x / bin_size;
        int col = particles[p].y / bin_size;
        grid[row][col].push_back(p);
    }


    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;
        //
        //  compute forces
        //
        for (int p = 0; p < n; p ++) {
            particles[p].ax = particles[p].ay = 0; // init acceleration
            // find where the particle is
            int row = particles[p].x / bin_size;
            int col = particles[p].y / bin_size;

            for (int r = max(0, row - 1); r <= min(row + 1, num_bins - 1); r ++) {
                for (int c = max(0, col - 1); c <= min(col + 1, num_bins - 1); c ++) {
                    for (int i = 0; i < grid[r][c].size(); i ++) {
                        int p_i = grid[r][c][i];
                        apply_force(particles[p], particles[p_i], &dmin, &davg, &navg);
                    }
                }
            }
        }

        //
        //  move particles
        //
        for( int i = 0; i < n; i++ )
            move( particles[i] );

        //
        // update grid / bins
        //
        for (int r = 0; r < num_bins; r ++) {
            for (int c = 0; c < num_bins; c ++) {
                // for each particle in the bin
                for (int i = 0; i < grid[r][c].size(); i ++) {
                    int p_idx = grid[r][c][i];
                    int new_r = particles[p_idx].x / bin_size;
                    int new_c = particles[p_idx].y / bin_size;
                    if (r != new_r || c != new_c) {
                        grid[new_r][new_c].push_back(p_idx);
                        grid[r][c].erase(grid[r][c].begin() + i);
                    }
                }
            }
        }

        if( find_option( argc, argv, "-no" ) == -1 )
        {
            //
            // Computing statistical data
            //
            if (navg) {
                absavg +=  davg/navg;
                nabsavg++;
            }
            if (dmin < absmin) absmin = dmin;

            //
            //  save if necessary
            //
            if( fsave && (step%SAVEFREQ) == 0 )
                save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

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
        fprintf(fsum,"%d %g\n",n,simulation_time);

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
