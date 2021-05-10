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

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"

#include <vector>
#include <cmath>
#include <signal.h>
#include <unistd.h>

using namespace std;

//
//  tuned constants
//
#define density 0.0005
#define cutoff  0.01

double bin_size;
int num_bins;


//
// Put the particle in an appropriate bin in the grid
//
void put_particle_in_bin(particle_t& particle, vector<vector<particle_t>>& grid)
{
    int row = particle.x / bin_size;
    int col = particle.y / bin_size;
    grid[row * num_bins + col].push_back(particle);
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 
 
    //
    //  process command line parameters
    //
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
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    if( rank == 0 )
        init_particles( n, particles );
     MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);

    //
    // Initialize grid and bins
    //
    double size = sqrt(n * density);
    bin_size = cutoff;
    num_bins = ceil(size / bin_size);

    vector<vector<particle_t>> grid;
    grid.resize(num_bins * num_bins);

    //
    // Put particles in the bins where they belong
    //
    for (int p = 0; p < n; p ++) {
        int row = particles[p].x / bin_size;
        int col = particles[p].y / bin_size;
        grid[row * num_bins + col].push_back(particles[p]);
    }

    //
    // Assign work to each prosessor
    //
    int rows_per_proc = num_bins / n_proc;
    // Although each proc has a full sized grid, it is only in charge of the rows of bins
    // between my_row_start and my_row_end
    int my_row_start = rank * rows_per_proc;
    int my_row_end = (rank + 1) * rows_per_proc;
    if (rank == n_proc - 1) {
        my_row_end = num_bins;
    }
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;
        
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( find_option( argc, argv, "-no" ) == -1 )
          if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
        
        //
        //  compute all forces
        //
        for (int r = my_row_start; r < my_row_end; r ++) { // Iterate through each bin
            for (int c = 0; c < num_bins; c ++) {
                vector<particle_t>& bin = grid[r * num_bins + c];
                for (int p = 0; p < bin.size(); p ++) { // Iterate through the particles in the bin
                    bin[p].ax = 0;
                    bin[p].ay = 0;
                    // Iterate through each neighbor bins, including itself
                    for (int new_r = max(0, r - 1); new_r <= min(r + 1, num_bins - 1); new_r ++) {
                        for (int new_c = max(0, c - 1); new_c <= min(c + 1, num_bins - 1); new_c ++) {
                            vector<particle_t>& neighbor_bin = grid[new_r * num_bins + new_c];
                            for (int new_p = 0; new_p < neighbor_bin.size(); new_p ++) {
                                apply_force(bin[p], neighbor_bin[new_p], &dmin, &davg, &navg);
                            }
                        }
                    }
                }
            }
        }
        
     
        if( find_option( argc, argv, "-no" ) == -1 )
        {
          
          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

 
          if (rank == 0){
            //
            // Computing statistical data
            //
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin) absmin = rdmin;
          }
        }

        //
        //  move particles, if particles are going to new bins, take them out
        //
        vector<particle_t> move_remote; // These particles will need to be sent to another proc
        vector<particle_t> move_local; // These particles are going to new bins within the same proc
        for (int r = my_row_start; r < my_row_end; r ++) {
            for (int c = 0; c < num_bins; c ++) {
                vector<particle_t>& bin = grid[r * num_bins + c];
                int bin_tail = bin.size();
                int p = 0;
                while (p < bin_tail) {
                    move(bin[p]);
                    int new_r = bin[p].x / bin_size;
                    int new_c = bin[p].y / bin_size;
                    if (my_row_start <= new_r && new_r < my_row_end) { // This particle is not going to another proc
                        if (r != new_r || c != new_c) { // This particle is going to a new bin
                            move_local.push_back(bin[p]);
                            bin_tail -= 1;
                            bin[p] = bin[bin_tail];
                        } else { // This particle stays in the current bin
                            p += 1;
                        }
                    } else { // This particle is going to another proc
                        move_remote.push_back(bin[p]);
                        bin_tail -= 1;
                        bin[p] = bin[bin_tail];
                    }
                }
                bin.resize(p);
            }
        }

        //
        // Move the particles locally if necessary
        //
        for (int p = 0; p < move_local.size(); p ++) {
            put_particle_in_bin(move_local[p], grid);
        }

        //
        // Prepare to update ghost region
        //

        // For updating upper ghost region
        if (rank != 0) {
            for (int r = my_row_start - 1, c = 0; c < num_bins; ++c) {
                // Clear its upper ghost region, will receive new ghost region data from preceding proc
                grid[r * num_bins + c].clear();
            }
            for (int r = my_row_start, c = 0; c < num_bins; ++c) {
                // Will send its upper row to its preceiding proc, so that the preceding proc
                // could update its bottom ghost region
                vector<particle_t>& bin = grid[r * num_bins + c];
                move_remote.insert(move_remote.end(), bin.begin(), bin.end());
                bin.clear();
            }
        }
        // For updating lower ghost region
        if (rank != n_proc - 1) {
            for (int r = my_row_end, c = 0; c < num_bins; ++c) {
                // Clear its lower ghost region, will receive new ghost region data from succeeding proc
                grid[r * num_bins + c].clear();
            }
            for (int r = my_row_end - 1, c = 0; c < num_bins; ++c) {
                // Will send its bottom row to its succeeding proc, so that the preceding proc
                // could update its upper ghost region
                vector<particle_t>& bin = grid[r * num_bins + c];
                move_remote.insert(move_remote.end(), bin.begin(), bin.end());
                bin.clear();
            }
        }

        //
        // Update ghost region
        //
        vector<particle_t> incoming_particles;
        int send_count = move_remote.size();
        int recv_counts[n_proc];

        // Send send_count to proc 0
        MPI_Gather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // Now proc 0 knows how much data it will receive

        int displs[n_proc];
        int total_count = 0;

        if (rank == 0) { // proc 0 calculates displacement, and total # of particles it will receive
            displs[0] = 0;
            for (int i = 1; i < n_proc; i ++) {
                displs[i] = displs[i - 1] + recv_counts[i - 1];
            }
            total_count = recv_counts[n_proc - 1] + displs[n_proc - 1];
            incoming_particles.resize(total_count);
        }

        // Send data (all particles that need to go to another proc) to proc 0
        MPI_Gatherv(move_remote.data(), send_count, PARTICLE, incoming_particles.data(), recv_counts, displs, PARTICLE, 0, MPI_COMM_WORLD);

        // Data to send out to each proc
        vector<vector<particle_t>> scatter_particles;
        scatter_particles.resize(n_proc);

        if (rank == 0) { // Proc 0 prepares data to send out to each proc
            for (int i = 0; i < incoming_particles.size(); i ++) {
                int row = incoming_particles[i].x / bin_size;
                // Which proc this particle belongs to
                int which_proc = min(row / rows_per_proc, n_proc - 1);
                scatter_particles[which_proc].push_back(incoming_particles[i]);
                // Which row the particle is at in its prosessor
                int row_proc = row % rows_per_proc;
                if (row_proc == 0 && which_proc != 0) {
                    // This particle is on the top row of some proc, send it to the proc's previous proc as ghost region
                    scatter_particles[which_proc - 1].push_back(incoming_particles[i]);
                }
                if (row_proc == rows_per_proc - 1 && which_proc != n_proc - 1) {
                    // This particle is on the buttom row of some proc, send it to the proc's secceed proc as ghost region
                    scatter_particles[which_proc + 1].push_back(incoming_particles[i]);
                }
            }
            // Update recv_counts (going to be send count when sending data to each proc)
            for (int i = 0; i < n_proc; i ++) {
                recv_counts[i] = scatter_particles[i].size();
            }
            // Update displacement
            displs[0] = 0;
            for (int i = 1; i < n_proc; i ++) {
                displs[i] = displs[i - 1] + recv_counts[i - 1];
            }
        }

        // Let each proc know how much data it will revceive
        send_count = 0;
        MPI_Scatter(recv_counts, 1, MPI_INT, &send_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // The received particels will be stored here
        vector<particle_t> recv_particles;
        recv_particles.resize(send_count);

        // Flatten scatter_particles
        vector<particle_t> scatter_particles_flat;
        for (int i = 0; i < scatter_particles.size(); i ++) {
            scatter_particles_flat.insert(scatter_particles_flat.end(), scatter_particles[i].begin(), scatter_particles[i].end());
        }

        // Send particles to their new proc
        MPI_Scatterv(scatter_particles_flat.data(), recv_counts, displs, PARTICLE, recv_particles.data(), send_count, PARTICLE, 0, MPI_COMM_WORLD);

        // Put particles in the bins they belong. The bin could be a normal bin, or it could be a ghost bin
        for(int i = 0; i < send_count; ++i) {
            particle_t &p = recv_particles[i];
            put_particle_in_bin(p, grid);
        }
    }
    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
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
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
  
    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
