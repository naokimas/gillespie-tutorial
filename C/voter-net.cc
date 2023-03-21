/* Numericala simulations of the voter model on networks using the direct method of Gillespie.
Undirected networks are assumed. */

#include <iostream>
using namespace std;
#include <fstream>
#include <cstdlib> // atoi
#include <cstring>
#include <cmath> // log

#include "mt19937ar.c"
/* Mersenne Twister to generate random variates on [0,1].
    mt19937ar.c is available at 
    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c

    If one prefers to use the in-built pseudo-random number generator,
    (genrand_int32()+0.5)/4294967296.0 should be replaced by (double)rand()/RAND_MAX
    and
    init_genrand(time(NULL)) should be replaced by, e.g., srand(time(NULL))
    However, we recommend against the use of the in-built pseudo-random number generator. */

#include "tools-readfile.cc" // to read input file and manage memory

int main (int argc, char **argv) {

    if (argc != 2) {
        cerr << "Usage: voter-net.out infile" << endl;
        cerr << "nV = number of nodes" << endl; 
        exit(8);
    }

    init_genrand(time(NULL)); // seed the random number generator

    int i, j, vs, ve, ind_min, ind_max;
    int nV; // # nodes
    int nE; // # bidirectional edges. If the network is undirected, nE = 2 * (number of edges)
    int *nei_list; // list of neighbors for each node
    int *k; // k[i] = degree of node [i]
    int *accum_k; // cumulative degree distribution
    
    readfile(argv[1], &nV, &nE, &nei_list, &k, &accum_k);

    double beta_B_to_A = 1.0; // strength of opinion A. The strength of opinion B is assumed to be 1. The unbiased voter model corresponds to r=1.
    int st[nV]; // node's opinion. 0:A, 1:B

    double t; // time
    int trials = 3; // # trials
    int tr;
    double ra; // random variate
    double total_rate; // total state-transition rate 
    double rate[nV]; // state-transition rate for a node
    double accum_rate;
    double rate_new;
    int n_opposite_neighbors[nV]; // n_opposite_neighbors[i] = # neighbors in the opposite opinion of node i
    int nA; // # nodes in opinion A

    for (tr=0 ; tr<trials ; tr++) {

    // initialization
    for (i=0 ; i<nV/2 ; i++) st[i] = 0; // opinion A
    for (i=nV/2 ; i<nV ; i++) st[i] = 1; // opinion B
    nA=nV/2; // # opinion A
    total_rate = 0.0;
    for (i=0 ; i<nV ; i++) {
        n_opposite_neighbors[i] = 0;
        ind_min = (i==0)? 0 : accum_k[i-1];
        ind_max = accum_k[i] - 1;
        for (j = ind_min ; j <= ind_max ; j++) {
                if (st[i] != st[nei_list[j]])
                    n_opposite_neighbors[i]++;
        }
        if (st[i]==0)
            rate[i] = (double)n_opposite_neighbors[i]; // opinion A -> opinion B
        else
            rate[i] = beta_B_to_A * n_opposite_neighbors[i]; // B -> A  
        total_rate += rate[i];
    }
    t=0.0; // initialize time

    if (tr >= 1)
        cout << "nan nan" << endl; // separating the results for different runs. This is useful for plotting the results by Python matplotlib
        // dynamics
        while (total_rate > 1e-6) { // There are still two opinions coexisting.

        t += - log((genrand_int32()+0.5)/4294967296.0) / total_rate; // increment in t

    	// determine the node to be updated
    	ra = (genrand_int32()+0.5)/4294967296.0 * total_rate; // ra \in [0, total_rate], uniformly distributed
        i = 0;
        accum_rate = rate[0];
        while (accum_rate < ra) { // linear search, which is not fast
            i++;
            accum_rate += rate[i];
        }

        st[i] = 1-st[i]; // flip the opinion from 0 to 1 or from 1 to 0
        if (st[i]==0)
            nA++;
        else
            nA--;
        
        // update rate[i]
        n_opposite_neighbors[i] = k[i] - n_opposite_neighbors[i];
        if (st[i]==0)
            rate_new = (double)n_opposite_neighbors[i]; // A -> B
        else
            rate_new = beta_B_to_A * n_opposite_neighbors[i]; // B -> A
        total_rate += rate_new - rate[i];
        rate[i] = rate_new;
        
        // update neighbors' event rates
        ind_min = (i==0)? 0 : accum_k[i-1];
        ind_max = accum_k[i] - 1;
            for (j = ind_min ; j <= ind_max ; j++) {
                if (st[i] == st[nei_list[j]]) // i's new state is the same as the j'th neighbor's state
                    n_opposite_neighbors[nei_list[j]]--;
                else    
                    n_opposite_neighbors[nei_list[j]]++;
                if (st[nei_list[j]]==0)
                    rate_new = (double)n_opposite_neighbors[nei_list[j]]; // A -> B
                else
                    rate_new = beta_B_to_A * n_opposite_neighbors[nei_list[j]];
                total_rate += rate_new - rate[nei_list[j]];
                rate[nei_list[j]] = rate_new;
            }

        cout << t << " " << (double)nA/nV << endl; // one state-transition event completed
    } // one trial completed
    } // all the trials completed

    free_memory(nei_list, k, accum_k); // free memory
    return 0;
}