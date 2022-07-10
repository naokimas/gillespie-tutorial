/* Numericala simulations of the voter model on networks using the direct method of Gillespie.
Undirected networks are assumed.
A binary search tree is used. */

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
        cerr << "Usage: voter-net-binary-tree.out infile" << endl;
        exit(8);
    }

    init_genrand(time(NULL)); // seed the random number generator

    int i, j, ii, jj, vs, ve, ind_min, ind_max;
    int nV; // # nodes
    int nE; // # bidirectional edges. If the network is undirected, nE = 2 * (number of edges)
    int *nei_list; // list of neighbors for each node
    int *k; // k[i] = degree of node [i]
    int *accum_k; // cumulative degree distribution
    
    readfile(argv[1], &nV, &nE, &nei_list, &k, &accum_k);

    int Np = nV; // number of channels 
    // prepare to construct a binary tree
    int log2_Np_ceil = (int)(log(Np-1e-6)/log(2)) + 1;
    // 2^{log2_Np_ceil-1} < Np <= 2^{log2_Np_ceil}
    int pow2[log2_Np_ceil+1]; // power[i] = 2^i.
    pow2[0]=1;
    for (i = 0 ; i < log2_Np_ceil ; i++)
        pow2[i+1] = 2 * pow2[i];

    // event rates and the binary tree
    double lambda[2*pow2[log2_Np_ceil]]; // rate of Poisson processes and associated binary tree
    double Dlambda; // working var for updating the binary tree
    int start[log2_Np_ceil+1]; // start[i] = j if the i-th layer of the binary tree starts with lambda[j]
    start[0] = 0;
    for (i = 0 ; i < log2_Np_ceil ; i++)
        start[i+1] = start[i] + pow2[log2_Np_ceil-i];

    /* --- layer 0 ---
        lambda[0], ..., lambda[Np-1] contain original lambda.
        lambda[Np], ..., lambda[2^{log2_Np_ceil}-1] are set to 0 (dummy).
        These elements of lambda[] are leaves of the binary tree and define layer 0, so start[0] = 0.

        --- layer 1 ---  
        lambda[2^{log2_Np_ceil}], .., lambda[2^{log2_Np_ceil}+2^{log2_Np_ceil-1}-1] contain the sum of two adjacent lambda[i]'s in layer 0,
        where 0 <= i < 2^{log2_Np_ceil}.
        Therefore, start[1] = 2^{log2_Np_ceil}.

        --- layer 2 ---
        lambda[2^{log2_Np_ceil}+2^{log2_Np_ceil-1}], ..., lambda[2^{log2_Np_ceil}+2^{log2_Np_ceil-1}+2^{log2_Np_ceil-2}-1] contain the sum of two adjacent lambda[i]'s in layer 1,
        where 2^{log2_Np_ceil} <= i < 2^{log2_Np_ceil}+2^{log2_Np_ceil-1}).
        Therefore, start[2] = 2^{log2_Np_ceil}+2^{log2_Np_ceil-1}.
        Similar for layers 3, 4, ... */
        
    double beta_B_to_A = 1.0; // strength of opinion A. The strength of opinion B is assumed to be 1. The unbiased voter model corresponds to r=1.
    int st[nV]; // node's opinion. 0:A, 1:B

    double t; // time
    int trials = 1; // # trials
    int tr;
    double ra; // random variate
    double rate_new;
    int n_opposite_neighbors[nV]; // n_opposite_neighbors[i] = # neighbors in the opposite opinion of node i
    int nA; // # nodes in opinion A

    for (tr=0 ; tr<trials ; tr++) {

    // initialization
    for (i=0 ; i<nV/2 ; i++) st[i] = 0; // opinion A
    for (i=nV/2 ; i<nV ; i++) st[i] = 1; // opinion B
    nA=nV/2; // # opinion A
    for (i=0 ; i<nV ; i++) {
        n_opposite_neighbors[i] = 0;
        ind_min = (i==0)? 0 : accum_k[i-1];
        ind_max = accum_k[i] - 1;
        for (j = ind_min ; j <= ind_max ; j++) {
                if (st[i] != st[nei_list[j]])
                    n_opposite_neighbors[i]++;
        }
        if (st[i]==0)
            lambda[i] = (double)n_opposite_neighbors[i]; // opinion A -> opinion B
        else
            lambda[i] = beta_B_to_A * n_opposite_neighbors[i]; // B -> A  
    }

    // initialize the binary search tree    
    for (i = Np ; i < 2*pow2[log2_Np_ceil] ; i++)
        lambda[i] = 0.0;
    for (i = 0 ; i < log2_Np_ceil ; i++) { // construct binary tree
        for (j = 0 ; j < pow2[log2_Np_ceil-i] ; j++)
            lambda[start[i+1]+j/2] += lambda[start[i]+j];
    }
    
    if (tr >= 1)
        cout << "nan nan" << endl; // separating the results for different runs. This is useful for plotting the results by Python matplotlib

    t=0.0; // initialize time
    // dynamics
    while (lambda[2*pow2[log2_Np_ceil]-2] > 1e-6) { // There are still two opinions coexisting.

        t += - log((genrand_int32()+0.5)/4294967296.0) / lambda[2*pow2[log2_Np_ceil]-2];

    	// determine the node to be updated
        ra = (genrand_int32()+0.5)/4294967296.0 * lambda[2*pow2[log2_Np_ceil]-2]; // ra is between 0 and lambda[2*pow2[log2_Np_ceil]-2]
        i = 0;
        for (jj = log2_Np_ceil-1 ; jj >= 0 ; jj--) { // fast search using binary tree
            i *= 2;
            if (ra > lambda[start[jj]+i]) {
                ra -= lambda[start[jj]+i]; 
                i++;
            }
        } // node i will change the state

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
        Dlambda = rate_new - lambda[i]; // increment in lambda[i]
        ii = i;
        for (jj = 0 ; jj <= log2_Np_ceil ; jj++) { // update the binary tree including leaf node i
            lambda[start[jj]+ii] += Dlambda;
            ii /= 2;
        }

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
                Dlambda = rate_new - lambda[nei_list[j]]; // increment in lambda[nei_list[j]]
                ii = nei_list[j];
                for (jj = 0 ; jj <= log2_Np_ceil ; jj++) { // update the binary tree including leaf node nei_list[j]
                    lambda[start[jj]+ii] += Dlambda; // update the binary tree
                    ii /= 2;
                }
            }

        cout << t << " " << (double)nA/nV << endl; // one state-transition event completed
    } // one trial completed
    } // all the trials completed

    free_memory(nei_list, k, accum_k); // free memory
    return 0;
}