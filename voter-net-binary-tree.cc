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


int main (int argc, char **argv) {

    if (argc != 2) {
        cerr << "Usage: voter-net-net.out infile" << endl;
        exit(8);
    }

    init_genrand(time(NULL)); // seed the random number generator

    int i, j, ii, jj, vs, ve;
    int nV; // # nodes
    int nE; // # bidirectional edges. If the network is undirected, nE = 2 * (number of edges)
    char dummy[8];

    ifstream fin(argv[1]); // edge list
    if (!fin) {
        cerr << "voter-net.out: cannot open input edge file" << endl;
        exit(8);
    }

    // *** I SHOULD ASSUME A MORE POPULAR DATA FORMAT AND ADJUST THE WRAPPER AROUND HERE.

    fin >> dummy >> nV >> nE; // The first line of the file has '#', the number of nodes, and the number of edges, space separated
    cerr << nV << " nodes, " << nE << " edges" << endl; // undirected net assumed
    int k[nV]; // k[i] = degree of node i
    for (i=0 ; i<nV ; i++)
        k[i] = 0; // initialization
    int *E; // list of edges
    E = (int *)malloc(4*nE*sizeof(int)); // allocate memory
//    cerr << "hello" << endl;
    for (i=0 ; i<nE ; i++) { // read all edges
        fin >> vs >> ve >> dummy; // Each line of the input file has the index of the two nodes connected by an edge in the 1st and 2nd columns, and edge weight, which we do not use in this code, in the 3rd column. 
    // Note that in this code the node index starts from 0, not from 1.
    // And in the input file the node index starts from 1.
        vs--; ve--; // node index corrected here. Now the node index starts from 0.
        E[4*i+3]=E[4*i]=vs; // i-th edge
        E[4*i+2]=E[4*i+1]=ve;
        k[vs]++; k[ve]++; // undirected network assumed
    }
    fin.close();

    int accum_k[nV]; // accum_k[i] = sum of degrees of nodes from 0 to i
    accum_k[0] = k[0];
    for (i=1 ; i<nV ; i++) accum_k[i] = accum_k[i-1] + k[i];

    qsort(E,2*nE,2*sizeof(int),compare);
    for (i=0 ; i<2*nE ; i++) 
        E[i] = E[2*i+1]; // reuse E[]
        // Now E[x] where x = accum_k[i-1], ..., accum_k[i]-1 contains the neighbors of node i, with the convention accum_k[-1] = 0.
        // Note that accum_k[nV-1] = nV-1
    
    /* loading network data and preprocessing done */

    int Np = nV; // number of channels 
    // preparation to construct a binary tree
    int log2_Np_upper = (int)(log(Np-1e-6)/log(2)) + 1;
    // 2^{log2_Np_upper-1} < Np <= 2^{log2_Np_upper}
    int pow2[log2_Np_upper+1]; // power[i] = 2^i.
    pow2[0]=1;
    for (i = 0 ; i < log2_Np_upper ; i++)
        pow2[i+1] = 2 * pow2[i];

    // event rates and the binary tree
    double lambda[2*pow2[log2_Np_upper]]; // rate of Poisson processes and associated binary tree
    double Dlambda; // working var for updating the binary tree
    int start[log2_Np_upper+1]; // start[i] = j if the i-th layer of the binary tree starts with lambda[j]
    start[0] = 0;
    for (i = 0 ; i < log2_Np_upper ; i++)
        start[i+1] = start[i] + pow2[log2_Np_upper-i];

    /* --- layer 0 ---
        lambda[0], ..., lambda[Np-1] contain original lambda.
        lambda[Np], ..., lambda[2^{log2_Np_upper}-1] are set to 0 (dummy).
        These elements of lambda[] are leaves of the binary tree and define layer 0, so start[0] = 0.

        --- layer 1 ---  
        lambda[2^{log2_Np_upper}], .., lambda[2^{log2_Np_upper}+2^{log2_Np_upper-1}-1] contain the sum of two adjacent lambda[i]'s in layer 0,
        where 0 <= i < 2^{log2_Np_upper}.
        Therefore, start[1] = 2^{log2_Np_upper}.

        --- layer 2 ---
        lambda[2^{log2_Np_upper}+2^{log2_Np_upper-1}], ..., lambda[2^{log2_Np_upper}+2^{log2_Np_upper-1}+2^{log2_Np_upper-2}-1] contain the sum of two adjacent lambda[i]'s in layer 1,
        where 2^{log2_Np_upper} <= i < 2^{log2_Np_upper}+2^{log2_Np_upper-1}).
        Therefore, start[2] = 2^{log2_Np_upper}+2^{log2_Np_upper-1}.
        Similar for layers 3, 4, ... */
        
    double beta_B_to_A = 1.0; // strength of opinion A. The strength of opinion B is assumed to be 1. The unbiased voter model corresponds to r=1.
    int trials = 1; // # trials
    int tr;

    double t; // time
    int st[nV]; // node's opinion
    double ra; // random variate
    double rate_new;
    
    int n_opposite_neighbors[nV]; // n_opposite_neighbors[i] = # neighbors in the opposite opinion of node i
    int nA; // # nodes in opinion A
    int ind_min, ind_max;

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
                if (st[i] != st[E[j]])
                    n_opposite_neighbors[i]++;
        }
        if (st[i]==0)
            lambda[i] = (double)n_opposite_neighbors[i]; // opinion A -> opinion B
        else
            lambda[i] = beta_B_to_A * n_opposite_neighbors[i]; // B -> A  
    }

    // initialize the binary search tree    
    for (i = Np ; i < 2*pow2[log2_Np_upper] ; i++)
        lambda[i] = 0.0;
    for (i = 0 ; i < log2_Np_upper ; i++) { // construct binary tree
        for (j = 0 ; j < pow2[log2_Np_upper-i] ; j++)
            lambda[start[i+1]+j/2] += lambda[start[i]+j];
    }
    
    t=0.0; // initialize time

    if (tr >= 1)
        cout << "nan nan" << endl; // separating the results for different runs. This is useful for plotting the results by Python matplotlib

    // dynamics
    while (lambda[2*pow2[log2_Np_upper]-2] > 1e-6) { // There are still two opinions coexisting.

        t += - log((genrand_int32()+0.5)/4294967296.0) / lambda[2*pow2[log2_Np_upper]-2]; // increment in t

    	// determine the node to be updated
        ra = (genrand_int32()+0.5)/4294967296.0 * lambda[2*pow2[log2_Np_upper]-2]; // ra is between 0 and lambda[2*pow2[log2_Np_upper]-2]
        i = 0;
        for (jj = log2_Np_upper-1 ; jj >= 0 ; jj--) { // fast search using binary tree
            i *= 2;
            if (ra > lambda[start[jj]+i]) {
                ra -= lambda[start[jj]+i]; 
                i++;
            }
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
        Dlambda = rate_new - lambda[i]; // increment in lambda[i]
        ii = i;
        for (jj = 0 ; jj <= log2_Np_upper ; jj++) { // update the binary tree including leaf node i
            lambda[start[jj]+ii] += Dlambda;
            ii /= 2;
        }

        // update neighbors' event rates
        ind_min = (i==0)? 0 : accum_k[i-1];
        ind_max = accum_k[i] - 1;
            for (j = ind_min ; j <= ind_max ; j++) {
                if (st[i] == st[E[j]]) // i's new state is the same as the j'th neighbor's state
                    n_opposite_neighbors[E[j]]--;
                else    
                    n_opposite_neighbors[E[j]]++;
                if (st[E[j]]==0)
                    rate_new = (double)n_opposite_neighbors[E[j]]; // A -> B
                else
                    rate_new = beta_B_to_A * n_opposite_neighbors[E[j]];
                Dlambda = rate_new - lambda[E[j]]; // increment in lambda[E[j]]
                ii = E[j];
                for (jj = 0 ; jj <= log2_Np_upper ; jj++) { // update the binary tree including leaf node E[j]
                    lambda[start[jj]+ii] += Dlambda; // update the binary tree
                    ii /= 2;
                }
            }

        cout << t << " " << (double)nA/nV << endl; // one state-transition event completed
    } // one trial completed
    } // all the trials completed

    free(E); // free memory
    return 0;
}