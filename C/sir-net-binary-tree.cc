/* Numerical simulations of the SIR on networks using the direct method of Gillespie.
Undirected networks are assumed.
A binary search tree is used. */

#include <iostream>
using namespace std;
#include <fstream>
#include <cstdlib> // atoi
#include <cstring>
#include <cmath> // sqrt

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
        cerr << "Usage: sir-net-binary-tree.out infile" << endl;
        exit(8);
    }

    init_genrand(time(NULL)); // seed the random number generator

    int i, j, ii, jj, ind_min, ind_max;
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

    double beta = 0.6; // infection rate
    double mu = 1.0; // recovery rate
    int st[nV]; // state of each node. 0:S, 1:I, 2:R
    int kI[nV]; // # infected neighbors. When the focal node is I or R, we set -1.
    int pt0; // index patient

    double t; // time
    int trials=1; // not multiplied by # nodes. 
    int tr;
    double ra; // random variate
    int nI; // # infected nodes
    int nR; // # recovered nodes

    for (tr=0 ; tr<trials ; tr++) {

    // initialization
    for (i=0 ; i<nV ; i++) {
        st[i] = 0; // S
        kI[i] = 0; // # infected heighbors
        lambda[i] = 0.0; // node's event rate
    }
    pt0 = 0; // index patient, i.e., ID of the single node that is initially infected
    nI = 1; // # I
    nR = 0; // # R 
    st[pt0] = 1; // I; initially infected node
    lambda[pt0] = mu; // recovery rate
    
    ind_min = (pt0==0)? 0 : accum_k[pt0-1];
    ind_max = accum_k[pt0] - 1;
    for (j = ind_min ; j <= ind_max ; j++) {
        kI[nei_list[j]] = 1; // Node nei_list[j] initially has 1 infected neighbor, which is node pt0.
        lambda[nei_list[j]] = beta;
    }

    // initialize the binary search tree    
    for (i = Np ; i < 2*pow2[log2_Np_ceil] ; i++)
        lambda[i] = 0.0;
    for (i = 0 ; i < log2_Np_ceil ; i++) { // construct binary tree
        for (j = 0 ; j < pow2[log2_Np_ceil-i] ; j++)
            lambda[start[i+1]+j/2] += lambda[start[i]+j];
    }

    if (tr >= 1)
        cout << "nan nan nan nan" << endl; // separating the results for different runs. This is useful for plotting the results by Python matplotlib

    // dynamics
    t = 0.0; // initialize time
    while (nI > 0) { // I exists 

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

        // The rates of neighbors of i may be affected
        ind_min = (i==0)? 0 : accum_k[i-1];
        ind_max = accum_k[i] - 1;

        if (st[i]==1) { // I -> R

            for (j = ind_min ; j <= ind_max ; j++) { // node nei_list[j] loses one infected neighbor, which is node i
                if (st[nei_list[j]]==0) {
                    kI[nei_list[j]]--;
                    ii = nei_list[j];
                    for (jj = 0 ; jj <= log2_Np_ceil ; jj++) { // update the binary tree including leaf node i
                        lambda[start[jj]+ii] -= beta;
                        ii /= 2;
                    }
                }
            }
            st[i] = 2; // I -> R
            ii = i;
            for (jj = 0 ; jj <= log2_Np_ceil ; jj++) { // update the binary tree including leaf node i
                lambda[start[jj]+ii] -= mu;
                ii /= 2;
            }
            nI--;
            nR++;
            
        } else if (st[i]==0) { // S -> I

            for (j = ind_min ; j <= ind_max ; j++) { // node nei_list[j] gains one infected neighbor, which is node i
                if (st[nei_list[j]]==0) {
                    kI[nei_list[j]]++;
                    ii = nei_list[j];
                    for (jj = 0 ; jj <= log2_Np_ceil ; jj++) { // update the binary tree including leaf node i
                        lambda[start[jj]+ii] += beta;
                        ii /= 2;
                    }
                }
            }
            st[i] = 1; // S -> I
            Dlambda = mu - lambda[i]; // increment in lambda[i] = mu - beta * kI[i]
            ii = i;
            for (jj = 0 ; jj <= log2_Np_ceil ; jj++) { // update the binary tree including leaf node i
                lambda[start[jj]+ii] += Dlambda;
                ii /= 2;
            }
            nI++;
        
        }

        cout << t << " " << (double)(nV-nI-nR)/nV << " " << (double)nI/nV << " " << (double)nR/nV << endl; // one state-transition event completed
     } // one trial completed

    } // all the trials completed

    free_memory(nei_list, k, accum_k); // free memory
    return 0;
}