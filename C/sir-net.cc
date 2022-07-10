/* Numerical simulations of the SIR on networks using the direct method of Gillespie.
Undirected networks are assumed.

    Copyright (C) 2020, Naoki Masuda, All rights reserved.
    Please cite the following paper when you use this code.
    Naoki Masuda, Christian L. Vestergaard. Gillespie algorithms for stochastic multiagent dynamics in populations and networks. arXiv:xxx, Cambridge Elements, in review (2020).
*/

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
        cerr << "Usage: sir-net.out infile" << endl;
        exit(8);
    }

    init_genrand(time(NULL)); // seed the random number generator

    int i, j, ind_min, ind_max;
    int nV; // # nodes
    int nE; // # bidirectional edges. If the network is undirected, nE = 2 * (number of edges)

    int *nei_list; // list of neighbors for each node
    int *k; // k[i] = degree of node [i]
    int *accum_k; // cumulative degree distribution
    
    readfile(argv[1], &nV, &nE, &nei_list, &k, &accum_k);

    double beta = 0.6; // infection rate
    double mu = 1.0; // recovery rate
    int st[nV]; // state of each node. 0:S, 1:I, 2:R
    int kI[nV]; // # infected neighbors. When the focal node is I or R, we set -1.
    int pt0; // index patient

    double t; // time
    int trials=1; // not multiplied by # nodes. 
    int tr;
    double ra; // random variate
    double rate[nV]; // event rate of each node
    double total_rate; // total state-transition rate 
    double accum_rate;
    int nI; // # infected nodes
    int nR; // # recovered nodes

    for (tr=0 ; tr<trials ; tr++) {

    // initialization
    for (i=0 ; i<nV ; i++) {
        st[i] = 0; // S
        kI[i] = 0; // # infected heighbors
        rate[i] = 0.0; // node's event rate
    }
    pt0 = 0; // index patient, i.e., ID of the single node that is initially infected
    nI = 1; // # I
    nR = 0; // # R 
    st[pt0] = 1; // I; initially infected node
    rate[pt0] = mu; // recovery rate
    
    ind_min = (pt0==0)? 0 : accum_k[pt0-1];
    ind_max = accum_k[pt0] - 1;
    for (j = ind_min ; j <= ind_max ; j++) {
        kI[nei_list[j]] = 1; // Node nei_list[j] initially has 1 infected neighbor, which is node pt0.
        rate[nei_list[j]] = beta;
    }

    total_rate = beta*k[pt0] + mu;

    if (tr >= 1)
        cout << "nan nan nan nan" << endl; // separating the results for different runs. This is useful for plotting the results by Python matplotlib

    // dynamics
    t = 0.0; // initialize time
    while (nI > 0) { // I exists 

        t += - log((genrand_int32()+0.5)/4294967296.0) / total_rate; 

    	// determine the node to be updated
    	ra = (genrand_int32()+0.5)/4294967296.0 * total_rate; // ra is uniformly distributed on (0, total_rate)
        i = 0;
        accum_rate = rate[0];
        while (accum_rate < ra) { // linear search, which is not fast
            i++;
            accum_rate += rate[i];
        } // Node i will change the state

        // The rates of neighbors of i may be affected
        ind_min = (i==0)? 0 : accum_k[i-1];
        ind_max = accum_k[i] - 1;

        if (st[i]==1) { // I -> R

            for (j = ind_min ; j <= ind_max ; j++) { // node nei_list[j] loses one infected neighbor, which is node i
                if (st[nei_list[j]]==0) {
                    kI[nei_list[j]]--;
                    rate[nei_list[j]] = beta * kI[nei_list[j]];
                    total_rate -= beta;
                }
            }
            st[i] = 2; // I -> R
            rate[i] = 0.0;
            total_rate -= mu;
            nI--;
            nR++;
            
        } else if (st[i]==0) { // S -> I

            for (j = ind_min ; j <= ind_max ; j++) { // node nei_list[j] gains one infected neighbor, which is node i
                if (st[nei_list[j]]==0) {
                    kI[nei_list[j]]++;
                    rate[nei_list[j]] = beta * kI[nei_list[j]];
                    total_rate += beta;
                }
            }
            st[i] = 1; // S -> I
            rate[i] = mu; // recovery rate
            total_rate += mu - beta * kI[i];
            nI++;
        
        }

        cout << t << " " << (double)(nV-nI-nR)/nV << " " << (double)nI/nV << " " << (double)nR/nV << endl; // one state-transition event completed
     } // one trial completed

    } // all the trials completed

    free_memory(nei_list, k, accum_k); // free memory
    return 0;
}