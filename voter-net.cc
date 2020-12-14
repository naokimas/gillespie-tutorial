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

    If one prefers to use the in-built random number generator,
    (genrand_int32()+0.5)/4294967296.0 should be replaced by (double)rand()/RAND_MAX
    and
    init_genrand(time(NULL)) should be replaced by, e.g., srand(time(NULL)) */

int compare (const void *a, const void *b) {
  return (*(int*)a-*(int*)b);
}

int main (int argc, char **argv) {

    if (argc != 2) {
        cerr << "Usage: voter-net-net.out infile" << endl;
        exit(8);
    }

    init_genrand(time(NULL)); // seed the random number generator

    int i,j,vs,ve;
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

    double r = 1.0; // strength of opinion A. The strength of opinion B is assumed to be 1. The unbiased voter model corresponds to r=1.
    int tr;
    int trials = 4; // # trials
    // If one wants to do x times for each index patient, set trials = x*n;

    double t, dt; // time
    int st[nV]; // node's opinion
    double ra; // random variate
    double total_rate; // total state-transition rate 
    double rate[nV]; // state-transition rate for a node
    double accum_rate;
    double rate_new;
    
    int n_opposite_neighbors[nV]; // n_opposite_neighbors[i] = # neighbors in the opposite opinion of node i
    int nA; // # nodes in opinion A
    int ind_min, ind_max;

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
                if (st[i] != st[E[j]])
                    n_opposite_neighbors[i]++;
        }
        if (st[i]==0)
            rate[i] =  (double)n_opposite_neighbors[i]; // opinion A -> opinion B
        else
            rate[i] = r * n_opposite_neighbors[i]; // B -> A  
        total_rate += rate[i];
    }
    t=0.0; // initialize time

    if (tr >= 1)
        cout << "nan nan" << endl; // separating the results for different runs. This is useful for plotting the results by Python matplotlib
        // dynamics
        while (total_rate > 1e-6) { // There are still two opinions coexisting.

        dt = -1.0 / total_rate * log ((genrand_int32()+0.5)/4294967296.0); // increment in t
        if ((int)((t+dt)/10000) > (int)(t/10000)) { // avoid potential accumulation of roundoff error
            total_rate = 0.0;
            for (i=0 ; i<nV ; i++)
                total_rate += rate[i];
        }
        t += dt;
    	// (double)rand()/RAND_MAX \in [0, 1], uniformly distributed


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
            rate_new = r * n_opposite_neighbors[i]; // B -> A
        total_rate += rate_new - rate[i];
        rate[i] = rate_new;
        
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
                    rate_new = r * n_opposite_neighbors[E[j]];
                total_rate += rate_new - rate[E[j]];
                rate[E[j]] = rate_new;
            }

        	cout << t << " " << (double)nA/nV << endl;
    } // one state-transition event completed
    
    } // all the trials completed

  free(E); // free memory
  return 0;
}