/* Numericala simulations of the voter model in well-mixed populations using the direct method of Gillespie. */

#include <iostream>
using namespace std;
#include <cstdlib> // atoi
#include <cstring>
#include <ctime>
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

    if (argc < 2) {
        cerr << "Usage: voter-wellmixed.out nV" << endl;
    exit(8);
    }

    init_genrand(time(NULL)); // seed the random number generator

    int nV = atoi(argv[1]); // # nodes
    double beta_B_to_A = 1.0; // strength of opinion A. The strength of opinion B is assumed to be 1. The unbiased voter model corresponds to r=1.

    double t; // time
    int trials = 4; // # trials
    int tr;
    int nA; // # nodes in opinion A
    double ra; // random variate
    double total_rate; // total state-transition rate 

    for (tr=0 ; tr<trials ; tr++) {

    // initialization
    nA=nV/2; // # opinion A
    t = 0.0;
    if (tr >= 1)
        cout << "nan nan" << endl; // separating the results for different runs. This is useful for plotting the results by Python matplotlib

    while (nA > 0 && nA < nV) { // There are still two opinions coexisting.

        total_rate = (1 + beta_B_to_A) * nA * (nV-nA); // rate(B->A) = beta_B_to_A * (nV-nA) * nA; rate(A->B) = 1 * nA * (nV-nA)
        t += - log((genrand_int32()+0.5)/4294967296.0) / total_rate;

    	// determine whether A->B or B->A occurs
    	ra = (genrand_int32()+0.5)/4294967296.0 * total_rate; // ra is uniformly distributed on (0, total_rate)
        if (ra < nA * (nV-nA)) // A -> B
            nA--;
        else // B -> A
            nA++;

        cout << t << " " << (double)nA/nV << endl; // one state-transition event completed
    } // one trial completed
    } // all the trials completed

    return 0;
} // end of main
