/* Numericala simulations of the SIR model in well-mixed populations using the direct method of Gillespie. */

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
        cerr << "Usage: SIR-wellmixed.out nV" << endl;
        cerr << "nV = number of nodes" << endl; 
    exit(8);
    }

    init_genrand(time(NULL)); // seed the random number generator

    int nV = atoi(argv[1]); // # nodes
    double beta = (double)1.2/nV; // infection rate
    double mu = 1.0; // recovery rate
    cerr << "infection rate = " << beta << endl;

    double t; // time
    int trials = 1; // # trials
    int tr;
    double ra; // random variate
    double total_rate; // total state-transition rate 
    int nI; // # infected nodes
    int nR; // # recovered nodes

    for (tr=0 ; tr<trials ; tr++) {

    // initialization
    nI = 1; // # infected nodes
    nR = 0; // # recovered nodes
    
    t = 0.0;
    if (tr >= 1)
        cout << "nan nan nan nan" << endl; // separating the results for different runs. This is useful for plotting the results by Python matplotlib

    // dynamics
    t = 0.0; // initialize time
    while (nI > 0) { // I exists 

        total_rate = beta * (nV-nI-nR)*nI + mu * nI; // nV - nI - nR = number of susceptible nodes
        t += - log((genrand_int32()+0.5)/4294967296.0) / total_rate;
    	ra = (genrand_int32()+0.5)/4294967296.0 * total_rate; // ra \in [0, total_rate], uniformly distributed
        if (ra < mu * nI) { // I -> R
            nI--;
            nR++;
        } else { // S -> I
            nI++;
        }
        cout << t << " " << (double)(nV-nI-nR)/nV << " " << (double)nI/nV << " " << (double)nR/nV << endl; // one state-transition event completed
    } // one trial completed
    } // all the trials completed

    return 0;
} // end of main
