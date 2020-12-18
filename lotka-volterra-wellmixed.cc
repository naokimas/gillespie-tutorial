/* Numericala simulations of the Lotka-Volterra model in well-mixed populations using the direct method of Gillespie. */

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

    if (argc < 1) {
        cerr << "Usage: lotka-volterra-wellmixed.out" << endl;
    exit(8);
    }

    init_genrand(time(NULL)); // seed the random number generator

    int i,j;
    double alpha = 30; // natural birth date for rabbits
    double beta = 0.1; // predation rate
    double mu = 30; // natural death rate for foxes
    int trials = 1; // # trials
    int tr;

    double t; // time
    double t_max = 50;
    int nRabbit, nFox; // # rabbits and # foxes
    int nRabbit_init = 80; // initial condition
    int nFox_init = 20; // initial condition
    double ra;
    double total_rate; // total state-transition rate 

    for (tr=0 ; tr<trials ; tr++) {

    // initialization
    nRabbit = nRabbit_init;
    nFox = nFox_init;
    t = 0.0;
    if (tr >= 1)
        cout << "nan nan nan" << endl; // separating the results for different runs. This is useful for plotting the results by Python matplotlib

    while (nRabbit > 0 && nFox > 0 && t < t_max) { // There are still two opinions coexisting.

        total_rate = alpha * nRabbit + beta * nRabbit * nFox + mu * nFox;
        t += - log((genrand_int32()+0.5)/4294967296.0) / total_rate;

    	// determine whether A->B or B->A occurs
    	ra = (genrand_int32()+0.5)/4294967296.0 * total_rate; // ra \in [0, total_rate], uniformly distributed
        if (ra < alpha * nRabbit) // natural birth of a rabbit
            nRabbit++;
        else {
            ra -= alpha * nRabbit; // now ra \in [0, beta * nRabbit * nFox + mu * nFox]
            if (ra < beta * nRabbit * nFox) { // predation
                nRabbit--;
                nFox++;
            } else { // natural death of a fox
                nFox--;
            }
        }

        cout << t << " " << nRabbit << " " << nFox << endl;
        // one state-transition event completed
    } // one trial completed    
    } // all the trials completed

    return 0;
} // end of main
