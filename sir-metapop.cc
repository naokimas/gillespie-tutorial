/* Numerical simulations of the SIR model on the metapopulation model network using the direct method of Gillespie. */

// The recovery rate is set to 1.

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
        cerr << "Usage: sir-metapop.out infile" << endl;
        exit(8);
    }

    init_genrand(time(NULL)); // seed the random number generator

    int i, j, vs, ve, upper, lower;
    int nV; // # nodes
    int nE; // # bidirectional edges. If the network is undirected, nE = 2 * (number of edges)

    int *nei_list; // list of neighbors for each node
    int *k; // k[i] = degree of node [i]
    int *accum_k; // cumulative degree distribution
    
    readfile(argv[1], &nV, &nE, &nei_list, &k, &accum_k);

    double beta = 0.04; // infection rate
    double mu = 1.0; // recovery rate
    double D= 5.0; // diffusion rate
    double rho=50; // population density

    double t;
    int trials = 1; // # trials
    int tr;

    double ra; // random variate
    int rai; // integer random variate 
    double total_rate; // total state-transition rate 
    int nP, nPtmp; // # individuals
    int nS[nV]; // # susceptible individuals in a subpopulation
    int nI[nV]; // # infected individuals in a subpopulation
    int nIsum; // # infected individuals
    int nRsum; // # recovered individuals
    int accum_n_notR[nV], accum_nI[nV], accum_SIpair[nV], inc_SIpair;

    for (tr=0 ; tr<trials ; tr++) {

    // initialize nP[i] such that nP[i] propto k[i]
    nP=(int)(nV*rho);
    nPtmp = nP;
    for (i=0 ; i<nV ; i++) {
        nS[i]=nPtmp*k[i]/(2*nE); // integer
        nPtmp -= nS[i];
    }
    //      cerr << "N remaining = " << nPtmp << endl;
    if (nPtmp<0) {
        cerr << "nPtmp<0 invalid" << endl;
        exit(8);
    }
    for (i=0 ; i<nPtmp ; i++) { 
        /* assign each remaining individual to a metapopulation using the bisection method */
        lower=0;
        upper=nV-1;
        vs=upper;
        rai = genrand_int32()%accum_k[nV-1];
        while (lower != vs || upper != vs) {
            vs = (lower+upper)/2;
            if (rai >= accum_k[vs]) lower = vs+1;
            else upper=vs;
        }
        nS[vs]++;
    }
    // all individuals are allocated to subpopulations

    // cumulative distribution of the number of individuals excluding recovered individuals
    accum_n_notR[0]=nS[0];
    for (i=1 ; i<nV ; i++) accum_n_notR[i] = accum_n_notR[i-1]+nS[i];
    // For now, accum_n_notR[nV-1]  = nP

    // start from a single infected individual
    for (i=0 ; i<nV ; i++) nI[i]=0;
    rai = genrand_int32() % accum_n_notR[nV-1];
    lower=0; // determine the subpopulation containing the index patient by the bisection method
    upper=nV-1;
    vs=upper;
    while (lower != vs || upper != vs) {
        vs = (lower+upper)/2;
        if (rai >= accum_n_notR[vs]) lower = vs+1;
        else upper = vs;
    }
    nS[vs]--;
    nI[vs] = 1;
    nIsum = 1;
    nRsum=0;
    accum_nI[0]=nI[0]; // cumulative distribution of the number of I individuals in the subpopulation
    accum_SIpair[0] = nS[0]*nI[0]; // cumulative distribution of the number of S-I pairs in the subpopulation
    for (i=1 ; i<nV ; i++) {
        accum_nI[i]=accum_nI[i-1]+nI[i];
        accum_SIpair[i]=accum_SIpair[i-1] + nS[i]*nI[i];
    }

    // initialization done

    // dynamics
    t=0.0;
    if (tr >= 1)
        cout << "nan nan" << endl; // separating the results for different runs. This is useful for plotting the results by Python matplotlib

    while (nIsum>0) { // I exists 

        total_rate = beta * accum_SIpair[nV-1] + mu * nIsum + D * accum_n_notR[nV-1];
        // accum_nP[nV-1] = # individuals excluding R individuals. We do not need to simulate mobility of R individuals.
        t += - log((genrand_int32()+0.5)/4294967296.0) / total_rate; 
        ra = (genrand_int32()+0.5)/4294967296.0 * total_rate;

        if (ra < D*accum_n_notR[nV-1]) { // one individual moves to a neighboring subpopulation
            rai = (int)(ra/D); // now 0 <= rai <= accum_nP[nV-1]-1
            lower = 0; // determine the subpopulation containing the individual that is moving
            upper = nV-1;
            vs = upper;
            while (lower != vs || upper != vs) {
                vs = (lower+upper)/2;
                if (rai >= accum_n_notR[vs]) lower = vs+1;
                else upper = vs;
            } // an individual in subpopulation vs moves to a neighboring subpopulation
            ve = nei_list[accum_k[vs] - genrand_int32()%k[vs] - 1]; // accum_k[vs-1] <= ve < accum_k[vs] with the convention accum_k[-1]=0. An individual in subpopulation vs moves to subpopulation ve.

            if (genrand_int32()%(nS[vs]+nI[vs]) < nI[vs]) { // one infected individual moves from subpopulation vs to subpopulation ve

                for (i=ve ; i<nV ; i++) { // ** old code was wrong. It had i<vs, but i<nV should be correct.
                    accum_nI[i]++;
                    accum_SIpair[i] += nS[ve];
                    // incdrease in the number of S-I pairs owing to the entry of an infected individual to subpopulation ve = nS[ve]*(nI[ve]+1) - nS[ve]*nI[ve] = nS[ve]
                }

                for (i=vs ; i<nV ; i++) {
                    accum_nI[i]--;
                    accum_SIpair[i] -= nS[vs];
                    // increase in the number of S-I pairs owing to the leave of an infected individual from subpopulation vs = nS[vs]*(nI[vs]-1) - nS[vs]*nI[vs] = -nS[vs]
                }    

                nI[vs]--;
                nI[ve]++;

            } else { // one susceptible individual moves from subpopulation vs to subpopulation ve

                for (i=vs ; i<nV ; i++) // ** old code was wrong. It had i<ve, but i<nV should be correct.
                    accum_SIpair[i] -= nI[vs];
                    // increase in the number of S-I pairs owing to the leave of a susceptible individual from subpopulation vs = (nS[vs]-1)*nI[vs] - nS[vs]*nI[vs] = -nI[vs]

                for (i=ve ; i<nV ; i++)
                    accum_SIpair[i] += nI[ve];
                    // increase in the number of S-I pairs owing to the entry of a susceptible individual to subpopulation ve = (nS[ve]+1)*nI[ve] - nS[ve]*nI[ve] = nI[ve]
               
                nS[vs]--;
                nS[ve]++;

            } 

            // one individual (S or I) anyways moves from vs to ve
            if (vs<ve) {
                for (i=vs ; i<ve ; i++) accum_n_notR[i]--;
            } else {
                for (i=ve ; i<vs ; i++) accum_n_notR[i]++;
            }

        } else { // no individual moves. So, either one individual either gets infected or recovers.
            ra -= D*accum_n_notR[nV-1]; // now 0 < ra < beta*accum_SIpair[nV-1] + nIsum

        	if (ra < mu * nIsum) { // one individual recovers
                // determine which individual recovers using the bisection method
                rai=(int)ra; // now 0 <= rai <= nIsum-1
                lower=0;
                upper=nV-1;
                vs=upper;
                while (lower != vs || upper != vs) {
                    vs = (lower+upper)/2;
                    if (rai >= accum_nI[vs]) lower = vs+1;
                    else upper=vs;
                } // an infected individual in subpopulation vs recovers
                inc_SIpair = -nS[vs]; // increase in the number of S-I pairs = nS[vs]*(nI[vs]-1) - nS[vs]*nI[vs];
                nI[vs]--;
                nIsum--;
                nRsum++;
                for (i=vs ; i<nV ; i++) {
                    accum_n_notR[i]--;
            	    accum_nI[i]--;
            	    accum_SIpair[i] += inc_SIpair;
                }

            } else { // one susceptible individual gets infected
                // determine which individual gets infected using the bisection method
                rai = (int)((ra - mu*nIsum)/beta); // now 0 <= rai <= accum_SIpair[nV-1]-1
                // determine a tentative starting I vertex
                lower=0;
                upper=nV-1;
                vs=upper;
                while (lower != vs || upper != vs) {
                    vs = (lower+upper)/2;
                    if (rai >= accum_SIpair[vs]) lower = vs+1;
                    else upper=vs;
                } // rai is no longer necessary
                // an individual in subpopulation vs gets infected
                inc_SIpair = nS[vs] - nI[vs] - 1; // increase in the number of S-I pairs = (nS[vs]-1)*(nI[vs]+1) - nS[vs]*nI[vs]
                nS[vs]--;
                nI[vs]++;
                nIsum++;
                for (i=vs ; i<nV ; i++) {
                    accum_nI[i]++;
                    accum_SIpair[i] += inc_SIpair;
                }
        	} // infection of one individual completed

    	} // either one recovery or infection event completed
    	
        cout << t << " " << (double)(nP-nIsum-nRsum)/nP << " " << (double)nIsum/nP << " " << (double)nRsum/nP << endl; // a single-transition event completed
    } // one trial completed
    } // all the trials completed

    free_memory(nei_list, k, accum_k); // free memory
    return 0;
}
