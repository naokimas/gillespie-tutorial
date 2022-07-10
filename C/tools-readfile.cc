int compare (const void *a, const void *b) { // to use qsort
  return (*(int*)a-*(int*)b);
}

int readfile (char *infile_name, int *nV, int *nE, int **nei_list, int **k, int **accum_k) {
/* k[i] = degree of node i
accum_k[i] = sum of degrees of nodes from 0 to i
*/

    char dummy[8];
    int i, j, vs, ve;

    ifstream fin(infile_name); // edge list
    if (!fin) {
        cerr << "readfile: Cannot open input edge file" << endl;
        exit(8);
    }
    
    fin >> dummy >> *nV >> *nE; // The first line of the file has '#', the number of nodes, and the number of edges, space separated
    cerr << *nV << " nodes, " << *nE << " edges" << endl; // undirected net assumed
    *k = (int *)malloc((*nV)*sizeof(int)); // allocate memory
    for (i = 0 ; i < *nV ; i++)
        (*k)[i] = 0; // initialization
    int *E; // list of edges
    E = (int *)malloc(4*(*nE)*sizeof(int)); // allocate memory
    for (i = 0 ; i < *nE ; i++) { // read all edges
        fin >> vs >> ve >> dummy; // Each line of the input file has the index of the two nodes connected by an edge in the 1st and 2nd columns, and edge weight, which we do not use in this code, in the 3rd column. 
    // Note that in this code the node index starts from 0, not from 1.
    // And in the input file the node index starts from 1.
        vs--; ve--; // node index corrected here. Now the node index starts from 0.
        E[4*i+3]=E[4*i]=vs; // i-th edge
        E[4*i+2]=E[4*i+1]=ve;
        (*k)[vs]++; (*k)[ve]++; // undirected network assumed
    }
    fin.close();

    *accum_k = (int *)malloc((*nV)*sizeof(int)); // allocate memory
    (*accum_k)[0] = (*k)[0];
    for (i = 1 ; i < *nV ; i++) (*accum_k)[i] = (*accum_k)[i-1] + (*k)[i];

    qsort(E, 2* *nE, 2*sizeof(int), compare);

    *nei_list = (int *)malloc(2 * *nE *sizeof(int)); // allocate memory
    for (i = 0 ; i < 2 * *nE ; i++) 
        (*nei_list)[i] = E[2*i+1];
        // Now nei_list[x] where x = accum_k[i-1], ..., accum_k[i]-1 contains the neighbors of node i, with the convention accum_k[-1] = 0.
        // Note that accum_k[nV-1] = nV-1

	free(E);
	return 0;
}

int free_memory(int *nei_list, int *k, int *accum_k) {
    free(nei_list);
    free(k);
    free(accum_k);
    return 0;
}