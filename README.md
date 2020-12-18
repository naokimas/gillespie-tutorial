# gillespie-tutorial

- mt19937ar.c: Mersenne Twister. Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved. See the comments in the beginning of the code for the full copyright statement.                         
- tools-readfile.cc: Functions to read network files
- rrg-N100k5.mat: A regular random graph with N=100 nodes and degree 5.

- plot-trajectories-gillespie-tutorial.py: Plotting trajectories of the stochastic dynamics, producing the figures in the article

### SIR dynmics ###
#### well-mixed populations ####

TO COME

#### networks ####

TO COME

### SIR dynamics on the metapopulation network model ###
- sir-metapop.cc  
Usage:
```
g++ sir-metapop.cc
./a.out rrg-N100k5.mat > result-sir-metapop
```

### Voter model ###
#### well-mixed populations ####
- voter-wellmixed.cc  
Usage:
```
g++ voter-wellmixed.cc
./a.out 100 > result-voter-wellmixed
```
The first argument (=100 in the example above) is the number of nodes.

#### networks ####
- voter-net.cc
Usage:
```
g++ voter-net.cc
./a.out infilename > result-voter-net
```

- voter-net-binary-tree.cc: same as voter-net.cc but using the binary search tree for faster speed. Usage is the same as voter-net.cc

### Lotka-Volterra model (well-mixed populations only) ###
- lotka-volterra-wellmixed.cc
Usage:
```
g++ lotka-volterra-wellmixed.cc
./a.out > result-lv-wellmixed
```
