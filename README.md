# gillespie-tutorial

- mt19937ar.c: Mersenne Twister. Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved. See the comments in the beginning of the code for the full copyright statement.                         
- tools-readfile.cc: Functions to read network files and manage memory.
- rrg-N100k5.mat: A regular random graph with N=100 nodes and each node's degree 5.
- rrg-N1000k5.mat: A regular random graph with N=1000 nodes and each node's degree 5.

- plot-trajectories-gillespie-tutorial.py: Plotting trajectories of the stochastic dynamics, producing the figures in the article.

### SIR model ###
#### well-mixed populations ####

- sir-wellmixed.cc
```
g++ voter-wellmixed.cc
./a.out 100 > result-sir-wellmixed
```
The first argument (=100 in the example above) is the number of individuals.

#### networks ####

- sir-net.cc  
```
g++ sir-net.cc
./a.out infilename > result-sir-net
```

- sir-net-binary-tree.cc: same as sir-net.cc but using the binary search tree for faster speed.  
Usage is the same as sir-net.cc


### SIR model on the metapopulation network model ###
- sir-metapop.cc  
```
g++ sir-metapop.cc
./a.out rrg-N100k5.mat > result-sir-metapop
```

### Voter model ###
#### well-mixed populations ####
- voter-wellmixed.cc  
```
g++ voter-wellmixed.cc
./a.out 100 > result-voter-wellmixed
```
The first argument (=100 in the example above) is the number of individuals.

#### networks ####
- voter-net.cc  
```
g++ voter-net.cc
./a.out infilename > result-voter-net
```

- voter-net-binary-tree.cc: same as voter-net.cc but using the binary search tree for faster speed.  
Usage is the same as voter-net.cc

### Lotka-Volterra model (well-mixed populations only) ###
- lotka-volterra-wellmixed.cc  
```
g++ lotka-volterra-wellmixed.cc
./a.out > result-lv-wellmixed
```
