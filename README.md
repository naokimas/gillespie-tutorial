# gillespie-tutorial

- mt19937ar.c: Mersenne Twister. Because this is the developer's code, we are probably not allowed to upload it. Need to sort it out.

### Voter model on networks ###
- voter-net.cc: main code, calling mt19937ar.c

Usage:
```
g++ voter-net.cc -o voter-net.out
./voter-net.out infilename > result-voter-net
```

- rrg-N100k5.mat: A regular random graph with N=100 nodes and degree 5.

- plot-voter.py: Plotting

