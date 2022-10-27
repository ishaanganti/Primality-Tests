# APR-CL-
A slightly scuffed implementation of the APR-CL algorithm in Python featuring timer methods to assess the algorithm's runtime. Parts of the implementation, including much of the Jacobi Sum class, are from other attempts of implementing this algorithm that I found online or translated from C to Python from FLINT's source code. This code was used to assess the runtime of APR-CL for increasing bit-sized inputs when compared to other deterministic primality tests (as well as probabilistic tests, such as the Miller-Rabin primality test). 

Sources/Inspiration:
https://www.math.leidenuniv.nl/~hwl/PUBLICATIONS/1984a/art.pdf
“Implementation of a New Primality Test” by H. Cohen and A.K. Lenstra
FLINT (Fast Library for Number Theory) source code (https://github.com/wbhart/flint2/tree/trunk/aprcl)

# Miller-Rabin-

A simple implementation of the Miller Rabin primality test featuring a timer to compare runtimes. 
