# aks

AKS primality test implementation in Python

Basd on the algorithm in https://en.wikipedia.org/wiki/AKS_primality_test

## how it works

The AKS primality test is based on the fact that `(a + b)^n mod (X^r - 1, n) = (a^n + b^n)` due to the fact that `nCk(n, k) mod n = 0` if k is not 1 or n.
The general algorithm is outlined in the Wikipedia article mentioned above.

## implementation details
- uses karatsuba's algorithm for multiplication
