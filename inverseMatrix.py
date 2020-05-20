#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the function Matrix.inv_mod()
This function is to inverse a matrix over a finite field

@para
n: matrix dimension
mod: field size

G is the matrix
H is the inverse of G over a finite field

Output
1) mse: the mean square error of GTH (If H is accurate, GTH = I)
2) Computation time to inverse G
    
@author: litang
"""

from sympy import Matrix
import numpy as np
import time

def powMod(base, power, mod):
    rst = 1
    for i in range(power):
        rst = (base * rst) % mod
    return rst

n = 15
mod = 2147483647

G = np.zeros((n, n)).astype(np.int)
for i in range(n):
    for j in range(n):
        G[i, j] = powMod(j+1, i, mod)
        
invS = time.time()
H = np.array(Matrix(G).inv_mod(mod)).astype(np.int)
invE = time.time()

check = np.matmul(G, H) % mod
ident = np.identity(n)
mse = np.sqrt(np.sum(abs((check - ident)**2))/np.sum(abs(ident**2)))

print("mse is %s"%mse)
print("time to compute is %s"%(invE-invS))# inverseoverfield
