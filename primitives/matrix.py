"""
    author: Suhas Vittal
    date:   8 November 2023

    All functions in this file implement matrix operations
    for matrices over F2.
"""

import numpy as np


def rref(M: np.ndarray) -> tuple[np.ndarray, list[int]]:
    A = M.copy()

    rr = min(A.shape[0], A.shape[1])
    pivots = []

    p = 0
    for k in range(rr):
        if A[p,k] == 0:
            # Need to get A[k, k] = 1 by swapping rows.
            for i in range(p+1, A.shape[0]):
                if A[i,k] == 1:
                    A[[p,i], :] = A[[i,p], :]
                    break
            else:
                continue
        pivots.append(k)
        # Need to set all entries below A[p,k] to 0.
        for i in range(p+1, A.shape[0]):
            if A[i,k] == 1:
                A[i,:] = np.logical_xor(A[i,:], A[k,:])
        p += 1
    # Now make A become reduced REF (RREF).
    for _k in range(0, p):
        k = p - _k - 1
        for i in range(0, k):
            if A[i,k] == 1:
                A[i,:] = np.logical_xor(A[i,:], A[k,:])
    return A, pivots

def row(M: np.ndarray) -> list[np.ndarray]:
    _, pivots = rref(M)
    return [ M[p] for p in pivots ]

def ker(M: np.ndarray) -> list[np.ndarray]:
    A, pivots = rref(M)
    basis = []
    for j in range(A.shape[1]):
        if j in pivots:
            continue
        v = np.zeros(A.shape[1])
        v[j] = 1
        for i in range(len(pivots)):
            v[i] = A[i,j]
        basis.append(v)
    return basis

def quotient(basis: list[np.ndarray]) -> list[np.ndarray]:
    n = basis[0].shape[0]
    X = np.zeros((n, len(basis)))
    for (j, v) in enumerate(basis):
        X[:,j] = v[:]
    Q = np.concatenate((X, np.eye(n)), axis=1)
    A, pivots = rref(Q)
    quo = []
    for j in pivots:
        if j < len(basis):
            continue
        quo.append(Q[:,j])
    return quo 

