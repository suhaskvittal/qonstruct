"""
    author: Suhas Vittal
    date:   21 January 2024
"""

from qonstruct.code_builder.base import *

import networkx as nx
import numpy as np

def binrref(m: np.ndarray) -> tuple[np.ndarray, list[int]]:
    # Returns RREF and pivot columns.
    u = m.copy()
    pivots = []
    r,p=0,0
    while r < m.shape[0] and p < m.shape[1]:
        if not u[r,p]:
            i = r+1
            while i < m.shape[0] and not u[i,p]:
                i += 1
            if i == m.shape[0]:
                p += 1
                continue
            u[[r,i],:] = u[[i,r],:]
        pivots.append(p)
        for i in range(r+1, m.shape[0]):
            if u[i,p]:
                u[i] = np.logical_xor(u[i], u[r])
        r += 1
        p += 1
    # u is REF(m). Now we must make it a RREF.
    for _i in range(len(pivots)):
        i = len(pivots)-_i-1
        j = pivots[i]
        for ii in range(i):
            if u[ii,j]:
                u[ii] = np.logical_xor(u[ii], u[i])
    return u, pivots

def row(m: np.ndarray) -> np.ndarray:
    _, pivots = binrref(m.T)
    return m[pivots]

def null(m: np.ndarray) -> np.ndarray:
    A, pivots = binrref(m)
    out = np.zeros((m.shape[1]-len(pivots), m.shape[1]), dtype=bool)
    k = 0
    for j in range(m.shape[1]):
        if j in pivots:
            continue
        out[k,j] = True
        k += 1
    for (i,jp) in enumerate(pivots):
        k = 0
        for j in range(m.shape[1]):
            if j in pivots:
                continue
            out[k,jp] = A[i,j]
            k += 1
    return out

def make_support_graph(tanner_graph: nx.Graph, op_type: str) -> nx.Graph:
    gr = nx.Graph()
    for x in tanner_graph.nodes():
        if tanner_graph.nodes[x]['node_type'] == op_type:
            gr.add_node(x)
    for x in gr.nodes():
        for y in gr.nodes():
            if x == y:
                continue
            if len(list(nx.common_neighbors(tanner_graph, x, y))) > 0:
                gr.add_edge(x, y)
    return gr
