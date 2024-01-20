"""
    author: Suhas Vittal
    date:   8 November 2023
"""

from qonstruct.primitives.matrix import *
from qonstruct.primitives import ldpc

from qonstruct.codes.base import *

import networkx as nx

def make(seed_code: nx.Graph) -> nx.Graph:
    """
        This function returns a quantum tanner graph.
        Additional properties:
            (1) (graph) prod_map --> maps bits/checks (i, j) to a vertex in the graph.
            (2) (nodes) prod --> a tuple (i, j) corresponding to vertices in the seed graph.
    """
    gr = tanner_init()
    prod_map = {}

    dk, ck = 0, 0
    for i in seed_code.nodes():
        i_is_bit = seed_code.nodes[i]['node_type'] == 'bit'
        for j in seed_code.nodes():
            j_is_bit = seed_code.nodes[j]['node_type'] == 'bit'
            if i_is_bit == j_is_bit:
                # This is a data qubit.
                add_data_qubit(gr, dk, prod=(i, j)) 
                prod_map[(i, j)] = dk
                dk += 1
            else:
                add_check(gr, ck, 'z' if i_is_bit else 'x', [], prod=(i,j))
                prod_map[(i, j)] = ck
                ck += 1
    set_graph_property(gr, 'prod_map', prod_map)
    # Add edges to the graph. A data qubit (i, j) is connected to a stabilizer (x, y) if either
    #   (1) (i, x) is an edge in the seed code and j == y, or
    #   (2) (j, y) is an edge in the seed code and i == x.
    for d in gr.graph['data_qubits']:
        i, j = get_data_qubit_property(gr, d, 'prod')
        for s in gr.graph['checks']['all']:
            x, y = get_check_property(gr, s, 'prod')
            if (seed_code.has_edge(i, x) and j == y) or (seed_code.has_edge(j, y) and i == x):
                add_to_support(gr, s, d)
    compute_observables(gr, seed_code)
    return gr
        
def compute_observables(gr: nx.Graph, seed_code: nx.Graph) -> None:
    """
        This function computes the logical observables of the HGP code and sets them.
    """
    H = ldpc.tanner_graph_to_parity_check_matrix(seed_code)
    Ht = H.T
    ker_H = ker(H)
    ker_Ht = ker(Ht)
    m, n = H.shape
    # X observables:    { ker(H) x Fn/row(H) | 0 } U { 0 | Fm/row(Ht) x ker(Ht) }
    # Z observables:    { Fn/row(H) x ker(H) | 0 } U { 0 | ker(Ht) x Fm/row(Ht) }
    nquo = quotient(row(H))
    mquo = quotient(row(Ht))

    def _cartprod(lhs, rhs, vertices):
        obs_list = []
        for v in lhs:
            for w in rhs:
                obs = []
                for i in range(v.shape[0]):
                    x = vertices[i]
                    if v[i] == 0:
                        continue
                    for j in range(w.shape[0]):
                        y = vertices[j]
                        if w[j] == 0:
                            continue
                        q = gr.graph['prod_map'][(x, y)]
                        obs.append(q)
                obs_list.append(obs)
        return obs_list

    bits, checks = seed_code.graph['bits'], seed_code.graph['checks']
    x_ops = _cartprod(ker_H, nquo, bits)
    x_ops.extend(_cartprod(mquo, ker_Ht, checks))
    z_ops = _cartprod(nquo, ker_H, bits)
    z_ops.extend(_cartprod(ker_Ht, mquo, checks))

    for x in x_ops:
        add_observable(gr, x, 'x')
    for x in z_ops:
        add_observable(gr, x, 'z')
