"""
    author: Suhas Vittal
    date:   8 November 2023
"""

from qonstruct.primitives.matrix import *

import networkx as nx

def tanner_graph_to_parity_check_matrix(gr: nx.Graph) -> np.ndarray:
    """
        Converts a Tanner graph to a parity check matrix. Each
        row encodes a check, and each column corresponds to a bit.
    """
    bits, checks = gr.graph['bits'], gr.graph['checks']
    m, n = len(checks), len(bits)

    bit_to_index = { b:j for (j,b) in enumerate(bits) }

    H = np.zeros((m, n))
    for (i, c) in enumerate(checks):
        for b in gr.neighbors(c):
            j = bit_to_index[b]
            H[i,j] = 1
    return H

def make_regular_tanner_graph(r: int, c: int, s: int) -> nx.Graph:
    """
        Makes a Tanner graph with c*s bits and r*s checks such that
            (1) Each bit is connected to r checks.
            (2) Each check is connected to c bits.
        and this attempts to achieve a high girth.

        The Tanner graph will be bipartite and contain keys:
            "bits" --> all bit vertices
            "checks" --> all check vertices
        Each vertex will also have the data:
            "node_type" --> either 'bit' or 'check'
    """
    gr = nx.Graph()
    bits = []
    checks = []
    # Create nodes.
    for i in range(c*s):
        gr.add_node(i, node_type='bit')
        bits.append(i)
    for i in range(c*s, (c+r)*s):
        gr.add_node(i, node_type='check')
        checks.append(i)
    gr.graph['bits'] = bits
    gr.graph['checks'] = checks
    # We will add edges as per the Progressive Edge Growth algorithm,
    # which yields graphs with good girth.
    for i in range(c*s):
        progressive_edge_growth(gr, i, r)
    return gr

def progressive_edge_growth(gr: nx.Graph, i: int, dg: int) -> None:
    checks = gr.graph['checks']

    def _terrace(curr: set[int]) -> tuple[set[int], set[int]]:
        # Computes next terrace.
        next = set()
        for x in curr:
            for y in gr.neighbors(x):
                # Note that y is a data qubit, so we want the neighbors of y instead.
                for z in gr.neighbors(y):
                    next.add(z)
        compl_next = set(c for c in checks if c not in next)
        return next, compl_next
    
    def _search() -> tuple[set[int], int]:
        # Returns potential edge candidates, and depth of search.
        curr = set(gr.neighbors(i))
        compl_curr = set(c for c in checks if c not in curr)
        d = 0
        while len(curr) < len(checks):
            next, compl_next = _terrace(curr)
            if len(compl_next) == 0 or len(next) == len(curr):
                break
            curr = next
            compl_curr = compl_next
            d += 1
        return compl_curr, d

    def _lookahead(j: int) -> int:
        # Figure out how deep we would search for a terrace after adding an edge between i and j.
        # The deeper, the better.
        gr.add_edge(i, j)
        _, d = _search()
        gr.remove_edge(i, j)
        return d

    for k in range(dg):
        if k == 0:
            candidates = checks
        else:
            candidates, _ = _search() 
        min_degree = min(gr.degree(c) for c in candidates)
        candidates = [c for c in candidates if gr.degree(c) == min_degree]
        j = max(candidates, key=lambda x: _lookahead(x))
        gr.add_edge(i, j)
