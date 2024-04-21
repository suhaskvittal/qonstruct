"""
    author: Suhas Vittal
    date:   1 December 2023

    Color code tanner graphs.
"""

from qonstruct.code_builder.base import *
from qonstruct.utils import *

import networkx as nx

from itertools import permutations
from collections import deque

def color_tanner_graph(gr: nx.Graph):
    xsgr = make_support_graph(gr, 'x')
    xcm = nx.coloring.greedy_color(xsgr, strategy='DSATUR')
    max_color = max(x for (_, x) in xcm.items())
    print(f'computed {max_color+1}-coloring')
    # Identify plaquettes and colors.
    for (plaq, (xch, c)) in enumerate(xcm.items()):
        # Find Z checks with same support as xch.
        xsupp = [x for x in gr.neighbors(xch)]
        for zch in gr.neighbors(xsupp[0]):
            if gr.nodes[zch]['node_type'] != 'z':
                continue
            zsupp = [x for x in gr.neighbors(zch)]
            if xsupp == zsupp:
                add_plaquette(gr, plaq, xsupp, c)
                for ch in [xch, zch]:
                    gr.nodes[ch]['color'] = c
                    gr.nodes[ch]['plaquette'] = plaq
                    set_plaquette(gr, ch, plaq)
                break

def make_hexagonal(d: int, both_at_once=True) -> nx.Graph:
    """
        This function will return a tanner graph for a hexagonal color
        code (weight-6 RGB plaquettes).
    """
    gr = tanner_init()

    # Some geometric properties of the code:
    offset, side_len = 2, (3*d - 1)//2
    # In our initial pass, we will go through the structure, and identify
    # each qubit.
    loc_map, check_locs, obs = {}, [], []
    n = 0
    for r in range(side_len):
        row_off = offset
        for c in range(r+1):
            if row_off == 0:
                check_locs.append((r, c))
            else:
                loc_map[(r, c)] = n
                add_data_qubit(gr, n)
                if c == 0:
                    # This is the left edge of the triangle, which is a
                    # X and Z operator.
                    obs.append(n)
                n += 1
            row_off = (row_off+1) % 3
        offset = (offset+1) % 3
    add_observable(gr, obs, 'x')
    add_observable(gr, obs, 'z')
    # Now, create all the checks and complete the Tanner graph.
    get_loc = lambda _r, _c: loc_map[(_r, _c)] if (_r, _c) in loc_map else None

    plaq = 0
    for (i, j) in check_locs:
        a = get_loc(i-1, j-1)
        b = get_loc(i-1, j)
        c = get_loc(i, j+1)
        d = get_loc(i+1, j+1)
        e = get_loc(i+1, j)
        f = get_loc(i, j-1)

        add_plaquette(gr, plaq, [a, b, c, d, e, f], i%3)
        for stabilizer in ['x', 'z']:
            # Create check vertex.
            if both_at_once:
                if stabilizer == 'z':
#                   schedule_order = [b, c, d, a, f, e, None]
                    schedule_order = [b, c, d, e, f, a, None, None, None, None, None, None]
                else:
#                   schedule_order = [None, b, a, f, c, d, e]
                    schedule_order = [None, None, None, None, None, None, b, c, d, e, f, a]
            else:
                schedule_order = [b, c, e, d, a, f]
            add_check(gr, n, stabilizer, [a, b, c, d, e, f],
                    color=i%3, plaquette=plaq, schedule_order=schedule_order)
            set_plaquette(gr, n, plaq)
            n += 1
        plaq += 1
    return gr

def make_semihyperbolic(gr: nx.Graph, period: int) -> nx.Graph:
    ngr = tanner_init()
    # Input gr should be a hyperbolic color code.
    n = len(gr.graph['data_qubits'])
    for x in gr.graph['data_qubits']:
        add_data_qubit(ngr, x)
    # Track replacements in the semihyperbolic code.
    replaced = {}
    not_replaced = set()
    # Go through each check and try to replace a qubit with Steane's code.
    plaquettes = []
    
    def get_qubits_sharing_an_edge(q):
        adj = []
        visited = set()
        for tmp in gr.neighbors(q):
            for _q in gr.neighbors(tmp):
                if _q in visited or q == _q:
                    continue
                comm = len(list(nx.common_neighbors(gr, q, _q)))
                bboth_have = [True, True, True]
                for x in gr.neighbors(q):
                    bboth_have[gr.nodes[x]['color']] = False
                for x in gr.neighbors(_q):
                    bboth_have[gr.nodes[x]['color']] = False
                comm += 2*sum(bboth_have)
                if comm == 4:
                    adj.append(_q)
                visited.add(_q)
        return adj
    
    # First, identify the qubits we are expanding.
    ctr_map = {0: 0}
    visited = set()
    bfs = deque([0])
    while len(bfs):
        q = bfs.popleft()
        if q in visited:
            continue
        # Check if any neighbors of q have been replaced.
        for _q in get_qubits_sharing_an_edge(q):
            bfs.append(_q)
            if _q not in ctr_map:
                if ctr_map[q] == 0:
                    ctr_map[_q] = period-1
                else:
                    ctr_map[_q] = ctr_map[q]-1
            else:
                if ctr_map[q] == 0:
                    ctr_map[_q] = min(ctr_map[_q], period-1)
                else:
                    ctr_map[_q] = min(ctr_map[_q], ctr_map[q]-1)
                print(f'{q} --> {_q}: ctr = {ctr_map[q]} --> {ctr_map[_q]}')
        visited.add(q)
        # Check if we should replace q with a Steane's code.
        if ctr_map[q] > 0:
            continue
        print('Replacing qubit', q)
        new_q = list(range(n, n+6))
        qa,qb,qc,qd,qe,qf = new_q[:]
        for _q in new_q:
            add_data_qubit(ngr, _q)
        new_q.append(q)
        n += 6
        replaced[q] = new_q
        # Add checks for Steane's code. q is in all checks (it is the central qubit).
        for (c, steane_supp) in enumerate([[qa, qb, qd, q], [qb, qc, q, qe], [qd, q, qe, qf]]):
            plaquettes.append((c, steane_supp))
    # Now make the checks
    for (_plaq, cyc) in gr.graph['plaquette_support_map'].items():
        c = gr.graph['plaquette_color_map'][_plaq]
        # Find first qubit that has already been replaced.
        big_supp = [] # Overall support after adding Steane's code.
        for q in cyc:
            if q in replaced:
                if c == 0:
                    qa,qb,qc = replaced[q][2], replaced[q][4], replaced[q][5]
                elif c == 1:
                    qa,qb,qc = replaced[q][0], replaced[q][3], replaced[q][5]
                else:
                    qa,qb,qc = replaced[q][0], replaced[q][1], replaced[q][2]
                big_supp.extend([qa,qb,qc])
            else:
                big_supp.append(q)
        plaquettes.append((c, big_supp))
    # Add plaquettes and checks to ngr.
    for (plaq, (c, supp)) in enumerate(plaquettes):
        add_plaquette(ngr, plaq, supp, c)
        for typ in ['x', 'z']:
            add_check(ngr, n, typ, supp, color=c, plaquette=plaq)
            set_plaquette(ngr, n, plaq)
            n += 1
    # Update observables.
    for typ in ['x', 'z']:
        for obs in gr.graph['obs_list'][typ]:
            new_obs = []
            for q in obs:
                if q in replaced:
                    new_obs.extend(replaced[q])
                else:
                    new_obs.append(q)
            add_observable(ngr, new_obs, typ)
    return ngr

