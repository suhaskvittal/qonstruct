"""
    author: Suhas Vittal
    date:   1 December 2023

    Color code tanner graphs.
"""

from qonstruct.code_builder.base import *
from qonstruct.utils import *

import networkx as nx

from itertools import permutations

def color_tanner_graph(gr: nx.Graph):
    xsgr = make_support_graph(gr, 'x')
    xcm = nx.coloring.greedy_color(xsgr, strategy='connected_sequential_bfs')
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

def make_semihyperbolic(gr: nx.Graph, period: list[int]) -> nx.Graph:
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
    
    def compute_cycle(supp):
        if len(supp) > 10:
            return supp
        # Adjacent qubits should have two plaquettes in common.
        comm_map = {}
        for (i, x) in enumerate(supp):
            for (j, y) in enumerate(supp):
                if i >= j:
                    continue
                comm = len(list(nx.common_neighbors(gr, x, y)))
                comm_map[(x,y)] = comm
                comm_map[(y,x)] = comm
        for cyc in permutations(supp):
            is_valid = True
            # Check if the cycle is good.
            for (i, q) in enumerate(cyc):
                _q = cyc[i-1]
                # Count boundaries as well.
                comm = comm_map[(q,_q)]
                bboth_have = [True, True, True]
                for x in gr.neighbors(q):
                    bboth_have[gr.nodes[x]['color']] = False
                for x in gr.neighbors(_q):
                    bboth_have[gr.nodes[x]['color']] = False
                comm += 2*sum(bboth_have)
                if comm < 4:
                    is_valid = False
                    break
            if is_valid:
                return cyc
        exit(1)
    # First, identify the qubits we are expanding.
    plaq_to_cycle = {}
    for (_plaq, supp) in gr.graph['plaquette_support_map'].items():
        c = gr.graph['plaquette_color_map'][_plaq]
        cyc = compute_cycle(supp)
        plaq_to_cycle[_plaq] = cyc
        for i in range(0, len(cyc), period[c]):
            q = cyc[i]
            if q in replaced:
                continue
            new_q = list(range(n, n+6))
            qa,qb,qc,qd,qe,qf = new_q[:]
            for _q in new_q:
                add_data_qubit(ngr, _q)
            new_q.append(q)
            n += 6
            replaced[q] = new_q
            # Add checks for Steane's code. q is in all checks (it is the central qubit).
            for (_c, steane_supp) in enumerate([[qa, qb, qd, q], [qb, qc, q, qe], [qd, q, qe, qf]]):
                plaquettes.append((_c, steane_supp))
    # Now make the checks
    for (_plaq, cyc) in plaq_to_cycle.items():
        c = gr.graph['plaquette_color_map'][_plaq]
        # Find first qubit that has already been replaced.
        print('Expanding plaquette', cyc)
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
        print(f'Check weight: {len(cyc)} --> {len(big_supp)}')
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

