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

def make_semihyperbolic(gr: nx.Graph, seed: str) -> nx.Graph:
    """
        Constructs a semihyperbolic color code using a hyperbolic color code
        and a seed planar code.
    """
    ngr = tanner_init()
    # Input gr should be a hyperbolic color code.
    n = len(gr.graph['data_qubits'])
    for x in gr.graph['data_qubits']:
        add_data_qubit(ngr, x)
    # Track replacements in the semihyperbolic code.
    replaced = {}
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
                    # Identify edge color between q and _q.
                    colors = { gr.nodes[x]['color'] for x in nx.common_neighbors(gr, q, _q) }
                    for (i,b) in enumerate(bboth_have):
                        if b:
                            colors.add(i)
                    if 0 not in colors:
                        c = 0
                    elif 1 not in colors:
                        c = 1
                    else:
                        c = 2
                    adj.append((c, _q))
                visited.add(_q)
        return adj
    
    # First, identify the qubits we are expanding.
    code_grp_map = {}

    visited = set()
    code_grp = 0
    for q in gr.graph['data_qubits']:
        if q in visited:
            continue
        if seed == '7_1_3':
            n = rf_7_1_3(ngr, n, q, replaced)
        elif seed == '12_2_3':
            # Find pair for q.
            qadj = get_qubits_sharing_an_edge(q)
            for (c,_q) in qadj:
                if _q in visited:
                    continue
                n = rf_12_2_3(ngr, n, q, _q, c, replaced)
                visited.add(_q)
                code_grp_map[_q] = code_grp
                break
        elif seed == '4_2_2':
            n = rf_4_2_2(ngr, n, q, replaced)
        print(q, '-->', replaced[q]['bulk'])

        visited.add(q)
        plaquettes.extend(replaced[q]['plaquettes'])
        code_grp_map[q] = code_grp
        code_grp += 1

    # Now make the checks
    for (_plaq, cyc) in gr.graph['plaquette_support_map'].items():
        c = gr.graph['plaquette_color_map'][_plaq]
        # Find first qubit that has already been replaced.
        big_supp = [] # Overall support after adding Steane's code.
        visited = set()
        for q in cyc:
            if q in code_grp_map and code_grp_map[q] in visited:
                continue
            if q in replaced:
                big_supp.extend(replaced[q]['border'][c])
                visited.add(code_grp_map[q])
            else:
                big_supp.append(q)
        plaquettes.append((c, big_supp))
        print(_plaq, cyc, big_supp, len(cyc), len(big_supp))
    # Add plaquettes and checks to ngr.
    for (plaq, (c, supp)) in enumerate(plaquettes):
        add_plaquette(ngr, plaq, supp, c)
        for typ in ['x', 'z']:
            add_check(ngr, n, typ, supp, color=c, plaquette=plaq)
            set_plaquette(ngr, n, plaq)
            n += 1
        for (_, _supp) in plaquettes:
            comm = [x for x in supp if x in _supp]
            if len(comm) % 2 == 1:
                print('Anticommutation:', supp, _supp, comm)
    # Update observables.
    for typ in ['x', 'z']:
        for obs in gr.graph['obs_list'][typ]:
            new_obs = []
            for q in obs:
                if q in replaced:
                    new_obs.extend(replaced[q]['bulk'])
                else:
                    new_obs.append(q)
            add_observable(ngr, new_obs, typ)
    return ngr

# Replacement function for Steane's code -- [[7, 1, 3]]
# Returns new "n"
def rf_7_1_3(ngr: nx.Graph, n: int, q: int, replaced: dict[int, dict]) -> int:
    qlist = list(range(n, n+6))
    a, b, c, d, e, f = qlist[:]
    for _q in qlist:
        add_data_qubit(ngr, _q)
    qlist.append(q)
    
    plaquettes = [(0, [a, b, d, q]), (1, [b, c, q, e]), (2, [d, q, e, f])]
    replaced[q] = {
        'bulk': qlist,
        'plaquettes': plaquettes,
        'border': {
            0: [c, e, f], 1: [a, d, f], 2: [a, b, c]
        }
    }
    return n+6

def rf_12_2_3(ngr: nx.Graph, n: int, q1: int, q2: int, merge_along_color: int, replaced: dict[int, dict]) -> int:
    qlist = list(range(n, n+10))
    a,b,c,d,e,\
            f,g,h,i,j = qlist[:]
    for _q in qlist:
        add_data_qubit(ngr, _q)
    qlist.extend([q1, q2])

    plaquettes = [
        (1, [e,q1,f,i]),
        (1, [g,q2,h,j]),
        (2, [a,b,e,q1]),
        (2, [c,d,q2,h]),
        (0, [b,c,q1,f,g,q2])
    ]
    # Update colors of plaquettes for merge.
    plaquettes = [ ( (c+merge_along_color % 3), p ) for (c,p) in plaquettes ]

    b0_1 = [a,e,i]
    b0_2 = [d,h,j]
    b1 = [a,b,c,d]
    b2 = [f,g,i,j]
    bord_q1 = {0: b0_1, 1: b1, 2: b2}
    bord_q2 = {0: b0_2, 1: b1, 2: b2}
    # Update border colors
    bord_q1 = { ((c+merge_along_color) % 3):b for (c,b) in bord_q1.items() }
    bord_q2 = { ((c+merge_along_color) % 3):b for (c,b) in bord_q2.items() }
    replaced[q1] = {
        'bulk': [a,b,e,f,i,q1],
        'plaquettes': plaquettes,
        'border': bord_q1
    }
    replaced[q2] = {
        'bulk': [c,d,g,h,j,q2],
        'plaquettes': plaquettes,
        'border': bord_q2
    }
    return n+10

def rf_4_2_2(ngr: nx.Graph, n: int, q: int, replaced: dict[int, dict]) -> int:
    a, b, c = n, n+1, n+2
    qlist = [a,b,c]
    for _q in qlist:
        add_data_qubit(ngr, _q)
    qlist.append(q)

    replaced[q] = {
        'bulk': qlist,
        'plaquettes': [(3, [a,b,c,q])],
        'border': {0: [a,b], 1: [a,c], 2: [b,c]}
    }
    return n+3
