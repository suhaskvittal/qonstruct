"""
    author: Suhas Vittal
    date:   1 December 2023

    Color code tanner graphs.
"""

from qonstruct.code_builder.base import *
from qonstruct.utils import make_support_graph

import networkx as nx

def color_tanner_graph(gr: nx.Graph):
    xsgr = make_support_graph(gr, 'x')
    xcm = nx.coloring.greedy_color(xsgr, strategy='DSATUR')
    max_color = max(x for (_, x) in xcm.items())
    #print(f'computed {max_color+1}-coloring')
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
