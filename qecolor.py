"""
    author: Suhas Vittal
    date:   1 December 2023

    Named qecolor to avoid conflicts with other packages named color.
"""

import networkx as nx

def make_hexagonal_tanner_graph(d: int, both_at_once=False) -> nx.Graph:
    """
        This function will return a tanner graph for a hexagonal color
        code (weight-6 RGB plaquettes). The graph will have properties:
            (1) data --> the data qubits of the code
            (2) x --> x stabilizers
            (3) z --> z stabilizers
            (4) checks --> all stabilizers
            (5) obs_list --> {'x': x observables, 'z': z observables}
            (6) plaquettes --> the plaquettes of the code
    """
    gr = nx.Graph()
    data_qubits, x_checks, z_checks, plaquettes = [], [], [], []
    # Some geometric properties of the code:
    offset = 2
    side_len = (3*d - 1)//2
    # In our initial pass, we will go through the structure, and identify
    # each qubit.
    loc_map = {}
    check_locs = []
    x_obs, z_obs = [], []

    n = 0
    for r in range(side_len):
        row_off = offset
        for c in range(r+1):
            if row_off == 0:
                check_locs.append((r, c))
            else:
                loc_map[(r, c)] = n
                gr.add_node(n, node_type='data')
                if c == 0:
                    # This is the left edge of the triangle, which is a
                    # X and Z operator.
                    x_obs.append(n)
                    z_obs.append(n)
                data_qubits.append(n)
                n += 1
            row_off = (row_off+1) % 3
        offset = (offset+1) % 3
    # Now, create all the checks and complete the Tanner graph.
    gr.graph['plaquette_check_map'] = {}
    gr.graph['plaquette_support_map'] = {}
    gr.graph['plaquette_color_map'] = {}
    plaq = n

    x_obs = data_qubits
    z_obs = x_obs
    
    get_loc = lambda _r, _c: loc_map[(_r, _c)] if (_r, _c) in loc_map else None

    for (i, j) in check_locs:
        # Compute support:
        #
        #     a   b           a  b
        #   f   P   c    -->  f  P  c
        #     e   d              e  d
        #
        a = get_loc(i-1, j-1)
        b = get_loc(i-1, j)
        c = get_loc(i, j+1)
        d = get_loc(i+1, j+1)
        e = get_loc(i+1, j)
        f = get_loc(i, j-1)

        gr.graph['plaquette_check_map'][plaq] = []
        gr.graph['plaquette_support_map'][plaq] = [b, c, e, d, a, f]
        gr.graph['plaquette_color_map'][plaq] = i % 3
        for stabilizer in ['x', 'z']:
            # Create check vertex.
            if both_at_once:
                if stabilizer == 'z':
                    schedule_order = [b, c, d, a, f, e, None]
                else:
                    schedule_order = [None, b, a, f, c, d, e]
            else:
                schedule_order = [b, c, e, d, a, f]
            gr.add_node(n,\
                        node_type=stabilizer,\
                        color=i%3, plaquette=plaq,\
                        schedule_order=schedule_order)
            for x in [a, b, c, d, e, f]:
                if x is not None:
                    gr.add_edge(x, n)
            gr.graph['plaquette_check_map'][plaq].append(n)
            if stabilizer == 'x':
                x_checks.append(n)
            else:
                z_checks.append(n)
            n += 1
        plaquettes.append(plaq)
        plaq += 1
    gr.graph['data'] = data_qubits
    gr.graph['x'] = x_checks
    gr.graph['z'] = z_checks
    gr.graph['checks'] = [*x_checks, *z_checks]
    gr.graph['plaquettes'] = plaquettes
    gr.graph['obs_list'] = {'x': [x_obs], 'z': [z_obs]}
    return gr

