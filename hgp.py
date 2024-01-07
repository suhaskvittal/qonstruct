"""
    author: Suhas Vittal
    date:   8 November 2023
"""

from qonstruct.primitives.matrix import *
from qonstruct import ldpc

import networkx as nx

def make_hgp_quantum_tanner_graph(seed_code: nx.Graph) -> nx.Graph:
    """
        This function returns a quantum tanner graph. The graph will
        have properties:
            (1) data --> the data qubits of the code
            (2) x --> x stabilizers
            (3) z --> z stabilizers
            (4) checks --> all stabilizers
            (5) obs_list --> {'x': x observables, 'z': z observables}
            (6) prod_map --> maps bits/checks (i, j) to a vertex in the graph.
        Each vertex will have properties:
            (1) node_type --> one of 'data', 'x', or 'z'
            (2) prod --> a tuple (i, j) corresponding to vertices in the seed graph.
    """
    gr = nx.Graph()
    data, xstab, zstab = [], [], []
    prod_map = {}

    k = 0
    for i in seed_code.nodes():
        i_is_bit = seed_code.nodes[i]['node_type'] == 'bit'
        for j in seed_code.nodes():
            j_is_bit = seed_code.nodes[j]['node_type'] == 'bit'
            if i_is_bit == j_is_bit:
                # This is a data qubit.
                gr.add_node(k, node_type='data', prod=(i, j))
                data.append(k)
            elif i_is_bit:
                # Z stabilizer
                gr.add_node(k, node_type='z', prod=(i, j))
                zstab.append(k)
            else:
                # X stabilizer
                gr.add_node(k, node_type='x', prod=(i, j))
                xstab.append(k)
            prod_map[(i, j)] = k
            k += 1
    # Add edges to the graph. A data qubit (i, j) is connected to a stabilizer (x, y) if either
    #   (1) (i, x) is an edge in the seed code and j == y, or
    #   (2) (j, y) is an edge in the seed code and i == x.
    all_stab = [*xstab, *zstab]
    for d in data:
        i, j = gr.nodes[d]['prod']
        for s in all_stab:
            x, y = gr.nodes[s]['prod']
            if (seed_code.has_edge(i, x) and j == y)\
            or (seed_code.has_edge(j, y) and i == x):
                gr.add_edge(d, s)
    gr.graph['data'] = data
    gr.graph['x'] = xstab
    gr.graph['z'] = zstab
    gr.graph['checks'] = all_stab
    gr.graph['prod_map'] = prod_map

    x_obs_list, z_obs_list = get_observables(gr, seed_code)
    gr.graph['obs_list'] = {'x': x_obs_list, 'z': z_obs_list}
    return gr
        
def get_observables(hgp_code: nx.Graph, seed_code: nx.Graph) -> tuple[list[list[int]], list[list[int]]]:
    """
        This function computes the logical observables of the HGP code.
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
                        q = hgp_code.graph['prod_map'][(x, y)]
                        obs.append(q)
                obs_list.append(obs)
        return obs_list

    bits, checks = seed_code.graph['bits'], seed_code.graph['checks']
    x_obs_list = _cartprod(ker_H, nquo, bits)
    x_obs_list.extend(_cartprod(mquo, ker_Ht, checks))
    z_obs_list = _cartprod(nquo, ker_H, bits)
    z_obs_list.extend(_cartprod(ker_Ht, mquo, checks))
    # We need to select as many observables as there are logical qubits.
    # We also need to make sure the corresponding X and Z parts share qubits.
    nlq = len(hgp_code.graph['data']) - len(hgp_code.graph['checks'])
    # Compute a basis of the logical operators.
    x_obs_list.sort(reverse=True, key=lambda x: len(x))
    z_obs_list.sort(reverse=True, key=lambda x: len(x))

    x_obs_list = x_obs_list[:nlq]
    z_obs_list = z_obs_list[:nlq]
    return x_obs_list, z_obs_list

if __name__ == '__main__':
    import parsing.cmd
    from sys import argv
    
    arg_list = parsing.cmd.parse(argv)
    output_file = parsing.cmd.try_get_string(arg_list, 'out')
    r = parsing.cmd.try_get_int(arg_list, 'r')
    c = parsing.cmd.try_get_int(arg_list, 'c')
    s = parsing.cmd.try_get_int(arg_list, 's')

    seed_code = ldpc.make_regular_tanner_graph(r, c, s)
    gr = make_hgp_quantum_tanner_graph(seed_code)
    write_tanner_graph_file(gr, output_file)

