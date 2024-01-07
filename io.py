"""
    author: Suhas Vittal
    date:   2 January 2024
"""

import networkx as nx

def write_tanner_graph_file(tanner_graph: nx.Graph, output_file: str) -> None:
    """
        This code will write a Tanner graph to an output file.
    """
    writer = open(output_file, 'w')
    # Map data qubit nodes to indices.
    data_to_index = { d : i for (i, d) in enumerate(tanner_graph.graph['data']) }
    # Write stabilizers to the file.
    x_ctr, z_ctr = 0, 0
    for s in tanner_graph.graph['checks']:
        if tanner_graph.nodes[s]['node_type'] == 'x':
            writer.write('X%d' % x_ctr)
            x_ctr += 1
        else:
            writer.write('Z%d' % z_ctr)
            z_ctr += 1
        for d in tanner_graph.neighbors(s):
            writer.write(',%d' % data_to_index[d])
        writer.write('\n')
    # Write logical observables of the code.
    x_obs_list, z_obs_list = tanner_graph.graph['obs_list']['x'], tanner_graph.graph['obs_list']['z']
    for (i, obs) in enumerate(x_obs_list):
        writer.write('OX%d' % i)
        for d in obs:
            writer.write(',%d' % data_to_index[d])
        writer.write('\n')
    for (i, obs) in enumerate(z_obs_list):
        writer.write('OZ%d' % i)
        for d in obs:
            writer.write(',%d' % data_to_index[d])
        writer.write('\n')
    writer.close()
