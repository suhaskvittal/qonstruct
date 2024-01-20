"""
    author: Suhas Vittal
    date:   2 January 2024
"""

from qonstruct.codes.base import *

import networkx as nx

def read_tanner_graph_file(input_file: str) -> nx.Graph:
    """
        This reads a Tanner graph from an input file. Each line should be
        as follows:
            <O?><X/Z><number>,<data-qubit>,<data-qubit>,...
        If a line begins with an 'O', it is considered an observable declaration.
        Otherwise, it is a check declaration. Example (of not a real code):

            Z0,0,1,2
            Z1,2,3,4
            X0,1,2,3
            X2,0,2,4
            OZ0,0,1,3
            OX0,0,2,3
    """
    reader = open(input_file, 'r')
    lines = reader.readlines()
    reader.close()

    gr = nx.Graph()

    data_qubits, x_checks, z_checks, x_obs_list, z_obs_list = [], [], [], [], []
    check_map = {}  # Add the checks after parsing
    for ln in lines:
        line_data = ln.split(',')
        ch_s = line_data[0]
        support = [int(x) for x in line_data[1:]]
        if ch_s[0] == 'O':
            # This is an observable.
            if ch_s[1] == 'X' or ch_s[1] == 'x':  # X obs.
                x_obs_list.append(support)
            else:
                z_obs_list.append(support)
        else:
            check_map[ch_s] = support
        # Create data qubits.
        for x in support:
            if not nx.has_node(x):
                gr.add_node(x, node_type='data')
                data_qubits.append(x)
    # Add the data qubits now:
    n = len(data_qubits)
    for (ch_s, support) in check_map.items():
        stabilizer = 'x' if ch_s[0] == 'X' or ch_s[0] == 'x' else 'z'
        gr.add_node(n, node_type=stabilizer) 
        if stabilizer == 'x':
            x_checks.append(n)
        else:
            z_checks.append(n)
        # Add edges to support.
        for x in support:
            gr.add_edge(x, n)
        n += 1
    return gr

def write_tanner_graph_file(tanner_graph: nx.Graph, output_file: str) -> None:
    """
        This code will write a Tanner graph to an output file.
    """
    writer = open(output_file, 'w')
    # Map data qubit nodes to indices.
    x_ctr, z_ctr = 0, 0
    for s in tanner_graph.graph['checks']['all']:
        if get_check_property(tanner_graph, s, 'node_type') == 'x':
            writer.write('X%d' % x_ctr)
            x_ctr += 1
        else:
            writer.write('Z%d' % z_ctr)
            z_ctr += 1
        for d in tanner_graph.neighbors(s):
            writer.write(',%d' % tanner_graph.nodes[d]['label'])
        writer.write('\n')
    # Write logical observables of the code.
    x_obs_list, z_obs_list = tanner_graph.graph['obs_list']['x'], tanner_graph.graph['obs_list']['z']
    for (i, obs) in enumerate(x_obs_list):
        writer.write('OX%d' % i)
        for d in obs:
            writer.write(',%d' % d)
        writer.write('\n')
    for (i, obs) in enumerate(z_obs_list):
        writer.write('OZ%d' % i)
        for d in obs:
            writer.write(',%d' % d)
        writer.write('\n')
    writer.close()

