"""
    author: Suhas Vittal
    date:   21 January 2024
"""

from qonstruct.code_builder.base import *

import networkx as nx
import numpy as np

def make_support_graph(tanner_graph: nx.Graph, op_type: str) -> nx.Graph:
    gr = nx.Graph()
    for x in tanner_graph.nodes():
        if tanner_graph.nodes[x]['node_type'] == op_type:
            gr.add_node(x)
    for x in gr.nodes():
        for y in gr.nodes():
            if x == y:
                continue
            if len(list(nx.common_neighbors(tanner_graph, x, y))) > 0:
                gr.add_edge(x, y)
    return gr
