"""
    author: Suhas Vittal
    date:   17 January 2024

    Utility functions for constructing codes.
"""

import networkx as nx

from collections import defaultdict

def tanner_init() -> nx.Graph:
    gr = nx.Graph()
    # Initialize structures for the graph.
    gr.graph['data_qubits'] = []
    gr.graph['checks'] = {'all': [], 'x': [], 'z': []}
    gr.graph['obs_list'] = {'x': [], 'z': []}
    # These are specific to color codes:
    gr.graph['plaquettes'] = []
    gr.graph['plaquette_support_map'] = {}
    gr.graph['plaquette_color_map'] = {}
    gr.graph['plaquette_check_map'] = defaultdict(list)
    # Hidden properties:
    gr.graph['_index'] = 0
    gr.graph['_data_to_index'] = {}
    gr.graph['_check_to_index'] = {}
    return gr

def add_data_qubit(gr: nx.Graph, q: int, **kwargs) -> int:
    n = gr.graph['_index']
    gr.add_node(n, node_type='data', label=q)
    for (k, v) in kwargs.items():
        gr.nodes[n][k] = v
    gr.graph['data_qubits'].append(q)

    gr.graph['_data_to_index'][q] = n
    gr.graph['_index'] += 1
    return n

def contains_data_qubit(gr: nx.Graph, q: int) -> bool:
    return q in gr.graph['_data_to_index']

def get_data_qubit_node(gr: nx.Graph, q: int) -> int:
    return gr.graph['_data_to_index'][q]

def add_check(gr: nx.Graph, check: int, check_type: str, support: list[int], **kwargs) -> int:
    n = gr.graph['_index']
    gr.add_node(n, node_type=check_type, label=check)
    for (k, v) in kwargs.items():
        gr.nodes[n][k] = v
    gr.graph['checks']['all'].append(check)
    gr.graph['checks'][check_type].append(check)

    gr.graph['_check_to_index'][check] = n
    gr.graph['_index'] += 1
    # Add edges to data qubits:
    for q in support:
        if contains_data_qubit(gr, q):
            x = get_data_qubit_node(gr, q)
        else:
            x = add_data_qubit(gr, q)
        gr.add_edge(x, n)
    return n

def add_to_support(gr: nx.Graph, check: int, q: int) -> None:
    n = gr.graph['_check_to_index'][check]
    x = gr.graph['_data_to_index'][q]
    gr.add_edge(n, x)
        
def contains_check(gr: nx.Graph, check: int) -> bool:
    return check in gr.graph['_check_to_index']

def get_check_node(gr: nx.Graph, check: int) -> int:
    return gr.graph['_check_to_index'][check]

def get_support(gr: nx.Graph, check: int) -> list[int]:
    return [gr.nodes[x]['label'] for x in gr.neighbors(get_check_node(gr, check))]

def add_observable(gr: nx.Graph, observable: list[int], obs_type: str) -> None:
    gr.graph['obs_list'][obs_type].append(observable)

def get_graph_property(gr: nx.Graph, property_key: str) -> any:
    return gr.graph['property_key']

def get_data_qubit_property(gr: nx.Graph, q: int, property_key: str) -> any:
    n = get_data_qubit_node(gr, q)
    return gr.nodes[n][property_key]

def get_check_property(gr: nx.Graph, check: int, property_key: str) -> any:
    n = get_check_node(gr, check)
    return gr.nodes[n][property_key]

def set_graph_property(gr: nx.Graph, property_key: str, property_value: any) -> None:
    gr.graph[property_key] = property_value

def set_data_qubit_property(gr: nx.Graph, q: int, property_key: str, property_value: any) -> None:
    n = get_data_qubit_node(gr, q)
    gr.nodes[n][property_key] = property_value

def set_check_property(gr: nx.Graph, check: int, property_key: str, property_value: any) -> None:
    n = get_check_node(gr, check)
    gr.nodes[n][property_key] = property_value

# Specific Code Functions:

def add_plaquette(gr: nx.Graph, plaquette: int, support: list[int], color: int) -> None:
    gr.graph['plaquettes'].append(plaquette)
    gr.graph['plaquette_support_map'][plaquette] = support
    gr.graph['plaquette_color_map'][plaquette] = color

def set_plaquette(gr: nx.Graph, check: int, plaquette: int) -> None:
    if plaquette not in gr.graph['plaquettes']:
        add_plaquette(gr, plaquette)
    gr.graph['plaquette_check_map'][plaquette].append(check)

def get_plaquette(gr: nx.Graph, check: int) -> int:
    n = get_check_node(gr, check)
    return gr.nodes[n]['plaquette']


