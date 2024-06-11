"""
    author: Suhas Vittal
    date:   24 April 2024

    This file contains code for laying out Tanner graphs with a bus architecture.
    Requires cplex.
"""

from qonstruct.code_builder.base import *

import networkx as nx

import docplex.mp
from docplex.mp.model import Model

def compute_layout(tanner_graph: nx.Graph, radius: int) -> dict[int, tuple[int,int]]:
    """
        Computes a grid layout for the given Tanner graph.
    """
    prog = Model(name='layout')

    q2x, q2y = {}, {}
    max_x, max_y = prog.integer_var(name='xmax'), prog.integer_var(name='ymax')
    for q in tanner_graph.nodes():
        x, y = prog.integer_var(name=f'x{q}'), prog.integer_var(name=f'y{q}')
        prog.add_constraint(x >= 0)
        prog.add_constraint(y >= 0)
        prog.add_constraint(max_x >= x)
        prog.add_constraint(max_y >= y)
        q2x[q] = x
        q2y[q] = y
    # Add uniqueness constraint between each xy pair.
    for q1 in tanner_graph.nodes():
        x1, y1 = q2x[q1], q2y[q1]
        for q2 in tanner_graph.nodes():
            if q1 > q2:
                continue
            x2, y2 = q2x[q2], q2y[q2]
            prog.add_constraint(x1 != x2)
            prog.add_constraint(y1 != y2)
    # Add distance constraints for each edge.
    opt_sum_array = []
    for (q1, q2) in tanner_graph.edges():
        x1, y1, x2, y2 = q2x[q1], q2y[q1], q2x[q2], q2y[q2]
        prog.add_constraint(x1**2 - 2*x1*x2 + x2**2 <= radius + 1000000*(x1==0) + 1000000*(x2==0))
        prog.add_constraint(y1**2 - 2*y1*y2 + y2**2 <= radius + 1000000*(y1==0) + 1000000*(y2==0))
        opt_sum_array.extend([x1**2 - 2*x1*x2 + x2**2, y1**2 - 2*y1*y2 + y2**2])
    prog.minimize(prog.sum(opt_sum_array))
    soln = prog.solve()
    
