"""
    author: Suhas Vittal
    date:   1 December 2023

    Pre-written generation code.
"""

from asm.manager import AsmManager

import qecolor

def asm_hexagonal_color_code(output_file: str, distance: int, rounds: int, memory='z'):
    code = qecolor.make_hexagonal_tanner_graph(distance)
    mgr = AsmManager(code, use_plaquettes_instead_of_checks=True)
    mgr.memory = memory
    # Create flag qubits.
    plaquettes = code.graph['plaquettes']
    for plq in plaquettes:
        support = code.graph['plaquette_support_map'][plq]
        for i in range(0, len(support), 2):
            q1 = support[i]
            q2 = support[i+1]
            if q1 is None or q2 is None:
                continue
            mgr.add_flags_to(q1, q2, plq)
    # Write file.
    mgr.fopen(output_file)
    mgr.write_memory_experiment(rounds)
    mgr.fclose()

if __name__ == '__main__':
    from sys import argv
    
    from parsing.cmd import *

    inputs = parse(argv)

    output_file = try_get_string(inputs, 'out')
    d = try_get_int(inputs, 'd')
    r = try_get_int(inputs, 'r')

    memory = 'x' if 'x' in inputs else 'z'

    asm_hexagonal_color_code(output_file, d, r, memory=memory)
