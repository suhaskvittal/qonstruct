import stim

def add_operator(circuit, arr, is_x):
    targets = []
    for (i, q) in enumerate(arr):
        targets.append(stim.target_x(q) if is_x else stim.target_z(q))
        if i < len(arr)-1:
            targets.append(stim.target_combiner())
    circuit.append("MPP", targets)

def measure_checks(circuit, gr):
    for ch in gr.graph['checks']['all']:
        add_operator(circuit, list(gr.neighbors(ch)), gr.nodes[ch]['node_type'] == 'x')

def measure_logical_qubits(circuit, gr, is_x):
    typ = 'x' if is_x else 'z'
    for (i, obs) in enumerate(gr.graph['obs_list'][typ]):
        add_operator(circuit, obs, is_x)
        circuit.append('OBSERVABLE_INCLUDE', [stim.target_rec(-1)], i)

def compute_distance(gr, is_x):
    circuit = stim.Circuit()
    measure_logical_qubits(circuit, gr, is_x)
    measure_checks(circuit, gr)
    # Inject error on the data qubits.
    if is_x:
        circuit.append('Z_ERROR', list(gr.graph['data_qubits']), 0.1)
    else:
        circuit.append('X_ERROR', list(gr.graph['data_qubits']), 0.1)
    # Measure checks and then the logical qubits.
    measure_checks(circuit, gr)
    n_checks = len(gr.graph['checks']['all'])
    t = 'x' if is_x else 'z'
    for (i, ch) in enumerate(gr.graph['checks']['all']):
        if gr.nodes[ch]['node_type'] == t:
            circuit.append('DETECTOR', [stim.target_rec(-(n_checks-i)), stim.target_rec(-(2*n_checks-i))])
    measure_logical_qubits(circuit, gr, is_x)
    errors = circuit.search_for_undetectable_logical_errors(
                dont_explore_edges_increasing_symptom_degree=False,
                dont_explore_detection_event_sets_with_size_above=3,
                dont_explore_edges_with_degree_above=3)
    return len(errors)