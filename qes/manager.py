"""
    author: Suhas Vittal
    date:   1 December 2023
"""

import networkx as nx

from collections import defaultdict

def concat(arr: list[int], delimiter=','):
    return delimiter.join(str(x) for x in arr)

class QesManager:
    def __init__(self, tanner_graph: nx.Graph):
        self.code = tanner_graph
        self.n = tanner_graph.number_of_nodes()

        # Public parameters:
        self.memory = 'z'

        # Simulation structures
        self.meas_ctr_map = {'curr': 0}
        self.event_ctr = 0
        self.obs_ctr = 0

        # Flag data structures:
        #   - flag_assignment_map[x] for some parity qubit x contains
        #       { data : flag } and an entry 'all' which contains all flags
        #       used in the syndrome extraction of x.
        self.flag_qubits = []
        self.flag_assignment_map = {}

    def add_flags_to(self, q1: int, q2: int, check: int) -> int:
        """
            Adds flags to qubits q1 and q2 during the measurement of the specified check.
            Returns the flag qubit id.
        """
        if check not in self.flag_assignment_map:
            self.flag_assignment_map[check] = {'all': []}
        # If everything checks out, update flag assignment map.
        self.flag_assignment_map[check]['all'].append(self.n)
        self.flag_assignment_map[check][q1] = self.n
        self.flag_assignment_map[check][q2] = self.n
        self.flag_qubits.append(self.n)
        self.n += 1

    def fopen(self, filename: str):
        self._writer = open(filename, 'w')

    def fclose(self):
        self._writer.close()

    def write_memory_experiment(self, rounds: int):
        data_qubits = self.code.graph['data_qubits']
        mem_checks = self.code.graph['checks'][self.memory]

        mem_flags = []
        for ch in self.flag_assignment_map:
            if self.code.nodes[ch]['node_type'] == self.memory:
                continue
            mem_flags.extend(self.flag_assignment_map[ch]['all'])
        #
        # PROLOGUE
        #
        self.big_comment('PROLOGUE')
        self.reset(data_qubits)
        if self.memory == 'x':
            self.h(data_qubits)
        #
        # BODY
        #
        self.big_comment('BODY')
        n_meas_per_round = len(mem_checks) + len(mem_flags)
        ectr = 0
        for r in range(rounds):
            prev_meas_ctr_map = self.meas_ctr_map.copy()
            self.annotation('timing_error')
            self._standard_syndrome_extraction()
            # Create detection events.
            for pq in mem_checks:
                if 'color' in self.code.nodes[pq]:
                    self.property('color', self.code.nodes[pq]['color'])
                base_ectr = ectr % n_meas_per_round
                if rounds > 1:
                    base_ectr += n_meas_per_round
                self.property('base', base_ectr)
                if r == 0:
                    m = self.meas_ctr_map[pq]
                    self.event(ectr, m)
                else:
                    m1, m2 = self.meas_ctr_map[pq], prev_meas_ctr_map[pq]
                    self.event(ectr, m1, m2)
                ectr += 1
            for fq in mem_flags:
                m = self.meas_ctr_map[fq]
                self.annotation('flag')
                self.event(ectr, m)
                ectr += 1
        #
        # EPILOGUE
        #
        self.big_comment('EPILOGUE')
        if self.memory == 'x':
            self.h(data_qubits)
        self.measure(data_qubits)
        # Create detection events and observables.
        for pq in mem_checks:
            support = self.code.nodes[pq]['schedule_order']
            meas_array = [self.meas_ctr_map[x] for x in support if x is not None]
            meas_array.append(self.meas_ctr_map[pq])
            base_ectr = ectr % n_meas_per_round
            if rounds > 1:
                base_ectr += n_meas_per_round
            self.property('base', base_ectr)
            if 'color' in self.code.nodes[pq]:
                self.property('color', self.code.nodes[pq]['color'])
            self.event(ectr, *meas_array)
            ectr += 1
        # Write observables.
        obs_list = self.code.graph['obs_list'][self.memory]
        for obs in obs_list:
            meas_array = [self.meas_ctr_map[x] for x in obs]
            self.auto_obs(*meas_array)

    def h(self, operands: list[int]):
        self._op('h', operands)

    def cx(self, operands: list[int]):
        self._op('cx', operands)

    def measure(self, operands: list[int], update_meas_ctr_map=True):
        self._op('measure', operands)
        k = self.meas_ctr_map['curr']
        for (i, x) in enumerate(operands):
            if update_meas_ctr_map:
                self.meas_ctr_map[x] = k+i
        self.meas_ctr_map['curr'] = k + len(operands)

    def reset(self, operands: list[int]):
        self._op('reset', operands)
    
    def event(self, event_no: int, *meas_operands: list[int]):
        self._op('event', [event_no, *meas_operands])

    def obs(self, obs_no: int, *meas_operands: list[int]):
        self._op('obs', [obs_no, *meas_operands])
    
    def auto_event(self, *meas_operands: list[int]):
        self._op('event', [self.event_ctr, *meas_operands])
        self.event_ctr += 1

    def auto_obs(self, *meas_operands: list[int]):
        self._op('obs', [self.obs_ctr, *meas_operands])
        self.obs_ctr += 1

    def _op(self, opname: str, operands: list[int]):
        if len(operands) == 0:
            return
        self._writer.write('%s %s;\n' % (opname, concat(operands)))

    def annotation(self, name: str):
        self._writer.write('@annotation %s\n' % name)

    def property(self, name: str, value: int|float):
        self._writer.write('@property %s %s\n' % (name, str(value)))

    def comment(self, text: str):
        self._writer.write('# %s\n' % text)

    def big_comment(self, text: str):
        self._writer.write('#\n# %s \n#\n' % text)

    def _standard_syndrome_extraction(self):
        checks = self.code.graph['checks']['all']
        x_checks = self.code.graph['checks']['x']

        self.reset(self.flag_qubits)
        self.reset(checks)

        self.h(x_checks)

        self.comment('FLAG SETUP')
        self._flag_hadamards(checks)
        self._flag_cnots(checks)

        self.comment('DATA CNOTS')
        self._data_cnots(checks)

        self.comment('FLAG TEARDOWN')
        self._flag_cnots(checks)
        self._flag_hadamards(checks)

        self.h(x_checks)
        self.measure(self.flag_qubits)
        self.measure(checks)

    def _flag_hadamards(self, checks: list[int]):
        depth = 0
        while True: 
            h_list = []
            for pq in checks:
                if pq not in self.flag_assignment_map:
                    continue
                flag_list = self.flag_assignment_map[pq]['all']
                if depth >= len(flag_list):
                    continue
                fq = flag_list[depth]
                s = self.code.nodes[pq]['node_type']
                if s == 'z':
                    h_list.append(fq)
            if len(h_list) == 0:
                break
            self.h(h_list)
            depth += 1

    def _flag_cnots(self, checks: list[int]):
        depth = 0
        while True:
            cx_list = []
            for pq in checks:
                if pq not in self.flag_assignment_map:
                    continue
                flag_list = self.flag_assignment_map[pq]['all']
                if depth >= len(flag_list):
                    continue
                fq = flag_list[depth]
                s = self.code.nodes[pq]['node_type']
                if s == 'x':
                    cx_list.extend([pq, fq])
                else:
                    cx_list.extend([fq, pq])
            if len(cx_list) == 0:
                break
            self.cx(cx_list)
            depth += 1

    def _data_cnots(self, checks: list[int]):
        depth = 0
        while True:
            cx_list = []
            for pq in checks:
                support = self.code.nodes[pq]['schedule_order']
                if depth >= len(support):
                    continue
                dq = support[depth]
                if dq is None:
                    continue
                cx_list.extend(self._get_cnot_targets(dq, pq))
            if len(cx_list) == 0:
                break
            self.cx(cx_list)
            depth += 1

    def _get_cnot_targets(self, dq: int, pq: int) -> tuple[int, int]:
        s = self.code.nodes[pq]['node_type']
        if pq not in self.flag_assignment_map or dq not in self.flag_assignment_map[pq]:
            fq = pq
        else:
            fq = self.flag_assignment_map[pq][dq]
        if s == 'x':
            return fq, dq
        else:
            return dq, fq

