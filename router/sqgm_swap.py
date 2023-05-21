# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2020.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
#
# SQGMSwap has been modified based on Qiskit 0.33.0's SabreSwap.

"""Routing via SWAP insertion using the SABRE method from Li et al."""

import logging
from collections import defaultdict
from copy import copy, deepcopy
import numpy as np

from qiskit.circuit.library.standard_gates import SwapGate
from qiskit.circuit.quantumregister import Qubit
from qiskit.transpiler.basepasses import TransformationPass
from qiskit.transpiler.exceptions import TranspilerError
from qiskit.transpiler.layout import Layout
from qiskit.dagcircuit import DAGOpNode

logger = logging.getLogger(__name__)

EXTENDED_SET_WEIGHT = 0.5  # Weight of lookahead window compared to front_layer
DECAY_RATE = 0.001  # Decay coefficient for penalizing serial swaps
DECAY_RESET_INTERVAL = 5  # How often to reset all decay rates to 1


class SQGMSwap(TransformationPass):

    def __init__(self, coupling_map, heuristic="basic", enhance_h=True, seed=None, fake_run=False):

        super().__init__()

        # Assume bidirectional couplings, fixing gate direction is easy later.
        if coupling_map is None or coupling_map.is_symmetric:
            self.coupling_map = coupling_map
        else:
            self.coupling_map = deepcopy(coupling_map)
            self.coupling_map.make_symmetric()

        self.heuristic = heuristic
        self.seed = seed
        self.fake_run = fake_run
        self.applied_predecessors = None
        self.qubits_decay = None
        self._bit_indices = None
        self.dist_matrix = None

        self.extended_set_size = self.coupling_map.size() # extended set size varies depending on physical circuit size
        self.enhance_h = enhance_h # if True, use H_sqgm (Eq 5), else use H_decay (Eq 4)
        self.pre_h_score = None # record the best H score of the previous action selection

    def run(self, dag):
        """Run SQGM on `dag`.

        Args:
            dag (DAGCircuit): the directed acyclic graph to be mapped.
        Returns:
            DAGCircuit: A dag mapped to be compatible with the coupling_map.
        Raises:
            TranspilerError: if the coupling map or the layout are not
            compatible with the DAG
        """
        if len(dag.qregs) != 1 or dag.qregs.get("q", None) is None:
            raise TranspilerError("Sabre swap runs on physical circuits only.")

        if len(dag.qubits) > self.coupling_map.size():
            raise TranspilerError("More virtual qubits exist than physical.")

        self.dist_matrix = self.coupling_map.distance_matrix

        rng = np.random.default_rng(self.seed)

        # Preserve input DAG's name, regs, wire_map, etc. but replace the graph.
        mapped_dag = None
        if not self.fake_run:  
            mapped_dag = dag._copy_circuit_metadata()

        canonical_register = dag.qregs["q"]
        current_layout = Layout.generate_trivial_layout(canonical_register)

        self._bit_indices = {bit: idx for idx, bit in enumerate(canonical_register)} # bit is a qubit

        # A decay factor for each qubit used to heuristically penalize recently
        # used qubits (to encourage parallelism).
        self.qubits_decay = {qubit: 1 for qubit in dag.qubits}
        
        # Introduce progress and gate buffer
        num_swap = 0 # num_swap is the number of swaps before we execute any gate
        phy_qubits = [current_layout._v2p[qubit] for qubit in dag.qubits] # list of qubit nums e.g. [0, 1, 2, 3, 4, 5]
        qct_progress = {idx: 0 for idx in phy_qubits} # store the progress of each qubit
        gate_storage = {idx: [] for idx in phy_qubits} # store the single qubit gates
                    
        # Start algorithm from the front layer and iterate until all gates done.
        num_search_steps = 0
        front_layer = dag.front_layer() # This front layer may contain 1qbg
        self.applied_predecessors = defaultdict(int)
        for _, input_node in dag.input_map.items(): # map from wire (reg,idx) to input nodes of the graph
            for successor in self._successors(input_node, dag): 
                self.applied_predecessors[successor] += 1

        while front_layer:
            execute_gate_list = []

            # Remove as many immediately applicable gates as possible
            for node in front_layer:

                if len(node.qargs) == 2:  
                    v0, v1 = node.qargs
                    # Accessing layout._v2p directly to avoid overhead from __getitem__ and a
                    # single access isn't feasible because the layout is updated on each iteration
                    if self.coupling_map.graph.has_edge(current_layout._v2p[v0], current_layout._v2p[v1]):
                        
                        if self.fake_run: continue
                        
                        # Release/apply 1qb gates in the buffer and realign the qct progress
                        for qubit in [current_layout._v2p[v0], current_layout._v2p[v1]]:
                            for qnode in gate_storage[qubit]:
                                self._apply_gate(mapped_dag, qnode, current_layout, canonical_register) # apply 1qb gate
                                qct_progress[qubit] += 1
                            gate_storage[qubit] = [] # reset

                        execute_gate_list.append(node)

                        p0, p1 = current_layout._v2p[v0], current_layout._v2p[v1]
                        temp_progress =  max(qct_progress[p0], qct_progress[p1]) + 1 # realign
                        qct_progress[p0] = temp_progress
                        qct_progress[p1] = temp_progress
                        num_swap = 0  # reset
                        
                else: # Single-qubit gates as well as barriers are free
                    
                    execute_gate_list.append(node)
                    v0 = node.qargs[0]
                    p0 = current_layout._v2p[v0]
                    num_swap = 0  # reset
                    gate_storage[p0].append(node)

            if execute_gate_list:

                self.pre_h_score = None # reset
                for node in execute_gate_list:

                    if self.fake_run or len(node.qargs) == 2: # only apply the single qubit gates on a fake run
                        self._apply_gate(mapped_dag, node, current_layout, canonical_register)

                    front_layer.remove(node)
                    for successor in self._successors(node, dag):
                        self.applied_predecessors[successor] += 1
                        if self._is_resolved(successor): 
                            front_layer.append(successor)

                    if node.qargs:
                        self._reset_qubits_decay()

                # Diagnostics
                logger.debug(
                    "free! %s",
                    [
                        (n.name if isinstance(n, DAGOpNode) else None, n.qargs)
                        for n in execute_gate_list
                    ],
                )
                logger.debug(
                    "front_layer: %s",
                    [(n.name if isinstance(n, DAGOpNode) else None, n.qargs) for n in front_layer],
                )
                continue
                
            # After all free gates are exhausted, heuristically find the best swap and insert it.
            # When two or more swaps tie for best score, pick one randomly.

            extended_set = self._obtain_extended_set(dag, front_layer)
            swap_candidates = self._obtain_swaps(front_layer, current_layout)
            sqgm_swap_scores = dict.fromkeys(swap_candidates, 0) # sqgm H score
            sabre_swap_scores = dict.fromkeys(swap_candidates, 0) # sabre H score
            
            for swap_qubits in sqgm_swap_scores:
                trial_layout = current_layout.copy()
                trial_layout.swap(*swap_qubits)
                score = self._score_heuristic(
                    self.heuristic, front_layer, extended_set, trial_layout, swap_qubits
                )
                sabre_swap_scores[swap_qubits] = score # sabre H score
                if self.enhance_h:
                    p0, p1 = current_layout[swap_qubits[0]], current_layout[swap_qubits[1]]
                    diff = max(qct_progress[p0], qct_progress[p1])
                    score += diff/self.coupling_map.size()
                    sqgm_swap_scores[swap_qubits] = score # sqgm H score
            
            if self.enhance_h:
                sqgm_min_score = min(sqgm_swap_scores.values()) # sqgm best score
            sabre_min_score = min(sabre_swap_scores.values()) # sabre best score

            if self.enhance_h and (self.pre_h_score is None or sqgm_min_score < self.pre_h_score):
                self.pre_h_score = sqgm_min_score
                best_swaps = [k for k, v in sqgm_swap_scores.items() if v == sqgm_min_score]   
            else: # if sqgm best score worsens then we retreat to sabre best score
                best_swaps = [k for k, v in sabre_swap_scores.items() if v == sabre_min_score]

            best_swaps.sort(key=lambda x: (self._bit_indices[x[0]], self._bit_indices[x[1]]))
            best_swap = rng.choice(best_swaps)

            # use 1qb gates in gate_storage to make the progress as even as possible
            p0, p1 = current_layout._v2p[best_swap[0]], current_layout._v2p[best_swap[1]]
            if not self.fake_run:
                while qct_progress[p0] < qct_progress[p1] and len(gate_storage[p0]) > 0:
                    self._apply_gate(mapped_dag, gate_storage[p0].pop(0), current_layout, canonical_register) # apply the single qubit gate
                    qct_progress[p0] += 1
                while qct_progress[p1] < qct_progress[p0] and len(gate_storage[p1]) > 0:
                    self._apply_gate(mapped_dag, gate_storage[p1].pop(0), current_layout, canonical_register) # apply the single qubit gate
                    qct_progress[p1] += 1
            

            swap_node = DAGOpNode(op=SwapGate(), qargs=best_swap)
            self._apply_gate(mapped_dag, swap_node, current_layout, canonical_register)

            # swap the gate_storages
            gate_storage[p0], gate_storage[p1] = gate_storage[p1], gate_storage[p0] # swap the 2 qubit storages to reflect the swap
            num_swap += 1
            if num_swap > self.coupling_map.size():
                print('@@', num_search_steps, num_swap, qct_progress)
                raise TranspilerError('Fallback is necessary!')
                
            current_layout.swap(*best_swap)
            
            # realign qct_progress after apply the selected swap
            v0, v1 = best_swap
            p0, p1 = current_layout._v2p[v0], current_layout._v2p[v1]
            temp_progress =  max(qct_progress[p0], qct_progress[p1]) + 3
            qct_progress[p0] = temp_progress
            qct_progress[p1] = temp_progress

            num_search_steps += 1
            if num_search_steps % DECAY_RESET_INTERVAL == 0:
                self._reset_qubits_decay()
            else:
                self.qubits_decay[best_swap[0]] += DECAY_RATE
                self.qubits_decay[best_swap[1]] += DECAY_RATE

            # Diagnostics
            logger.debug("SWAP Selection...")
            logger.debug("extended_set: %s", [(n.name, n.qargs) for n in extended_set])
            logger.debug("swap scores: %s", sqgm_swap_scores if self.enhance_h else sabre_swap_scores)
            logger.debug("best swap: %s", best_swap)
            logger.debug("qubits decay: %s", self.qubits_decay)
            
        self.property_set["final_layout"] = current_layout
        
        if not self.fake_run:
            # Apply 1qb gates left in the gate_storage to finish the transformation
            for qubit in phy_qubits:
                for node in gate_storage[qubit]:
                    self._apply_gate(mapped_dag, node, current_layout, canonical_register) # apply the remaining gates
            return mapped_dag

        return dag

    def _apply_gate(self, mapped_dag, node, current_layout, canonical_register):
        if self.fake_run:
            return
        new_node = _transform_gate_for_layout(node, current_layout, canonical_register)
        mapped_dag.apply_operation_back(new_node.op, new_node.qargs, new_node.cargs)

    def _reset_qubits_decay(self):
        """Reset all qubit decay factors to 1 upon request (to forget about
        past penalizations).
        """
        self.qubits_decay = {k: 1 for k in self.qubits_decay.keys()}
        
    def _successors(self, node, dag):
        for _, successor, edge_data in dag.edges(node):
            if not isinstance(successor, DAGOpNode):
                continue
            if isinstance(edge_data, Qubit):
                yield successor

    def _is_resolved(self, node):
        """Return True if all of a node's predecessors in dag are applied."""
        return self.applied_predecessors[node] == len(node.qargs)

    def _obtain_extended_set(self, dag, front_layer):
        """Populate extended_set by looking ahead a fixed number of gates.
        For each existing element add a successor until reaching limit.
        """
        extended_set = []
        incremented = []
        tmp_front_layer = front_layer
        done = False
        while tmp_front_layer and not done:
            new_tmp_front_layer = []
            for node in tmp_front_layer:
                for successor in self._successors(node, dag):
                    incremented.append(successor)
                    self.applied_predecessors[successor] += 1  
                    if self._is_resolved(successor):
                        new_tmp_front_layer.append(successor)
                        if len(successor.qargs) == 2:
                            extended_set.append(successor)
                if len(extended_set) >= self.extended_set_size: 
                    done = True
                    break
            tmp_front_layer = new_tmp_front_layer
        for node in incremented:
            self.applied_predecessors[node] -= 1 
        return extended_set

    def _obtain_swaps(self, front_layer, current_layout):
        """Return a set of candidate swaps that affect qubits in front_layer.

        For each virtual qubit in front_layer, find its current location
        on hardware and the physical qubits in that neighborhood. Every SWAP
        on virtual qubits that corresponds to one of those physical couplings
        is a candidate SWAP.

        Candidate swaps are sorted so SWAP(i,j) and SWAP(j,i) are not duplicated.
        """
        candidate_swaps = set()
        for node in front_layer:
            for virtual in node.qargs: 
                physical = current_layout[virtual]
                for neighbor in self.coupling_map.neighbors(physical):
                    virtual_neighbor = current_layout[neighbor]
                    swap = sorted([virtual, virtual_neighbor], key=lambda q: self._bit_indices[q])
                    candidate_swaps.add(tuple(swap))

        return candidate_swaps

    def _compute_cost(self, layer, layout):
        cost = 0
        layout_map = layout._v2p
        for node in layer:
            if len(node.qargs) == 1: 
                raise TranspilerError('The layer contains no 1-qubit gates!', node.qargs, len(node.qargs))

            cost += self.dist_matrix[layout_map[node.qargs[0]], layout_map[node.qargs[1]]]
        return cost

    def _score_heuristic(self, heuristic, front_layer, extended_set, layout, swap_qubits=None):
        """Return a heuristic score for a trial layout.

        Assuming a trial layout has resulted from a SWAP, we now assign a cost
        to it. The goodness of a layout is evaluated based on how viable it makes
        the remaining virtual gates that must be applied.
        """
        first_cost = self._compute_cost(front_layer, layout)
        if heuristic == "basic":
            return first_cost

        first_cost /= len(front_layer)         
        second_cost = 0
        if extended_set:
            second_cost = self._compute_cost(extended_set, layout) / len(extended_set)

        total_cost = first_cost + EXTENDED_SET_WEIGHT * second_cost
        if heuristic == "lookahead":
            return total_cost

        if heuristic == "decay":
            return (
                max(self.qubits_decay[swap_qubits[0]], self.qubits_decay[swap_qubits[1]])
                * total_cost
            )

        raise TranspilerError("Heuristic %s not recognized." % heuristic)

def _transform_gate_for_layout(op_node, layout, device_qreg):
    """Return node implementing a virtual op on given layout."""
    mapped_op_node = copy(op_node)

    premap_qargs = op_node.qargs
    mapped_qargs = map(lambda x: device_qreg[layout._v2p[x]], premap_qargs)
    mapped_op_node.qargs = list(mapped_qargs)

    return mapped_op_node