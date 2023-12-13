import csv
import json
import os
import sys
import time
import qiskit

from collections import defaultdict
from numpy import Infinity
from statistics import mean
from typing import List, Tuple, Union
from util import ARCHGRAPHS, ROUTERS

from qiskit import QuantumCircuit
from qiskit.circuit.library.standard_gates import SwapGate
from qiskit.converters import circuit_to_dag
from qiskit.quantum_info.synthesis.two_qubit_decompose import TwoQubitBasisDecomposer
from qiskit.transpiler import CouplingMap, PassManager
from qiskit.transpiler.passes import (
    OptimizeSwapBeforeMeasure,
    RemoveDiagonalGatesBeforeMeasure,
    FullAncillaAllocation,
    EnlargeWithAncilla,
    ApplyLayout,
    BarrierBeforeFinalMeasurements,
    Unroll3qOrMore,
    Optimize1qGatesDecomposition,
    RemoveResetInZeroState,
    Decompose,
    Depth,
    FixedPoint,
    CommutativeCancellation
)
from qiskit.transpiler.passes.synthesis.unitary_synthesis import (
    _choose_euler_basis,
    _choose_kak_gate
)

'''ARCHITECTURE GRAPHS'''
from ag import *

'''MAPPERS'''
from qiskit.transpiler.passes import SabreLayout                # SABRE39
# from mapper.sabre0330_layout import SabreLayout               # SABRE33

'''ROUTERS'''
from qiskit.transpiler.passes import SabreSwap as router_0394   # SABRE39
from router.sabre0330_swap import SabreSwap as router_0330      # SABRE33
from router.sqgm_swap import SQGMSwap as router_sqgm            # SQGM
from router.nassc_swap import NASSCSwap as router_nassc         # NASSC
from qiskit.transpiler.passes import StochasticSwap as router_stoch

'''CPU COUNT'''
from qiskit.tools.parallel import CPU_COUNT


def get_non_single_qg_names(non_single_qg):
    """From the provided non-single-qubit gates, get all of the gate names.

    Returns:
        A list of the names of all non-single-qubit gates.
    """
    non_single_qg_names = set()
    for gate in non_single_qg:
        if gate.name == "cx":
            continue
        non_single_qg_names.add(gate.name)
    return list(non_single_qg_names)


def initialise_routers(coupling_map, base_routers, basis_gates, heuristics, inc_cc, inc_naive_h):
    """Initialise routers, including +CC and naive-H variations.

    Returns:
        A dictionary of the initialised routers.
    """
    routers = dict()

    for name in base_routers:

        Router = ROUTERS[name]
        for h in heuristics:

            if Router == router_nassc:
                routers[f"{name}-{h}"] = Router(
                    coupling_map,
                    heuristic=h,
                    enable_factor_block=True,
                    enable_factor_commute_0=True,
                    enable_factor_commute_1=True,
                    factor_block=1,
                    factor_commute_0=1,
                    factor_commute_1=1,
                    decomposer2q=TwoQubitBasisDecomposer(
                        gate=_choose_kak_gate(basis_gates),
                        euler_basis=_choose_euler_basis(basis_gates)
                    ),
                    approximation_degree=None
                )
            elif Router == router_stoch:
                routers[f"{name}"] = Router(
                    coupling_map
                )
            else:
                routers[f"{name}-{h}"] = Router(
                    coupling_map,
                    heuristic=h
                )
                # Include CC or naive H if required
                if inc_cc:
                    routers[f"{name}+cc-{h}"] = Router(
                        coupling_map,
                        heuristic=h
                    )
                if Router == router_sqgm and inc_naive_h:
                    routers[f"{name}+nh-{h}"] = Router(
                        coupling_map,
                        heuristic=h,
                        enhance_h=False
                    )
    
    return routers


def initialise_circuit(filepath, verbose):
    """Initialise quantum circuit, and decompose if needed.

    Returns:
        The initialised (and decomposed, if needed) circuit, along with
        the starting and ending time of the decomposition (if any).
    """
    # Initialise quantum circuit and extract metrics
    cir = QuantumCircuit.from_qasm_file(filepath)
    dag = circuit_to_dag(cir)
    two_qg_list = dag.two_qubit_ops()
    mul_qg_list = dag.multi_qubit_ops()
    non_single_qg_names = get_non_single_qg_names(two_qg_list + mul_qg_list)
    if verbose:
        if non_single_qg_names:
            print(f"GATES TO DECOMPOSE: {non_single_qg_names}")
        else:
            print("NO GATES TO DECOMPOSE")
    
    init_pm = PassManager()
    init_pm.append([
        Unroll3qOrMore(),
        RemoveResetInZeroState(),
        OptimizeSwapBeforeMeasure(),
        RemoveDiagonalGatesBeforeMeasure()
    ])
    if non_single_qg_names: # If there are non-single-qubit gates which are not CNOTs, decompose those gates before routing
        init_pm.append(Decompose(non_single_qg_names))
    
    init_start = time.time()
    init_cir = init_pm.run(cir)
    init_end = time.time()

    return init_cir, init_start, init_end


def build_routing_passes(name, Router, basis_gates, verbose):
    """Build the appropriate router passes for Router.
    """
    router_pm = PassManager()

    # Add pre-routing optimisations and router
    if isinstance(Router, router_nassc):
        if verbose:
            print("Enabling pre-routing optimisation")
        router_pm.append([
            Optimize1qGatesDecomposition(basis=basis_gates),
            BarrierBeforeFinalMeasurements(),
            Router,
            CommutativeCancellation(basis_gates=basis_gates),
            RemoveResetInZeroState()
        ])
    else:
        if verbose:
            print("No pre-routing optimisation")
        router_pm.append(Router)
    
    # Add decomposition after routing
    router_pm.append(Decompose(SwapGate))
    
    # Add post-routing optimisations
    if isinstance(Router, router_nassc):
        if verbose:
            print("Enabling post-routing optimisation")
        router_pm.append([
            Depth(),
            FixedPoint("depth"),
            Optimize1qGatesDecomposition(basis=basis_gates),
            CommutativeCancellation(basis_gates=basis_gates)
            ], do_while=lambda property_set: not property_set["depth_fixed_point"]
        )
    elif name.split('-')[0].endswith("cc"):
        if verbose:
            print("Enabling simple commutative cancellation")
        router_pm.append([
            Depth(),
            FixedPoint("depth"),
            CommutativeCancellation()
            ], do_while=lambda property_set: not property_set["depth_fixed_point"]
        )
    else:
        # Do nothing
        if verbose:
            print("No post-routing optimisation")

    return router_pm


def run(
        filepath: str,
        arch_graph: str,
        base_routers: Tuple[str],
        reps: int = 5,
        heuristics: Tuple[str] = ("decay",),
        mode: str = "best",
        objective: str = "depth",
        init_mapping: QuantumCircuit = None,
        inc_naive_h: bool = False,
        inc_cc: bool = False,
        verbose: bool = False
    ) -> Tuple[
        Tuple[Union[router_0394, router_sqgm, router_nassc]],
        int,
        Union[int,List[int]],
        int,
        Union[int,List[int]],
        float,
        QuantumCircuit,
        QuantumCircuit
    ]:
    """Run experiments on provided routers and configurations.
    
    If any error is raised during a run, simply rerun.
    
    Args:
        `filepath` (str): The path to the logical circuit to be transformed.
        `arch_graph`: The target backend.
        `reps` (int): The number of reps.
        `base_routers` (tuple): The routers to be run.
        `heuristics` (str): Heuristic(s) to be used for each router.
        `mode` (int): Whether to record the best result (0) or every result
            (1) after specified reps.
        `objective` (str): Whether to determine the best result based on
            "[gate] size" or "[circuit] depth".
        `init_mapping` (QuantumCircuit): A well-decomposed mapped circuit to run
            the routers on. If None, a mapping will be generated using SabreLayout.
        `inc_cc` (bool): If True, each base router will be accompanied with
            a "+cc" variation with post-routing commutative cancellation.
        `inc_naive_h` (bool): If True, SQGM will will be accompanied with a
            "+nh" variation which uses the naive H_decay instead of H_sqgm.
        `verbose` (bool): If True, debugging messages are displayed.
    Returns:
        `routers` (tuple): A list of routers, including both the base
            routers and any of their variations which may have been created.
        `num_2qbg_in` (int): The number of 2-qubit gates in the input circuit.
        `num_2qbg_out` (Union[int,List[int]]): The number of 2-qubit gates in
            the output circuit(s).
        `depth_in` (int): The depth of the input circuit.
        `depth_out` (Union[Int,List[int]]): The depth of the output circuit(s).
        `time_used` (float): The running time(s) of each run, which includes both
            initial mapping generation (if applicable) and routing.
        `input_cir` (QuantumCircuit): The input circuit after initial mapping
            has been generated and applied (if applicable).
        `outpit_cir` (QuantumCircuit): The output circuit following routing.
    """
    # Initialise passers
    coupling_map = CouplingMap(couplinglist=ARCHGRAPHS[arch_graph].edges())
    coupling_map.make_symmetric()
    try: # ensuring compatibility with older versions of Qiskit
        Mapper = SabreLayout(coupling_map, skip_routing=True)
    except TypeError:
        Mapper = SabreLayout(coupling_map)

    # For NASSCSwap
    basis_gates = ["x", "h", "u", "u1", "u2", "u3", "cx"]

    # Initialise router and result containers
    routers = initialise_routers(
        coupling_map,
        base_routers,
        basis_gates,
        heuristics,
        inc_cc,
        inc_naive_h
    )
    if mode == "best":
        depth_out = defaultdict(lambda: Infinity)
        num_2qg_out = defaultdict(lambda: Infinity)
        time_used = defaultdict(float)
        input_cir = dict()
        output_cir = dict()
    if mode == "all":
        depth_out = defaultdict(list)
        num_2qg_out = defaultdict(list)
        time_used = defaultdict(list)
        input_cir = defaultdict(list)
        output_cir = defaultdict(list)

    # Initialise quantum circuit
    init_cir, init_start, init_end = initialise_circuit(filepath, verbose)

    num_2qg_in = len(circuit_to_dag(init_cir).two_qubit_ops())
    depth_in = init_cir.depth()

    # If a well-decomposed initial mapping has been supplied, we use that instead
    if init_mapping:
        if len(circuit_to_dag(init_mapping).two_qubit_ops()) < num_2qg_in:
            raise ValueError("Supplied mapped circuit is not well-decomposed.")
        if verbose:
            print("USING SUPPLIED INITIAL MAPPING")
        init_cir = init_mapping

    cir_in_prev = None

    i = 1
    # Run the experiment for provided reps
    while i <= reps:
        print(i, end='\r')
        try:
            if init_mapping:
                
                map_start = 0
                map_end = 0
                cir_in = init_cir
            
            else:
                # Reset seed to ensure a different layout (initial mapping) is generated for every rep
                Mapper.seed = None

                # Perform and run mapping passes
                mapper_pm = PassManager()
                mapper_pm.append([
                    Mapper,
                    FullAncillaAllocation(coupling_map),
                    EnlargeWithAncilla(),
                    ApplyLayout()
                ])
                map_start = time.time()
                cir_in = mapper_pm.run(init_cir)
                map_end = time.time()

                # Compare current with previous initial mapping
                if cir_in_prev is not None:
                    assert cir_in != cir_in_prev, "Layout has not been reset!"
                cir_in_prev = cir_in
            
            # Within each rep, run each router on the same initial mapping
            for name, Router in routers.items():

                # Build routing passes
                router_pm = build_routing_passes(name, Router, basis_gates, verbose)
                
                # Begin routing
                route_start = time.time()
                cir_out = router_pm.run(cir_in)
                route_end = time.time()

                # Collect results
                dag_out = circuit_to_dag(cir_out)
                assert not dag_out.named_nodes(SwapGate), "Not all SWAPs have been decomposed!"
                depth = dag_out.depth()
                num_2qg = len(dag_out.two_qubit_ops())

                if verbose:
                    print(f"[{i}] {name}")
                    print(f"    depth: {depth}")
                    print(f"    num_2qg: {num_2qg}")
                    print(f"    mapper.seed: {Mapper.seed}")
                    print(f"    router.seed: {Router.seed}")
                
                run_time = (init_end - init_start) + (map_end - map_start) + (route_end - route_start)
                if mode == "best": # Record best rep
                    if objective == "depth" and depth >= depth_out[name]:
                        continue
                    if objective == "size" and num_2qg >= num_2qg_out[name]:
                        continue
                    depth_out[name] = depth
                    num_2qg_out[name] = num_2qg
                    time_used[name] = run_time
                    input_cir[name] = cir_in
                    output_cir[name] = cir_out
                if mode == "all": # Record all reps
                    depth_out[name].append(depth)
                    num_2qg_out[name].append(num_2qg)
                    time_used[name].append(run_time)
                    input_cir[name].append(cir_in)
                    output_cir[name].append(cir_out)
            i += 1
        
        except Exception as e: # If any error is raised during routing, simply re-run
            print(e)
            continue
    
    return routers, num_2qg_in, num_2qg_out, depth_in, depth_out, time_used, input_cir, output_cir, init_cir


def export_qasm(
        qasm_path: str,
        objective: str,
        mode: str,
        routers: dict[str, Union[router_0394, router_sqgm, router_nassc]],
        input_cirs: defaultdict,
        output_cirs: defaultdict,
        init_cirs: defaultdict,
        num_2qg_outs: defaultdict,
        depth_outs: defaultdict,
    ) -> None:

    for filename in num_2qg_outs.keys():
        init_cirs[filename].qasm(filename=f"{qasm_path}{filename.replace('.qasm', '')}-init.qasm")

        for router_name in routers.keys():
            
            if mode == "best": # Record best rep
                if objective == "size":
                    min_num_2qg_cir_in = input_cirs[filename][router_name]
                    min_num_2qg_cir_out = output_cirs[filename][router_name]
                    min_num_2qg_cir_in.qasm(filename=f"{qasm_path}{filename.replace('.qasm', '')}-{router_name.split('-')[0]}-2qg_min-in.qasm")
                    min_num_2qg_cir_out.qasm(filename=f"{qasm_path}{filename.replace('.qasm', '')}-{router_name.split('-')[0]}-2qg_min-out.qasm")
                if objective == "depth":
                    min_depth_cir_in = input_cirs[filename][router_name]
                    min_depth_cir_out = output_cirs[filename][router_name]
                    min_depth_cir_in.qasm(filename=f"{qasm_path}{filename.replace('.qasm', '')}-{router_name.split('-')[0]}-depth_min-in.qasm")
                    min_depth_cir_out.qasm(filename=f"{qasm_path}{filename.replace('.qasm', '')}-{router_name.split('-')[0]}-depth_min-out.qasm")
            
            if mode == "all": # Record all reps
                if objective == "size":
                    min_num_2qg_idx = num_2qg_outs[filename][router_name].index(min(num_2qg_outs[filename][router_name]))
                    min_num_2qg_cir_in = input_cirs[filename][router_name][min_num_2qg_idx]
                    min_num_2qg_cir_out = output_cirs[filename][router_name][min_num_2qg_idx]
                    min_num_2qg_cir_in.qasm(filename=f"{qasm_path}{filename.replace('.qasm', '')}-{router_name.split('-')[0]}-2qg_min-in.qasm")
                    min_num_2qg_cir_out.qasm(filename=f"{qasm_path}{filename.replace('.qasm', '')}-{router_name.split('-')[0]}-2qg_min-out.qasm")
                if objective == "depth":
                    min_depth_idx = depth_outs[filename][router_name].index(min(depth_outs[filename][router_name]))
                    min_depth_cir_in = input_cirs[filename][router_name][min_depth_idx]
                    min_depth_cir_out = output_cirs[filename][router_name][min_depth_idx]
                    min_depth_cir_in.qasm(filename=f"{qasm_path}{filename.replace('.qasm', '')}-{router_name.split('-')[0]}-depth_min-in.qasm")
                    min_depth_cir_out.qasm(filename=f"{qasm_path}{filename.replace('.qasm', '')}-{router_name.split('-')[0]}-depth_min-out.qasm")


def export_csv(
        csv_path: str,
        heuristics: Tuple[str],
        reps: int,
        mode: str,
        rep_splits: Tuple[int],
        routers: dict[str, Union[router_0394, router_sqgm, router_nassc]],
        num_2qg_ins: defaultdict,
        depth_ins: defaultdict,
        num_2qg_outs: defaultdict,
        depth_outs: defaultdict,
        times_used: defaultdict
    ) -> None:

    if not rep_splits:
        rep_splits = tuple([n for n in range(1, reps + 1)])
    
    with open(csv_path, "w", newline='') as file:

        writer = csv.writer(file)
        colnames = ['circuit']
        colnames.extend(['2qg_in', 'depth_in'])
        
        for router_name in routers.keys():
            
            if len(heuristics) == 1:
                prefix = router_name.split('-')[0]
            if len(heuristics) > 1:
                match router_name.split('-'):
                    case router, 'basic':
                        prefix = f'{router}-bs'
                    case router, 'lookahead':
                        prefix = f'{router}-la'
                    case router, 'decay':
                        prefix = f'{router}-dc'
            
            if mode == "best":
                colnames.append(f'{prefix}_min_2qg_out')
                colnames.append(f'{prefix}_min_depth_out')
                colnames.append(f'{prefix}_min_2qg_ratio')
                colnames.append(f'{prefix}_min_depth_ratio')
                colnames.append(f'{prefix}_avg_time')
            if mode == "all":
                colnames.extend([f'{prefix}_2qg_out_{i}' for i in rep_splits])
                colnames.extend([f'{prefix}_depth_out_{i}' for i in rep_splits])
                colnames.append(f'{prefix}_avg_time') # Only average time is recorded
        
        writer.writerow(colnames)
        for filename in num_2qg_ins.keys():

            record = [filename]
            record.extend([num_2qg_ins[filename], depth_ins[filename]])
            
            for router_name in routers.keys():
                
                if mode == "best":
                    record.append(num_2qg_outs[filename][router_name])
                    record.append(depth_outs[filename][router_name])
                    record.append(round(num_2qg_outs[filename][router_name]/num_2qg_ins[filename], 3))
                    record.append(round(depth_outs[filename][router_name]/depth_ins[filename], 3))
                    record.append(round(times_used[filename][router_name], 2))
                if mode == "all":
                    record.extend([min(num_2qg_outs[filename][router_name][:i]) for i in rep_splits])
                    record.extend([min(depth_outs[filename][router_name][:i]) for i in rep_splits])
                    record.append(round(mean(times_used[filename][router_name]), 2))
            
            writer.writerow(record)


def main(
        bench_path: str,
        bench_filter: str,
        qasm_path: str,
        csv_path: str,
        arch_graph: str,
        reps: int,
        objective: str,
        mode: str,
        rep_splits: Tuple[int],
        base_routers: Tuple[str],
        heuristics: Tuple[str],
        init_mapping: QuantumCircuit,
        inc_cc: bool,
        inc_naive_h: bool,
        exp_qasm: bool,
        exp_csv: bool,
        verbose: bool,
    ) -> None:
    """Run experiments and export results.

    Args:
        `verbose` (bool): Toggle display detailed experiment metrics and information.
        `export` (bool): Toggle export of the results to qasm and csv files.
    """
    
    '''CHECK FOR INVALID SETTINGS'''
    if reps < 0:
        raise ValueError("`reps` must be a non-zero positive integer")
    if objective not in ("size", "depth"):
        raise ValueError("`objective` must be either 'size' or 'depth'")
    if mode not in ("best", "all"):
        raise ValueError("`mode` must be either 'best' or 'all'")
    for heuristic in heuristics:
        if heuristic not in ("basic", "lookahead", "decay"):
            raise ValueError("`heuristic` must be one of the following: 'basic', 'lookahead', 'decay'")
    
    '''INITIALISE RESULT CONTAINERS'''
    num_2qg_ins = defaultdict(int)
    num_2qg_outs = defaultdict(dict)
    depth_ins = defaultdict(int)
    depth_outs = defaultdict(dict)
    times_used = defaultdict(dict)
    input_cirs = defaultdict(dict)
    output_cirs = defaultdict(dict)
    init_cirs = defaultdict(QuantumCircuit)

    print("+-----------------------+")
    print(f"BACKEND\t\t{arch_graph}")
    print(f"BENCHMARK\t{bench_path}")
    print(f"REPS\t\t{reps}")
    print(f"OBJECTIVE\t{objective}")
    print(f"MODE\t\t{'record best result' if mode == 'best' else 'record all results'}")
    print(f"ROUTERS\t\t{', '.join(base_routers)}")
    print(f"HEURISTICS\t{', '.join(heuristics)}")
    print(f"INC CC?\t\t{'yes' if inc_cc else 'no'}")
    print(f"INC NAIVE H?\t{'yes' if inc_naive_h else 'no'}")
    print(f"EXP QASM?\t{'yes' if exp_qasm else 'no'}")
    print(f"EXP CSV?\t{'yes' if exp_csv else 'no'}")
    print("+-----------------------+")
    
    '''RUN EXPERIMENT'''
    for filename in os.listdir(bench_path):

        if not filename.endswith(".qasm"): continue
        if bench_filter not in filename: continue
        print(f"[{filename}]")
        
        routers, num_2qg_in, num_2qg_out, depth_in, depth_out, time_used, input_cir, output_cir, init_cir = run(
            bench_path+filename,
            arch_graph,
            base_routers,
            reps,
            heuristics,
            mode,
            objective,
            init_mapping,
            inc_naive_h,
            inc_cc,
            verbose
        )
        for router_name in routers.keys():
            num_2qg_outs[filename][router_name] = num_2qg_out[router_name]
            depth_outs[filename][router_name] = depth_out[router_name]
            times_used[filename][router_name] = time_used[router_name]
            input_cirs[filename][router_name] = input_cir[router_name]
            output_cirs[filename][router_name] = output_cir[router_name]

        num_2qg_ins[filename] = num_2qg_in
        depth_ins[filename] = depth_in
        init_cirs[filename] = init_cir
        
        for router_name in routers.keys():

            print(router_name)
            if mode == "best": # Record best rep
                print(f"    2qg out/in ratio: {round(num_2qg_outs[filename][router_name]/num_2qg_in, 3)}")
                print(f"    depth out/in ratio: {round(depth_outs[filename][router_name]/depth_in, 3)}")
                print(f"    used time (s): {round(times_used[filename][router_name], 2)}")
            if mode == "all": # Record all reps
                print(f"    2qg out/in ratios: {', '.join([str(round(num_2qg_/num_2qg_in, 3)) for num_2qg_ in num_2qg_outs[filename][router_name]])}")
                print(f"    depth out/in ratios: {', '.join([str(round(depth_/depth_in, 3)) for depth_ in depth_outs[filename][router_name]])}")
                print(f"    used times (s): {', '.join([str(round(time_, 2)) for time_ in times_used[filename][router_name]])}")

        print("+-----------------------+")
    
    if mode == "best": # Record best rep
        total_num_2qg_in = sum(num_2qg_ins.values())
        total_depth_in = sum(depth_ins.values())
        for router_name in routers.keys():
            print(router_name)
            total_num_2qg_out = sum([num_2qg_outs[filename][router_name] for filename in num_2qg_outs.keys()])
            total_depth_out = sum([depth_outs[filename][router_name] for filename in depth_outs.keys()])
            avg_time_used = mean([times_used[filename][router_name] for filename in times_used.keys()])
            print(f"    average 2qg out/in ratio: {round(total_num_2qg_out/total_num_2qg_in, 3)}")
            print(f"    average depth out/in ratio: {round(total_depth_out/total_depth_in, 3)}")
            print(f"    average used time (s): {round(avg_time_used, 2)}")
    if mode == "all": # Record all reps
        for filename in num_2qg_ins.keys():
            print(f'[{filename}]')
            for router_name in routers.keys():
                print(router_name)
                print(f"    average 2qg out/in ratio: {round(mean(num_2qg_outs[filename][router_name])/num_2qg_ins[filename], 3)}")
                print(f"    average depth out/in ratio: {round(mean(depth_outs[filename][router_name])/depth_ins[filename], 3)}")
                print(f"    average used time (s): {round(mean(times_used[filename][router_name]), 2)}")

    if exp_qasm:
        export_qasm(
            qasm_path,
            objective,
            mode,
            routers,
            input_cirs,
            output_cirs,
            init_cirs,
            num_2qg_outs,
            depth_outs
        )
    if exp_csv:
        export_csv(
            csv_path,
            heuristics,
            reps,
            mode,
            rep_splits,
            routers,
            num_2qg_ins,
            depth_ins,
            num_2qg_outs,
            depth_outs,
            times_used
        )


if __name__ == "__main__":

    if len(sys.argv) < 2 or 4 < len(sys.argv):
        print("Expected: python main.py <test_config.json> [init-mapping.qasm] [--verbose or -v]")
        sys.exit(1)
    
    verbose = False
    if len(sys.argv) > 2 and sys.argv[-1] in ("--verbose", "-v"):
        verbose = True
    elif len(sys.argv) == 4 and sys.argv[3] not in ("--verbose", "-v"):
        print("Only supported flags are '--verbose' or '-v'")
        sys.exit(1)

    exp = open(sys.argv[1])
    config = json.load(exp)
    exp.close()

    if len(sys.argv) == 4:
        init_mapping = QuantumCircuit.from_qasm_file(sys.argv[2])
    else:
        init_mapping = None

    print("+-----------------------+")
    print(f"Qiskit {qiskit.__version__}")
    print(f"CPU count: {CPU_COUNT}")

    main(
        bench_path=config["bench_path"],
        bench_filter=config["bench_filter"],
        qasm_path=config["qasm_path"],
        csv_path=config["csv_path"],
        arch_graph=config["ag"],
        reps=config["reps"],
        objective=config["objective"],
        mode=config["mode"],
        rep_splits=tuple(config["rep_splits"]),
        base_routers=tuple(config["base_routers"]),
        heuristics=tuple(config["heuristics"]),
        init_mapping=init_mapping,
        inc_cc=bool(config["inc_cc"]),
        inc_naive_h=bool(config["inc_naive_h"]),
        exp_qasm=bool(config["exp_qasm"]),
        exp_csv=bool(config["exp_csv"]),
        verbose=verbose
    )
