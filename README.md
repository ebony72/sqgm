# SQGM: Optimising Quantum Circuit Depth in Qubit Mapping

This repository provides the source code and reproducible experiments presented in Li et al.'s *Single-Qubit Gates Matter for Optimising Quantum Circuit Depth in Qubit Mapping*.

## Motivation

Most existing quantum circuit transformation (QCT, or qubit mapping) algorithms prioritise minimising circuit size, potentially overlooking the impact of single-qubit gates on circuit depth. In the past several years, we have seen the size of quantum computer increases [from 5 to 433 qubits](https://newsroom.ibm.com/2022-11-09-IBM-Unveils-400-Qubit-Plus-Quantum-Processor-and-Next-Generation-IBM-Quantum-System-Two), but qubit coherence time in NISQ devices remains very short. This implies that we can only run a very limited number of quantum operations on each qubit; in other words, we cannot extract meaningful information from very deep quantum circuits on NISQ devices. Thus, minimising the depth of transformed circuits is perhaps a more important objective.

SQGM (Single-Qubit Gates Matter) promises a simple and effective method that takes into account the impact of single-qubit gates on circuit depth. Our method can be integrated with existing QCT algorithms. In this source code, we demonstrate the effectiveness of SQGM by embedding it in Qiskit's  [`SabreSwap` module](https://qiskit.org/documentation/stubs/qiskit.transpiler.passes.SabreSwap.html), which is an implementation of the state-of-the-art [SABRE (SWAP-based Bidirectional Heuristic Search Algorithm)](https://dl.acm.org/doi/abs/10.1145/3297858.3304023).

## Installation

SQGM uses modules from the Qiskit libraries, in particular, [Qiskit Terra](https://github.com/Qiskit/qiskit-terra). Information on installing Qiskit can be found [here](https://qiskit.org/documentation/getting_started.html). The version used for experiments was 0.39.4, the latest version at the time. In addition, the [NetworkX](https://www.osti.gov/biblio/960616) library is used to construct architecture graphs. Information on installing networkx can be found [here](https://networkx.org/documentation/stable/install.html).

## Architecture graphs

All experiments conducted used architecture graphs of the following quantum devices as target backends:

- `ourense`: IBM Q Ourense (5 qubits),
- `tokyo`: IBM Q Tokyo (20 qubits),
- `rochester`: IBM Q Rochester (53 qubits),
- `sycamore53`, `sycamore54`: Google Sycamore (53 and 54 qubits).

Note that Google Sycamore was originally designed with 54 qubits, but ended up with 1 bad qubit; both the designed and actual versions were used in the experiments. Architecture graphs were modelled using NetworkX, and can be found in `ag.py`. A visualisation of these graphs can be found in `display.ipynb`.

## Mappers and routers

All experiments conducted used Qiskit 0.39.4's `SabreLayout` module to generate initial mappings. Nonetheless, `/mapper` provides the Qiskit 0.33.0's `SabreLayout` module for reference and use, if needed.

The source code for `SQGMSwap` (written on top of `SabreSwap`) as well as Qiskit 0.33.0's `SabreSwap` module are provided in `/router`. The experiments also used Qiskit 0.39.4's `SabreSwap` module and `NASSCSwap` module, the original implementation of the [NASSC (Not All SWAPs Have the Same Cost)](https://ieeexplore.ieee.org/abstract/document/9773196/?casa_token=XzV4yy5W3D8AAAAA:ioe4xkNhEWtNZyiW0eWFsBf7WGRfpfAY7fBC5hwCRA4nzsTH2OvUG6OChXmQdbo_sU_aNiuc), retrieved from the author's [repository](https://github.com/peiyi1/nassc_code) solely for experimental purposes. [TOQM (Time-Optimal Qubit Mapping)](https://dl.acm.org/doi/10.1145/3445814.3446706) was also used for comparisons; its source code (written in C++) is available on the author's [repository](https://github.com/time-optimal-qmapper/TOQM).

## Benchmarks

All benchmarks used in experiments are provided in `/benchmark`. Contained in the directory are selected circuits from [QUEKO](https://ieeexplore.ieee.org/abstract/document/9140293?casa_token=sEb0o07UvpAAAAAA:TC3mMhNOUM6m0MmWxhikvqPxrsfcIKpjMjpLwsutmrf7UlcWBnfHhqNvWfg7JYpI7UTE1CoI), [QUEKNO](https://arxiv.org/abs/2301.08932), and [MQT Bench](https://www.cda.cit.tum.de/mqtbench/). In addition, `/example_4q` contains the original circuit, and initial mapping, of the running example in the paper (refer to Fig. 1).

## Experiments

Original experiment results are archived in `/sqgm_data`. Within the directory:

- `/example_4q-out` provides the output transformed circuits of the original `example_4q` circuit found in `/benchmark`.
- `/mqtbench-init` provides the initial mappings of each of the MQT Bench circuits found in `/benchmark`.
- `/exps` provides the results of the original experiments.
- `sqgm-data.xlsx` provides a summary of all results presented in the paper.

Details of the original experiments, in JSON format, can be found in `/exp-in`. This includes the validity test of the running example `example_4q`, whose result is provided in `display.ipynb`.

### Experiment configurations

All of the original experiments can be reproduced via their respective JSON configuration files in `/exp-in`. Each file contains the following information:

| Field         | Description 
|:-             |:-
|`bench_path`   | File path to the benchmarks to be used as input circuits.
|`bench_filter` | Filter for specific circuit within the benchmark specified by `bench_path`. Only circuits containing `bench_filter` in their names are run.
|`qasm_path`    | File path to the directory where the exported output circuit QASM files are exported.
|`csv_path`     | File path to the CSV file in which numerical results are exported.
|`ag`           | Architecture graph of the target backend.
|`reps`         | Number of reps.
|`objective`    | Indicate whether `depth` (transformed circuit depth) or `size` (number of 2-qubit gates in the transformed circuit) should be optimised.
|`mode`         | Indicate whether to record the `best` result, or `all` results amongst the specified number of reps.
|`rep_splits`   | Rep numbers to be reported in the exported CSV file. Each split corresponds to a column in the CSV, which reports the best result *up to* that split. This field is only relevant for `all` mode. 
|`base_routers` | Routers (amongst those provided by `/router`) to be run.
|`heuristics`   | Heuristics to be used.
|`inc_cc`       | Indicate whether or not to include, on top of each base router, an augmented router with Qiskit's [commutative gate cancellation pass](https://qiskit.org/documentation/stubs/qiskit.transpiler.passes.CommutativeCancellation.html) during post-routing optimisations. This field is only relevant when some router other than NASSC is included in `base_routers` (since `NASSCSwap` includes this pass by default).
|`inc_naive_h`  | Indicate whether or not to include, on top of the original SQGM, a variation where SABRE's original $H_\textsf{decay}$ is used instead of SQGM's enhanced $H_\textsf{sqgm}$ when selecting SWAPs. This field is only relevant when SQGM is included in `base_routers`.
|`exp_qasm`     | Indicate whether or not to export the output circuit QASM files. Note that if a QASM file of the same name exists, it will be overwritten.
|`exp_csv`      | Indicate whether or not to export numerical results to the CSV file. Note that if a CSV of the same name exists, it will be overwritten.

### Running an experiment

To run an experiment, simply run `main.py` and supply the path to the respective JSON configuration file as a command line argument. Additionally, a `--verbose`, or its shorthand `-v`, can be included to toggle the display of detailed information on each rep. For example, to run `exp1` in verbose mode:
```bash
$ python main.py ./exp-in/exp1.json --verbose
```

### Customising and defining an experiment

It is simple to customise, or define a new, experiment: simply edit the respective, or create a new, JSON configuration file as appropriate. Note that:

- All path fields, i.e., `bench_path`, `qasm_path`, and `csv_path`, must be suffixed with `/`.
- If `bench_filter` is left blank, i.e., `""`, all circuits in `bench_path` will be run.
- The value of `ag` must be one of the keys of `ARCHGRAPHS`, found in `util.py`.
- `reps` must be a non-zero positive integer.
- If `rep_splits` is left empty in `all` mode, there will be a column for *every* rep number in the exported CSV file.
- Each value in `base_routers` must be one of the keys of `ROUTERS`, found in `util.py`.

An empty JSON configuration file `exp_template.json` has been provided, and may be used to define new experiments.

### Adding an architecture graph and router

It is possible to add new architecture graphs and routers to be used for new experiments.

To add a new architecture graph, define a new function in `ag.py` which constructs and returns a `Graph` object (from NetworkX). Next, in `util.py`, add a new key-value pair for `ARCHGRAPHS`, where key is the name and value is a call to the function. Then, experiments can be conducted with the new architecture graph by specifying its key in the `ag` field of the respective JSON configuration file.

To add a new router, add the source code of that router in `/router`. Note that the router must be compatible in nature with Qiskit's `SabreSwap`, i.e., its constructor must have have at least the following parameters:

- `coupling_map` (`CouplingMap`);
- `heuristics` (`str`), which must be either `basic`, `lookahead`, or `decay`;
- `seed` (`int`); and
- `fake_run` (`bool`).

Next, in `util.py`, import the router's module and add a new key-value pair for `ROUTERS`, where key is the name and value is the relevant module. Then, experiments can be conducted with the new router by including its key in the `base_routers` field of the respective JSON configuration file.