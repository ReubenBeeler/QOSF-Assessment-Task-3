from QuantumComputer import *

# global params are "default"
# example:
# my_circuit = [
#   {"gate": "u3",
#    "params": {"theta": "default", "phi": 3.1415, "lambda": "default"}}
# ]

# global params are defined here:
defaults = {"theta": 3.1415, "phi": 1.5708, "lambda": -3.1415}

my_circuit = [
  {"gate": "u3", "params": {"theta": 1, "phi": 2, "lambda": -1}, "target": [0]},
  {"gate": "cx", "target": [0, 1]}
]

# Initial state of 4 qubits: [1, 0, 0, ..., 0]
initial_state = get_ground_state(4)

# State after applying circuit to initial state
final_state = run_program(initial_state, my_circuit, defaults)

counts = get_counts(final_state, 1000)

# results
print(counts)
