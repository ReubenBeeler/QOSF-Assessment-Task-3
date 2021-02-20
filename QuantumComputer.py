import numpy as np
import math
from random import uniform as randFloat

def get_ground_state(num_qubits):
    if num_qubits < 1: return
    else:
        ground_state = [1, 0]
    for i in range(1, num_qubits):
        ground_state = np.kron(ground_state, [1, 0])
    return ground_state

def U3(stringsThetaPhiLambda):
    theta = float(stringsThetaPhiLambda[0])
    phi = float(stringsThetaPhiLambda[1])
    lambda_ = float(stringsThetaPhiLambda[2])
    _00 = math.cos(theta / 2)
    _01 = -math.e**(complex(0, lambda_)) * math.sin(theta / 2)
    _10 = math.e**(complex(0, phi)) * math.sin(theta / 2)
    _11 = math.e**(complex(0, lambda_ + phi)) * math.cos(theta / 2)
    return np.array([
        [_00, _01],
        [_10, _11]
        ])

def get_operator(total_qubits, gate_unitary, target_qubits):
    # If the gate is a U3 gate, gate_unitary will be the string
    # "u3([theta], [phi], [lambda])", e.g. "u3(3.14, 0, -1.57)"
    gate = gate_unitary.lower()
    operators = []
    I = np.array([
        [1, 0],
        [0, 1]
        ])
    X = np.array([
        [0, 1],
        [1, 0]
        ])
    if gate in ["x", "h"] or gate[:2] == "u3":
        if gate == "x":
            U = X
        elif gate == "h":
            ir2 = 2**(-0.5)
            U = np.array([
                [ir2, ir2],
                [ir2,-ir2]
                ])
        else:
            U = U3(gate[3:len(gate)-1].split(", "))
        for i in range(0, total_qubits):
            if i in target_qubits: 
                operators.append(U)
            else: 
                operators.append(I)
        operator = operators[0]
        for i in range(1, len(operators)):
            operator = np.kron(operator, operators[i])
        return operator
    elif gate == "cx" or gate == "cnot":
        early_ops = [[], []]
        for i in range(0, total_qubits):
            if i == target_qubits[0]:
                early_ops[0].append(np.array([
                    [1, 0],
                    [0, 0]
                    ]))
            else:
                early_ops[0].append(I)
        for i in range(0, total_qubits):
            if i == target_qubits[0]:
                early_ops[1].append(np.array([
                    [0, 0],
                    [0, 1]
                    ]))
            elif i == target_qubits[1]:
                early_ops[1].append(X)
            else:
                early_ops[1].append(I)
        early_op = [early_ops[0][0], early_ops[1][0]]
        for i in range(1, total_qubits):
            early_op[0] = np.kron(early_op[0], early_ops[0][i])
            early_op[1] = np.kron(early_op[1], early_ops[1][i])
        return early_op[0] + early_op[1]
    else:
        return

def run_program(initial_state, program, defaults):
    state = initial_state
    q = int(math.log(len(state), 2))
    for event in program:
        if event["gate"] == "u3":
            rawTheta = event["params"]["theta"]
            rawPhi = event["params"]["phi"]
            rawLambda = event["params"]["lambda"]
            # Ternary operator to set numerical value to theta, phi, or lambda
            # if the value is the string "default", according to the dictionary
            # parameter "defaults"
            theta = defaults["theta"] if rawTheta == "default" else rawTheta
            phi = defaults["phi"] if rawPhi == "default" else rawPhi
            lambda_ = defaults["lambda"] if rawLambda == "default" else rawLambda
            state = np.dot(get_operator(
                q,
                "u3(" + str(theta) + ", " + str(phi) + ", " + str(lambda_) + ")",
                event["target"]
                ), state)
        else:
            state = np.dot(get_operator(q, event["gate"], event["target"]), state)
    return state

def measure_all(state_vector):
    # I used my own method of finding a random index from a list of weights.
    levels = []
    # each "level" in levels represents the probabilities of the state vector, but stacked.
    # Take the following quantum state: [sqrt(0.5), 0, sqrt(0.25), sqrt(0.25)].
    # levels would be [0.5, 0.5 + 0, 0.5 + 0 + 0.25, 0.5 + 0 + 0.25 + 0.25], or [0.5, 0.5, 0.75, 1].
    # I then generate a random float from 0 to the total weight (1). The first level (from left
    # to right) that is greater than the random float is the chosen index. If the random float is
    # 0.6, then the index would be 0b10 (2).
    levels.append(abs(state_vector[0])**2)
    for i in range(1, len(state_vector)):
        levels.append(levels[i-1] + (abs(state_vector[i])**2))
    val = randFloat(0, 1)
    for i in range(0, len(levels)):
        if (val < levels[i]):
            return i

def dictionaryToString(dictionary):
    string = "{"
    for tag in dictionary:
        string += "\n\t\"" + tag + "\": " + str(dictionary[tag]) + ","
    return (string[:len(string)-1] + "\n}")

def get_counts(state_vector, shots):
    counts = {}
    q = int(math.log(len(state_vector), 2))
    for i in range(0, len(state_vector)):
        bin = "{0:b}".format(i)
        counts[(q - len(bin)) * "0" + bin] = 0 # appends 0s at the beginning
        # of the string to show values for each qubit. Otherwise, index 3 in
        # a three-qubit quantum state would be "11" instead of "011". "counts"
        # is a dictionary mapping the strings of the indices to their counts.
    for i in range(0, shots):
        bin = "{0:b}".format(measure_all(state_vector))
        counts[(q - len(bin)) * "0" + bin] += 1 # adds one to the count of
        # the index that is returned by the weighted randomizer in measure_all(vector)
    return dictionaryToString(counts)
