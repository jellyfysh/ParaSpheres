"""
This program inputs a starting configuration and a precomputed constraint graph 
and then runs single-threaded ECMC with local time (Algorithm 4 in [Li2020]).
This program supports k active particles, with k = 1,2,4,... 8192.
The initial indices of active particles are read from file.
"""

import math
import time
import sys
import os.path
import copy
import random


def b_contact(sphere_1, sphere_2, dirc, radius):
    d_perp = abs(sphere_2[not dirc] - sphere_1[not dirc])
    d_perp = min(d_perp, 1.0 - d_perp)
    if d_perp > 2.0 * radius:
        return -float("inf")
    else:
        d_para = math.sqrt(4.0 * radius ** 2 - d_perp ** 2)
        return d_para


def advance(dirc, config, in_active_spheres, liftings, breakpoint_advance):
    """ 
    Fast-forward a given configuration config, by a time interval breakpoint_advance,
    given a reference set of liftings. Presently, t=0 is implemented.
    """
    new_active_spheres = in_active_spheres[:]
    new_config = copy.deepcopy(config)
    n_liftings = 0
    old_time = 0.0
    while True:
        t_lifting = liftings[n_liftings][0]
        time_of_flight = min(t_lifting - old_time, breakpoint_advance - old_time)
        for a in new_active_spheres:
            new_config[a][dirc] += time_of_flight
        if t_lifting < breakpoint_advance:
            new_active_spheres.remove(liftings[n_liftings][1])
            new_active_spheres.append(liftings[n_liftings][2])
            n_liftings += 1
            old_time = t_lifting
        else:
            break
    return n_liftings, new_config, new_active_spheres


N = 1
sigma = 0
validation_chain_length = 0
k_active_max = 0
k_active = int(sys.argv[1])
with open('../Data.run/info.dat', 'r') as input_info:
    for line in input_info:
        line_split = line.split()
        if line_split[0] == 'N':
            N = int(line_split[1])
        if line_split[0] == 'sigma':
            sigma = float(line_split[1])
        if line_split[0] == 'ChainLength':
            validation_chain_length = float(line_split[1])
        if line_split[0] == 'LogKmax':
            k_active_max = 2 ** int(line_split[1])
epsilon = 0.000000001  # this is the value to decide whether two reals are the same.

with open("../Data.run/configuration.dat", "r") as input_config:
    configuration = []
    local_times = []
    for line in input_config:
        line = [float(i) for i in line.split()]
        configuration.append(line)
        local_times.append(0)  # starting_time
with open('../Data.run/constraints.dat', 'r') as input_constraint:
    forward_nodes = []
    for line in input_constraint:
        line = [int(i) for i in line.split()]
        forward_nodes.append(line)
initial_active_spheres = []
with open('../Data.run/active' + str(k_active) + '.dat', 'r') as input_actives:
    for line in input_actives:
        initial_active_spheres.append(int(line))
if k_active != len(initial_active_spheres):
    sys.exit("problem in k_active")

reference_liftings = []
with open('../Data.run/lifting' + str(k_active) + '.dat', 'r') as benchmark_trajectory_file:
    for line in benchmark_trajectory_file:
        line = line.split()
        line = (float(line[0]), int(line[1]), int(line[2]))
        reference_liftings.append(line)

direction = 0
new_breakpoint = validation_chain_length / k_active
#
#    contact-matrix computation
#
b_matrix = {}
for i in range(N):
    for j_tilde in forward_nodes[i]:
        b_matrix[(i, j_tilde)] = b_contact(configuration[i], configuration[j_tilde], direction, sigma)

t = 0
breakpoint = random.uniform(0, new_breakpoint * 0.9)
n_liftings_offset, configuration, initial_active_spheres = \
    advance(direction, configuration, initial_active_spheres, reference_liftings, breakpoint)
print('distance_to_go, Delta_t', new_breakpoint, breakpoint)

run_output = open('output.dat', 'w')

shifted_new_breakpoint = new_breakpoint - breakpoint
active_spheres = list(initial_active_spheres[:])  # not yet at new_horizon: event_horizon > 0
stalled_spheres = set([])  # at new_horizon: event_horizon = 0
all_liftings = []
earliest_horizon_violation = float("inf")
lifting_divergence = float("inf")
while active_spheres:  # != set():
    i = random.choice(active_spheres)
    active_spheres.remove(i)
    j = i  # this allows the program to run for N = 1 disks.
    distance = shifted_new_breakpoint - local_times[i]  # this, by definition of active thread, is > 0.
    tau = distance  # this is the time-of-flight to the next breakpoint
    for j_tilde in forward_nodes[i]:
        tau_ij_tilde = configuration[j_tilde][direction] - configuration[i][direction]
        if tau_ij_tilde < 0.0:
            tau_ij_tilde += 1.0
        tau_ij_tilde = tau_ij_tilde - b_matrix[(i, j_tilde)]  # can be slightly negative!
        if local_times[j_tilde] > local_times[i] + tau_ij_tilde:
            earliest_horizon_violation = min(earliest_horizon_violation, local_times[i], local_times[j_tilde])
        if tau_ij_tilde < tau:
            j = j_tilde
            tau = tau_ij_tilde
    configuration[i][direction] = configuration[i][direction] + tau
    if configuration[i][direction] > 0.5:
        configuration[i][direction] -= 1.0
    lifting_time = local_times[i] + tau
    local_times[i] = lifting_time
    if j not in active_spheres:
        if tau < distance:
            all_liftings.append((lifting_time, i, j))
            i = j
            local_times[i] = lifting_time
            active_spheres.append(i)
        else:
            stalled_spheres.add(i)
    else:
        if tau < distance:
            active_spheres.append(i)
        else:
            stalled_spheres.add(i)
all_liftings.sort()
for k in range(len(all_liftings)):
    koff = k + n_liftings_offset
    run_output.write(
        "{:12.10E} {} {}\n".format(all_liftings[k][0] + breakpoint, all_liftings[k][1], all_liftings[k][2]))
    if (abs(all_liftings[k][0] - reference_liftings[koff][0] + breakpoint) > epsilon or
            all_liftings[k][1] != reference_liftings[koff][1] or
            all_liftings[k][2] != reference_liftings[koff][2]):
        lifting_divergence = min(all_liftings[k][0], reference_liftings[koff][0])
        break
if earliest_horizon_violation <= lifting_divergence:
    print("OK")
else:
    print("Not OK")
