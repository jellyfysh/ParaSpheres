# JeLLyFysh/ParaSpheres - Multithreaded event-chain Monte Carlo with local times
# - https://github.com/jellyfysh/paraspheres
# Copyright (C) 2020 The JeLLyFysh organization
# (see the AUTHORS file on jellyfysh/paraspheres for the full list of authors)
#
# This file is part of JeLLyFysh/ParaSpheres.
#
# JeLLyFysh/ParaSpheres is free software: you can redistribute it and/or modify
# it under the terms of the GNU General
# Public License as published by the Free Software Foundation, either > version
# 3 of the License, or (at your option)
# any later version.
#
# JeLLyFysh/ParaSpheres is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even 
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# JeLLyFysh/ParaSpheres in the LICENSE file. If not, see
# <https://www.gnu.org/licenses/>.
#
# If you use JeLLyFysh/ParaSpheres in published work, please cite the following
# reference
# (see [Li2020] in References.bib):
# Botao Li, Synge Todo, A. C. Maggs, Werner Krauth
# Multithreaded event-chain Monte Carlo with local times,
# arXiv e-prints: 2004.11040 (2020), https://arxiv.org/abs/2004.11040
#

"""
This program inputs a starting configuration and a precomputed constraint graph 
and then performs a global time ECMC computation: After each
event, all active positions are time-sliced, so that conflicts are avoided.
This program supports k active spheres. k will be given as a external parameter.
The indices of active spheres are read from file.
If the file of active spheres is absent, the indices of the active spheres are generated.
This program outputs the set of liftings without the positions.
The output is written into a file 'lifting(number of active particles).dat'
(The total length is the sum of the chain length for each active sphere.)
The format of output: time active_sphere target_sphere
"""

import math
import time
import sys
import os.path


def b_contact(sphere_1, sphere_2, dirc, radius):
    d_perp = abs(sphere_2[not dirc] - sphere_1[not dirc])
    d_perp = min(d_perp, 1.0 - d_perp)
    if d_perp > 2.0 * radius:
        return -float("inf")
    else:
        d_para = math.sqrt(4.0 * radius ** 2 - d_perp ** 2)
        return d_para


k_active = int(sys.argv[1])
N = 1
sigma = 0
validation_chain_length = 0
k_active_max = 0
with open('info.dat', 'r') as input_info:
    for line in input_info:
        line_split = line.split()
        if line_split[0] == 'N':
            N = int(line_split[1])
        if line_split[0] == 'sigma':
            sigma = float(line_split[1])
        if line_split[0] == 'ChainLength':
            validation_chain_length = float(line_split[1])
        if line_split[0] == 'LogKmax':
            k_active_max = 2**int(line_split[1])
print(N, sigma, validation_chain_length, k_active_max)
if k_active > k_active_max:
    sys.exit('Too many active spheres!')
if k_active >= N:
    sys.exit('Inconsistent number of active spheres!')

with open("configuration.dat", "r") as input_config:
    configuration = []
    for line in input_config:
        line = [float(a) for a in line.split()]
        configuration.append(line)
input_config.close()
with open('constraints.dat', 'r') as input_constraint:
    forward_nodes = []
    for line in input_constraint:
        line = [int(a) for a in line.split()]
        forward_nodes.append(line)
input_constraint.close()
active_spheres = []
active_file = 'active' + str(k_active) + '.dat'
if os.path.exists(active_file):
    print('Existing active list')
    with open(active_file, 'r') as input_active:
        for line in input_active:
            active_spheres.append(int(line))
else:
    print('Creating active list')
    for line in range(k_active):
        active_spheres.append(int(line * N / k_active))
    output_actives = open(active_file, 'w')
    for line in active_spheres:
        output_actives.write(str(line)+'\n')
    output_actives.close()
print('Number of active spheres: ', k_active)
# print(active_spheres)
output_file = 'lifting' + str(k_active) + '.dat'
if os.path.exists(output_file):
    sys.exit("output file already exists, cannot be overwritten ")
lifting_output = open(output_file, 'w')
print(validation_chain_length, ' total chain length')
direction = 0
distance = validation_chain_length / k_active
#
#    contact-matrix computation
#
b_matrix = {}
for a in range(N):
    for b in forward_nodes[a]:
        b_matrix[(a, b)] = b_contact(configuration[a], configuration[b], direction, sigma)

number = 0
t_0 = time.time()

global_time = 0.0
old_a = -1
new_a = -1
while distance > 0.0:
    number += 1
    event_min = distance
    for a in active_spheres:
        for b in forward_nodes[a]:
            if b not in active_spheres:
                event_b = configuration[b][direction] - configuration[a][direction]
                if event_b < 0:
                    event_b += 1.0
                event_b -= b_matrix[(a, b)]
                event_b = max(event_b, 0.00000000001)
                if event_b < event_min:
                    old_a = a
                    new_a = b
                    event_min = event_b
    for a in active_spheres:
        configuration[a][direction] = configuration[a][direction] + event_min
        if configuration[a][direction] > 0.5:
            configuration[a][direction] -= 1.0
    global_time += event_min
    if event_min < distance:
        active_spheres.remove(old_a)
        active_spheres.append(new_a)
        lifting_output.write("{:12.10E} {} {}\n".format(global_time, old_a, new_a))
    distance -= event_min
t_end = time.time()
lifting_output.close()
active_spheres.sort()
number -= 1  # End of the chain is not counted as event
print("Benchmark takes {} seconds".format(t_end - t_0))
print("{} events".format(number))
print("{:.2E} estimated event per hour".format(number / (t_end - t_0) * 3600))
