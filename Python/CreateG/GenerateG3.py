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

import math
import copy
import sys

#
# This program computes the constraint graph G3
#

def dist(a, b):
    d_x = abs(a[0] - b[0])
    d_x = min(d_x, 1.0 - d_x)
    d_y = abs(a[1] - b[1])
    d_y = min(d_y, 1.0 - d_y)
    return math.sqrt(d_x**2 + d_y**2)


def dist_perp_para(a, b, dirc):
    perp = (b[not dirc] - a[not dirc] + 1.0) % 1.0
    if perp > 0.5:
        perp -= 1.0
    para = (b[dirc] - a[dirc] + 1.0) % 1.0
    return perp, para


N = 1
sigma = 0
with open('info.dat', 'r') as input_info:
    for line in input_info:
        line_split = line.split()
        if line_split[0] == 'N':
            N = int(line_split[1])
        elif line_split[0] == 'sigma':
            sigma = float(line_split[1])
print('Number of spheres {}'.format(N), 'Sphere radius {}'.format(sigma))
with open("configuration.dat", "r") as input_config:
    configuration = []
    for line in input_config:
        line = [float(a) for a in line.split()]
        configuration.append(line)
print(' read in finished')
direction = 0  # assume active spheres moving in x direction


def main():
    forward_nodes = []
    backward_nodes = []
    #
    #   treat the three lanes separately (before pruning).
    #
    for i in range(N):
        if i % 1000 == 0:
            print(i, '  spheres treated')
        lane_lower = []
        lane_upper = []
        lane_center = []
        for j in range(N):
            if j != i:
                dist_perp, dist_para = dist_perp_para(configuration[i], configuration[j], direction)
                if -2.0 * sigma < dist_perp < - sigma:
                    lane_lower.append((dist_para, j))
                elif - sigma < dist_perp < sigma:
                    lane_center.append((dist_para, j))
                elif sigma < dist_perp < 2.0 * sigma:
                    lane_upper.append((dist_para, j))
        lane_lower.sort()
        lane_upper.sort()
        lane_center.sort()
        forward_nodes.append(set([]))
        backward_nodes.append(set([]))
        if lane_lower:
            forward_nodes[i].add(lane_lower[0][1])
            backward_nodes[i].add(lane_lower[-1][1])
        if lane_center:
            forward_nodes[i].add(lane_center[0][1])
            backward_nodes[i].add(lane_center[-1][1])
        if lane_upper:
            forward_nodes[i].add(lane_upper[0][1])
            backward_nodes[i].add(lane_upper[-1][1])

    output = open('G3_1.dat', 'w')
    number_edges = 0
    for nodes in forward_nodes:
        string = ' '
        for elements in nodes:
            number_edges += 1
            string += str(elements) + ' '
        string += '\n'
        output.write(string)
    output.close()
    print('{} edges in forward G^(3)'.format(number_edges))
    output = open('G3_-1.dat', 'w')
    number_edges = 0
    for nodes in backward_nodes:
        string = ' '
        for elements in nodes:
            number_edges += 1
            string += str(elements) + ' '
        string += '\n'
        output.write(string)
    output.close()
    print('{} edges in backward G^(3)'.format(number_edges))


if __name__ == '__main__':
    main()
