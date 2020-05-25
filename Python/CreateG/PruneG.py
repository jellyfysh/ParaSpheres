import math
import copy
import GenerateG3
from GenerateG3 import configuration

#
# This program prunes the constraint graph G3
#


def b_contact(sphere_1, sphere_2, dirc, radius):
    d_perp = abs(sphere_2[not dirc] - sphere_1[not dirc])
    d_perp = min(d_perp, 1.0 - d_perp)
    if d_perp > 2.0 * radius:
        return -float("inf")
    else:
        d_para = math.sqrt(4.0 * radius ** 2 - d_perp ** 2)
        return d_para


direction = 0
sigma = GenerateG3.sigma
# Order of pruning
# Will appear in the name of the output file after 'G_pruned'
order = 4
# Whether doing symmetrization
# If true, 's' will appear in the name of the output file before '.dat'
symmetrize = True
forward_nodes = []
backward_nodes = []
with open("G3_1.dat", "r") as input_config:
    for line in input_config:
        line = [int(a) for a in line.split()]
        forward_nodes.append(set(line))
print('Read in forward finished')
if symmetrize:
    with open("G3_-1.dat", "r") as input_config:
        for line in input_config:
            line = [int(a) for a in line.split()]
            backward_nodes.append(set(line))
    print('Read in backward finished')

contact_matrix = {}
for a in range(GenerateG3.N):
    for b in forward_nodes[a]:
        contact_matrix[(a, b)] = b_contact(configuration[a], configuration[b], 0, sigma)
back_contact_matrix = {}
if symmetrize:
    for a in range(GenerateG3.N):
        for b in backward_nodes[a]:
            back_contact_matrix[(a, b)] = b_contact(configuration[a], configuration[b], 0, sigma)


def find_path(graph, in_paths, remaining_order):

    out_paths = []
    for path in in_paths:
        for node in graph[path[-1]]:
            if node not in path:
                path_temp = path[:]
                path_temp.append(node)
                out_paths.append(path_temp)
    remaining_order -= 1
    return out_paths, remaining_order


def directed_distance(velocity, particle1, particle2):
    # assume particle2 is in front of particle1 when v=1
    if velocity == 1:
        return (particle2[0] - particle1[0] + 1.0) % 1
    elif velocity == -1:
        return (particle1[0] - particle2[0] + 1.0) % 1
    else:
        print('Bad velocity')
        return None


def shielded(graph, contact, velocity, index, start, end, remaining_order):
    paths = [[start]]
    while remaining_order > 0 and len(paths) > 0:
        new_paths = []
        for path in paths:
            if directed_distance(velocity, configuration[index], configuration[end]) > \
                    directed_distance(velocity, configuration[index], configuration[path[-1]]):
                new_paths.append(path)
            else:
                continue
            b_end = b_contact(configuration[end], configuration[path[-1]], direction, sigma)
            b_total = 0
            dx = 0
            if len(path) == 1:
                if contact[(index, start)] + b_end > contact[(index, end)]:
                    return True
            else:
                for p in range(len(path) - 1):
                    b_total += contact[(path[p], path[p + 1])]
                    dx += directed_distance(velocity, configuration[path[p]], configuration[path[p + 1]])
                dx += directed_distance(velocity, configuration[path[-1]], configuration[end])
                dx += directed_distance(velocity, configuration[index], configuration[start])
                if contact[(index, start)] + b_total + b_end > contact[(index, end)] and dx < 1:
                    return True
        paths, remaining_order = find_path(graph, new_paths, remaining_order)


def pruning(graph, contact, velocity, index, remaining_order):
    pruned = set()
    neighbour_list = list(graph[index])
    if len(graph[index]) == 3:
        list_dist = [directed_distance(velocity, configuration[index], configuration[neighbour_list[0]]),
                     directed_distance(velocity, configuration[index], configuration[neighbour_list[1]]),
                     directed_distance(velocity, configuration[index], configuration[neighbour_list[2]])]
        close_index = list_dist.index(min(list_dist))
        close = neighbour_list[close_index]
        far_index = list_dist.index(max(list_dist))
        far = neighbour_list[far_index]
        middle = list(neighbour_list)
        middle.remove(close)
        middle.remove(far)
        middle = middle[0]
        if shielded(graph, contact, velocity, index, close, middle, remaining_order):
            pruned.add(middle)
        if shielded(graph, contact, velocity, index, close, far, remaining_order):
            pruned.add(far)
        elif shielded(graph, contact, velocity, index, middle, far, remaining_order):
            pruned.add(far)
    elif len(graph[index]) == 2:
        list_dist = [directed_distance(velocity, configuration[index], configuration[neighbour_list[0]]),
                     directed_distance(velocity, configuration[index], configuration[neighbour_list[1]])]
        close_index = list_dist.index(min(list_dist))
        close = neighbour_list[close_index]
        far_index = list_dist.index(max(list_dist))
        far = neighbour_list[far_index]
        if shielded(graph, contact, velocity, index, close, far, remaining_order):
            pruned.add(far)
    return pruned


forward_nodes_copy = copy.deepcopy(forward_nodes)
backward_nodes_copy = copy.deepcopy(backward_nodes)
#
# Forward pruning
# Backward pruning if symmetrize
#
for i in range(GenerateG3.N):
    forward_nodes[i].difference_update(pruning(forward_nodes_copy, contact_matrix, 1, i, order))
if symmetrize:
    for i in range(GenerateG3.N):
        backward_nodes[i].difference_update(pruning(backward_nodes_copy, back_contact_matrix, -1, i, order))
#
# Symmetrization of forward and backward constraint graph
#
if symmetrize:
    for i in range(GenerateG3.N):
        temp_forward_node = set([])
        for j in forward_nodes[i]:
            if i in backward_nodes[j]:
                temp_forward_node.add(j)
        forward_nodes[i] = copy.copy(temp_forward_node)

if symmetrize:
    output = open('G_pruned' + str(order) + 's.dat', 'w')
else:
    output = open('G_pruned' + str(order) + '.dat', 'w')
number_edges = 0
for nodes in forward_nodes:
    string = ' '
    for elements in nodes:
        number_edges += 1
        string += str(elements) + ' '
    string += '\n'
    output.write(string)
output.close()
print('{} edges in forward G'.format(number_edges))

if symmetrize:
    for i in range(GenerateG3.N):
        for j in forward_nodes[i]:
            if i not in backward_nodes[j]:
                print('Not symmetric')
