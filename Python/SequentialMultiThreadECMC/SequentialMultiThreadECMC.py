# JeLLyFysh/ParaSpheres - Multithreaded event-chain Monte Carlo with local times
# - https://github.com/jellyfysh/paraspheres
# Copyright (C) 2020 The JeLLyFysh organization
# (see the AUTHORS file on jellyfysh/paraspheres for the full list of authors)
#
# This file is part of JeLLyFysh/ParaSpheres.
#
# JeLLyFysh/ParaSpheres is free software: you can redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation, either > version 3 of the License, or (at your option)
# any later version.
#
# JeLLyFysh/ParaSpheres is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# JeLLyFysh/ParaSpheres in the LICENSE file. If not, see <https://www.gnu.org/licenses/>.
#
# If you use JeLLyFysh/ParaSpheres in published work, please cite the following reference
# (see [Li2020] in References.bib):
# Botao Li, Synge Todo, A. C. Maggs, Werner Krauth
# Multithreaded event-chain Monte Carlo with local times,
# arXiv e-prints: 2004.11040 (2020), https://arxiv.org/abs/2004.11040
#
#
# Multithreaded ECMC (Sequential-consistency model).
#
import copy
import math
import random
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

epsilon = 0.000000001


def equal(a, b):
    if a is float and b is float:
        return abs(a - b) <= epsilon
    else:
        return a == b


class State:
    def __init__(self, folder, number_actives, debug=False):
        self.debug = debug
        self.spheres = []
        self.threads = []
        self.abort = False

        file_config = folder + "/configuration.dat"
        file_info = folder + "/info.dat"
        file_active = folder + "/active" + str(number_actives) + ".dat"

        self.active_list = []
        with open(file_active, 'r') as active_list:
            for line in active_list:
                self.active_list.append(int(line))

        with open(file_info, 'r') as info:
            for line in info:
                temp = line.split()
                if temp[0] == 'sigma':
                    self.sigma = float(temp[1])
                if temp[0] == 'ChainLength':
                    self.breakpoint = float(temp[1]) / len(self.active_list)

        with open(file_config, 'r') as initial_positions:
            for line in initial_positions:
                temp = [float(p) for p in line.split()]
                temp_dict = {'tag': 'static', 'x-value': temp[0], 'y-value': temp[1], 't-value': 0}
                self.spheres.append(copy.deepcopy(temp_dict))
        print('Next breakpoint: {}'.format(self.breakpoint))
        print('Sphere radius: {}'.format(self.sigma))

        self.number_spheres = len(self.spheres)
        self.number_threads = len(self.active_list)

        self.b = [[0 for _ in range(self.number_spheres)] for _ in range(self.number_spheres)]
        for i in range(self.number_spheres):
            self.b[i][i] = -float('Inf')
            for j in range(i):
                if self.spheres[i]['x-value'] < self.spheres[j]['x-value']:
                    self.b[i][j] = self.contact(self.spheres[i], self.spheres[j])
                    self.b[j][i] = -float('inf')
                else:
                    self.b[j][i] = self.contact(self.spheres[i], self.spheres[j])
                    self.b[i][j] = -float('inf')
        print(r'b_{ij}: ')
        print(self.b)

        for active in enumerate(self.active_list):
            self.threads.append({'distance': self.breakpoint, 'i': active[1], 'j': active[1], 'x': None, 'tau': None,
                                 'j_tilde': 0, 'x_tilde': None, 'tau_tilde': None, 'buffer': 1})
            self.spheres[active[1]]['tag'] = active[0]

        print('Initial State: ')
        self.output()

    def contact(self, sphere1, sphere2):
        distance_y = abs(sphere1['y-value'] - sphere2['y-value'])
        b_square = 4 * self.sigma ** 2 - distance_y ** 2
        if b_square > 0:
            return math.sqrt(b_square)
        else:
            return -float('Inf')

    # Implementation of algorithm 5
    def next_step(self, switch):
        if self.debug:
            for sphere in self.spheres:
                print(sphere)
            for iota in self.threads:
                print(iota)
        terminate = True
        for thread in self.threads:
            if not thread['buffer'] == 26:
                terminate = False
        if terminate or self.abort:
            if self.debug:
                sys.exit()
            return None
        current_thread = self.threads[switch]
        if current_thread['buffer'] == 1:
            current_thread['tau'] = current_thread['distance']
            current_thread['x'] = None
            current_thread['j'] = current_thread['i']
            current_thread['buffer'] = 2
        elif current_thread['buffer'] == 2:
            if current_thread['j_tilde'] >= self.number_spheres:
                current_thread['j_tilde'] = 0
                current_thread['buffer'] = 10
            elif current_thread['j_tilde'] == current_thread['i'] \
                    or self.b[current_thread['i']][current_thread['j_tilde']] < 0:
                current_thread['j_tilde'] += 1
                current_thread['buffer'] = 2
            else:
                current_thread['buffer'] = 3
        elif current_thread['buffer'] == 3:
            current_thread['x_tilde'] = self.spheres[current_thread['j_tilde']]['x-value']
            current_thread['buffer'] = 4
        elif current_thread['buffer'] == 4:
            current_thread['tau_tilde'] = current_thread['x_tilde'] - self.spheres[current_thread['i']]['x-value'] - \
                                          self.b[current_thread['i']][current_thread['j_tilde']]
            current_thread['buffer'] = 5
        elif current_thread['buffer'] == 5:
            if self.spheres[current_thread['i']]['t-value'] + current_thread['tau_tilde'] < \
                    self.spheres[current_thread['j_tilde']]['t-value'] + epsilon:
                self.abort = True
            else:
                current_thread['buffer'] = 6
        elif current_thread['buffer'] == 6:
            if current_thread['tau_tilde'] < current_thread['tau']:
                current_thread['buffer'] = 7
            else:
                current_thread['j_tilde'] += 1
                current_thread['buffer'] = 2
        elif current_thread['buffer'] == 7:
            current_thread['j'] = current_thread['j_tilde']
            current_thread['buffer'] = 8
        elif current_thread['buffer'] == 8:
            current_thread['x'] = current_thread['x_tilde']
            current_thread['buffer'] = 9
        elif current_thread['buffer'] == 9:
            current_thread['tau'] = current_thread['tau_tilde']
            current_thread['j_tilde'] += 1
            current_thread['buffer'] = 2
        elif current_thread['buffer'] == 10:
            if self.spheres[current_thread['j']]['tag'] == 'static':
                self.spheres[current_thread['j']]['tag'] = switch
            current_thread['buffer'] = 11
        elif current_thread['buffer'] == 11:
            if self.spheres[current_thread['j']]['tag'] == switch:
                current_thread['buffer'] = 12
            else:
                current_thread['buffer'] = 1
        elif current_thread['buffer'] == 12:
            if not equal(current_thread['tau'], current_thread['distance']):
                current_thread['buffer'] = 13
            else:
                current_thread['buffer'] = 21
        elif current_thread['buffer'] == 13:
            if equal(current_thread['x'], self.spheres[current_thread['j']]['x-value']):
                current_thread['buffer'] = 14
            else:
                current_thread['buffer'] = 20
        elif current_thread['buffer'] == 14:
            self.spheres[current_thread['j']]['t-value'] = \
                self.spheres[current_thread['i']]['t-value'] + current_thread['tau']
            current_thread['buffer'] = 15
        elif current_thread['buffer'] == 15:
            self.spheres[current_thread['i']]['t-value'] += current_thread['tau']
            current_thread['buffer'] = 16
        elif current_thread['buffer'] == 16:
            self.spheres[current_thread['i']]['x-value'] += current_thread['tau']
            current_thread['buffer'] = 17
        elif current_thread['buffer'] == 17:
            self.spheres[current_thread['i']]['tag'] = 'static'
            current_thread['buffer'] = 18
        elif current_thread['buffer'] == 18:
            current_thread['distance'] -= current_thread['tau']
            current_thread['buffer'] = 19
        elif current_thread['buffer'] == 19:
            current_thread['i'] = current_thread['j']
            current_thread['buffer'] = 25
        elif current_thread['buffer'] == 20:
            self.spheres[current_thread['j']]['tag'] = 'static'
            current_thread['buffer'] = 1
        elif current_thread['buffer'] == 21:
            self.spheres[current_thread['i']]['t-value'] += current_thread['tau']
            current_thread['buffer'] = 22
        elif current_thread['buffer'] == 22:
            self.spheres[current_thread['i']]['x-value'] += current_thread['tau']
            current_thread['buffer'] = 23
        elif current_thread['buffer'] == 23:
            current_thread['distance'] = 0
            current_thread['buffer'] = 24
        elif current_thread['buffer'] == 24:
            self.spheres[current_thread['i']]['tag'] = 'stalled'
            current_thread['buffer'] = 25
        elif current_thread['buffer'] == 25:
            if current_thread['distance'] > 0 + epsilon:
                current_thread['buffer'] = 1
            else:
                current_thread['buffer'] = 26
        elif current_thread['buffer'] == 26:
            current_thread['buffer'] = 26
        else:
            print('Bad buffer')

    def to_string(self):
        string = 'state'
        for s in range(self.number_spheres):
            for attribute in self.spheres[0].keys():
                if self.spheres[s][attribute] is float:
                    string += '{:.4}'.format(self.spheres[s][attribute])
                else:
                    string += str(self.spheres[s][attribute])
        for s in range(self.number_threads):
            for attribute in self.threads[0].keys():
                if self.threads[s][attribute] is float:
                    string += '{:.4}'.format(self.threads[s][attribute])
                else:
                    string += str(self.threads[s][attribute])
        return string

    def output(self):
        for sphere in self.spheres:
            print(sphere)
        for iota in self.threads:
            print(iota)
        print('==========================================================================')


def main():
    initial = State('../Data.run', 2)
    thread_list = [i for i in range(len(initial.active_list))]

    pocket_states = {initial}
    state_list_complete = [initial]
    state_dictionary = {initial.to_string(): 0}
    while len(pocket_states) > 0:
        state = pocket_states.pop()
        for thread in thread_list:
            next_state = copy.deepcopy(state)
            next_state.next_step(thread)
            temp = next_state.to_string()
            if temp not in state_dictionary.keys() and not next_state.abort:
                state_dictionary[temp] = len(state_dictionary.keys())
                state_list_complete.append(next_state)
                pocket_states.add(next_state)

    print('Number of total states: {}'.format(len(state_list_complete)))
    number_states = {}
    print('Terminate state(s): ')
    for state in state_list_complete:
        terminate = True
        temp = ' '
        for thread in state.threads:
            temp += str(thread['buffer']) + ' '
            if thread['buffer'] != 26:
                terminate = False
        if terminate:
            state.output()
        if temp in number_states.keys():
            number_states[temp] += 1
        else:
            number_states[temp] = 1

    print('Number of none empty buffer: {}'.format(len(number_states.keys())))
    print('Number of terminate states: {}'.format(number_states[' 26 26 ']))

    label_x = set([])
    label_y = set([])
    abort_numbers = {}
    abort_state = len(state_list_complete)
    transition_matrix = np.zeros((abort_state + 1, abort_state + 1), dtype=float)
    for state in enumerate(state_list_complete):
        for switch in [0, 1]:
            flow_out = copy.deepcopy(state[1])
            flow_out.next_step(switch)
            if flow_out.abort:
                temp = ' ' + str(state[1].threads[0]['buffer']) + ' ' + str(state[1].threads[1]['buffer']) + ' '
                label_x.add(state[1].threads[0]['buffer'])
                label_y.add(state[1].threads[1]['buffer'])
                if temp in abort_numbers.keys():
                    abort_numbers[temp] += 1
                else:
                    abort_numbers[temp] = 1
                transition_matrix[abort_state][state[0]] += 0.5
            else:
                transition_matrix[state_dictionary[flow_out.to_string()]][state[0]] += 0.5
    transition_matrix[abort_state][abort_state] = 1.0

    print('Buffers lead to abort: ')
    print('Buffer      Fully/Partly         number of states not leading to abort')
    for key in abort_numbers.keys():
        if number_states[key] - abort_numbers[key] == 0:
            print(key, '        Fully              ', number_states[key] - abort_numbers[key])
        else:
            print(key, '        Partly             ', number_states[key] - abort_numbers[key])

    # Figure 4 of [Li2020]
    rc('font', size=10)
    rc('text', usetex=True)
    plot_array = np.zeros((26, 26), dtype=int)
    for key in number_states.keys():
        temp = key.split()
        number_list = []
        for string in temp:
            number_list.append(int(string))
        plot_array[number_list[1] - 1][number_list[0] - 1] = number_states[key]

    fig, ax = plt.subplots(1, 1)
    ax.set_aspect(aspect=1)
    cmap = plt.cm.Blues
    cmap.set_bad(color='white')
    c = ax.pcolor(plot_array, edgecolors='k', cmap=cmap, vmax=5)
    ax.set_xlabel('buffer for thread $a$')
    ax.set_ylabel('buffer for thread $b$')
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    label_y.add(1)
    label_y.add(26)
    label_y = list(label_y)
    label_y.sort()
    label_y_string = ['$' + str(y) + '$' for y in label_y]

    label_x.add(1)
    label_x.add(26)
    label_x = list(label_x)
    label_x.sort()
    label_x_string = ['$' + str(x) + '$' for x in label_x]
    ax.set_xticks([x - 0.5 for x in label_x])
    ax.set_xticklabels(label_x_string)
    ax.set_yticks([x - 0.5 for x in label_y])
    ax.set_yticklabels(label_y_string)
    ax.tick_params(axis='both', which='major', labelsize=8)
    for key in ['top', 'right', 'bottom', 'left']:
        ax.spines[key].set_visible(False)
    cbaxes = fig.add_axes([0.05, 0.15, 0.02, 0.7])  # This is the position for the colorbar
    cbar = fig.colorbar(c, cax=cbaxes, ticks=[0, 1, 2, 3, 4, 5])
    cbar.ax.set_ylabel(r'number of states in each buffer')
    cbar.ax.yaxis.set_label_position("left")
    cbar.ax.set_yticklabels(['$0$', '$1$', '$2$', '$3$', '$4$', r'$\ge 5$'])  # vertically oriented colorbar
    fig.savefig('degenerate_table.pdf')

    print('Absorbing Markov chain analysis')
    counter = 1
    non_transient_states = set([])
    for state in state_list_complete:
        print('Processing state {}'.format(counter), end='\r')
        counter += 1
        temp_state_list = [copy.deepcopy(state)]
        exit_state_loop = False
        # If some states fail to reach the absorbing states, try increase the number of iterations
        for i in range(1000):
            if i == 999:
                state.output()
            new_state_list = []
            new_string_list = []
            for temp_state in temp_state_list:
                for direction in thread_list:
                    next_step_state = copy.deepcopy(temp_state)
                    next_step_state.next_step(direction)
                    temp_string = next_step_state.to_string()
                    if temp_string in non_transient_states:
                        non_transient_states.add(state.to_string())
                        exit_state_loop = True
                        break
                    if temp_string not in new_string_list:
                        new_string_list.append(temp_string)
                        new_state_list.append(next_step_state)
                    if not next_step_state.abort:
                        terminate = True
                        for thread in next_step_state.threads:
                            if thread['buffer'] != 26:
                                terminate = False
                        if terminate:
                            non_transient_states.add(state.to_string())
                            exit_state_loop = True
                            break
                    else:
                        non_transient_states.add(state.to_string())
                        exit_state_loop = True
                        break
                if exit_state_loop:
                    break
            if exit_state_loop:
                break
            temp_state_list = new_state_list
    print('Number of transients: {}'.format(len(state_list_complete) - len(non_transient_states)))


if __name__ == '__main__':
    main()
