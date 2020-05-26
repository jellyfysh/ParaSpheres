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

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mtick
import matplotlib
import sys


sample_size = int(sys.argv[2])
number_points = int(sys.argv[1])
filename = "EPH.dat"
x = np.arange(1, number_points + 1, 1)
j = 0
samples = []
i = 0
samples_temp = []
with open(filename) as file:
    for line in file:
        i += 1
        samples_temp.append(float(line[5:]))
        if i == sample_size:
            j += 1
            i = 0
            samples.append(samples_temp)
            samples_temp = []
EPH = []
EPH_median = []
errorbar = []

for k in range(j):
    EPH.append(np.mean(samples[k]))
    EPH_median.append(np.median(samples[k]))
    errorbar.append(np.sqrt(np.var(samples[k])))
    print(EPH[k])

matplotlib.rcParams.update({'font.size': 12})

fig, ax1 = plt.subplots()
ax1.set_xlabel('# of physical threads (OMP_NUM_THREADS)')
ax1.set_ylabel('Absolute speed (events/hour)')

bp1 = ax1.boxplot(samples, whis=100)
ax1.set_xlim(0, x[-1] + 1)
ax1.ticklabel_format(axis='y', style='sci', useMathText=True, scilimits=(0, 0))
ax1.set_ylim(0, EPH[-1]*1.1)
ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
ax1.spines["top"].set_visible(False)
ax1.tick_params(axis='y')

ax2 = ax1.twinx()

ax2.plot(x, [e / EPH[0] for e in EPH], alpha=0)
ax2.set_ylabel('Relative speed')
ax2.tick_params(axis='y')
ax2.spines["top"].set_visible(False)
ax2.set_ylim(0, EPH[-1]*1.1 / EPH[0])
ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
fig.tight_layout(pad=0.1, w_pad=0.1, h_pad=0.1)
fig.savefig('EPH.pdf')
