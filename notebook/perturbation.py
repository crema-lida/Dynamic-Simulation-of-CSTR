import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axisartist.axislines import AxesZero

fig = plt.figure()
ax = fig.add_subplot(axes_class=AxesZero)
for direction in ["xzero", "yzero"]:
    # adds arrows at the ends of each axis
    ax.axis[direction].set_axisline_style("-|>")

    # adds X and Y-axis from the origin
    ax.axis[direction].set_visible(True)

for direction in ["left", "right", "bottom", "top"]:
    # hides borders
    ax.axis[direction].set_visible(False)

colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown']
ends = [-1, -1, -1, -1, 6, 40]

"""exponential"""
# t = np.linspace(0, 10, 100)
# A = -np.exp(-1 * t) + 2 * np.exp(-0.2 * t)
# B = (np.exp(-0.2 * t) + np.exp(-1.5 * t)) / 2
# C = - np.exp(-0.2 * t) + 2 * np.exp(-2 * t)
# D = (np.exp(-0.05 * t) + np.exp(-0.05 * t)) / 2
# E = (np.exp(t) + np.exp(t)) / 2
# F = -np.exp(0.2 * t) + 2 * np.exp(-0.1 * t)
#
# for x, label, color, end in zip([A, B, C, D, E, F], ['A', 'B', 'C', 'D', 'E', 'F'], colors, ends):
#     ax.plot(t[:end], x[:end], color=f'tab:{color}')
#     ax.annotate("", xy=(t[end], x[end]), xytext=(t[end - 1], x[end - 1]),
#                 arrowprops=dict(width=0, headwidth=4, headlength=8, color=f'tab:{color}'))
#     ax.annotate(label, xy=(t[end] + 0.1, x[end] - 0.05), fontsize=13)
#
# ax.set(xlim=(0, 10.5), ylim=(-1, 2), xticks=[], yticks=[],
#        ylabel='x or y')
# ax.annotate('t', xy=(1, 0.31), ha='left', va='top',
#             xycoords='axes fraction')
# plt.savefig('figures/pertub_exp.jpg', dpi=300)

"""sine"""
t = np.linspace(0, 10, 100)
A = 0.5 * np.exp(0.1 * t) * (np.cos(t) + np.sin(t))
B = 0.5 * (np.cos(t) + np.sin(t))
C = 0.5 * np.exp(-0.2 * t) * (np.cos(2 * t) + np.sin(2 * t))

for x, label, color, end in zip([A, B, C], ['A', 'B', 'C'], colors, ends):
    ax.plot(t[:end], x[:end], color=f'tab:{color}')
    ax.annotate("", xy=(t[end], x[end]), xytext=(t[end - 1], x[end - 1]),
                arrowprops=dict(width=0, headwidth=4, headlength=8, color=f'tab:{color}'))
    ax.annotate(label, xy=(t[end] + 0.1, x[end] - 0.05), fontsize=13)

ax.set(xlim=(0, 10.5), ylim=(-2, 2), xticks=[], yticks=[],
       ylabel='x or y')
ax.annotate('t', xy=(1, 0.48), ha='left', va='top',
            xycoords='axes fraction')

# plt.savefig('figures/pertub_sine.jpg', dpi=300)

plt.show()
