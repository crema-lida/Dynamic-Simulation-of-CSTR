import os

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt


k1 = lambda T: 8.976446823e5 * np.exp(-41.4235531 / 8.314e-3 / T)  # 1/s
k2 = lambda T: 4.845990778e10 * np.exp(-112.9733266 / 8.314e-3 / T)  # 1/s
dH1 = -55e6  # J/kmol
dH2 = -71.5e6  # J/kmol

CA0 = 0.3  # kmol/m3
CP0 = 0.0
VR = 0.01  # m3
v = 0.01667  # m3/s
T0 = 280  # K
rho = 100  # kg/m3
Cp = 600  # J/(kg*K)
Tc = 330  # K
UA = 666.67  # J/(s*K)
# UA = 315.


def Q_gen(CA, CP, T):
    return (-dH1 * k1(T) * CA - dH2 * k2(T) * CP) * VR


def Q_rem(T, T0):
    return (T - T0) * v * rho * Cp + UA * (T - Tc)


text_color = '#333333'
plt.rcParams.update({
    'text.color': text_color,
    'font.size': 14,
    'axes.labelcolor': text_color,
    'axes.edgecolor': text_color,
    'xtick.color': text_color,
    'ytick.color': text_color,
})

fig = plt.figure(1, layout='tight')
ax = fig.add_subplot()

Tx = np.linspace(250, 1000, 100)
# Q_gen_steady = lambda T: -dH1 * VR * k1(T) * CA0 / (1 + k1(T) * VR / v)
Q_gen_steady = lambda T: -dH1 * VR * k1(T) * CA0 / (1 + k1(T) * VR / v) - \
               dH2 * VR * k2(T) * (v * CP0 + (k1(T) * CA0 * v * VR / (k1(T) * VR + v))) / (k2(T) * VR + v)
Q_rem_steady = Q_rem(Tx, T0)
ax.plot(Tx, Q_gen_steady(Tx), color='#FF731D', label='$Q_{gen}$')
ax.plot(Tx, Q_rem_steady, color='#5F9DF7', label='$Q_{rem}$')

residual = lambda T: Q_gen_steady(T) - Q_rem(T, T0)
TA = opt.root(residual, x0=300).x[0]
ax.scatter(TA, Q_rem(TA, T0), color='tab:blue', zorder=2)
TB = opt.root(residual, x0=360).x[0]
ax.scatter(TB, Q_rem(TB, T0), color='tab:red', zorder=2)
TC = opt.root(residual, x0=450).x[0]
ax.scatter(TC, Q_rem(TC, T0), color='tab:blue', zorder=2)
TD = opt.root(residual, x0=580).x[0]
ax.scatter(TD, Q_rem(TD, T0), color='tab:red', zorder=2)
TE = opt.root(residual, x0=700).x[0]
ax.scatter(TE, Q_rem(TE, T0), color='tab:blue', zorder=2)
ax.annotate('A', xy=(TA - 20, Q_rem(TA, T0) + 1e4))
ax.annotate('B', xy=(TB - 20, Q_rem(TB, T0) + 1e4))
ax.annotate('C', xy=(TC - 20, Q_rem(TC, T0) + 1e4))
ax.annotate('D', xy=(TD - 20, Q_rem(TD, T0) + 1e4))
ax.annotate('E', xy=(TE - 20, Q_rem(TE, T0) + 1e4))

ax.legend(loc='lower right')
ax.set(xlabel='T (K)', ylabel='Q (J/s)',
       xlim=(250, 750), ylim=(0, 8e5))
ax.ticklabel_format(axis='y', scilimits=[-3, 3])

# plt.savefig(f'QT2.jpg', dpi=600, bbox_inches='tight')

plt.show()
