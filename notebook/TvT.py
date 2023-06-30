import os

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

"""cascade reaction A→P→S"""

k1 = lambda T: 8.976446823e5 * np.exp(-41.4235531 / 8.314e-3 / T)  # 1/s
k2 = lambda T: 4.845990778e10 * np.exp(-112.9733266 / 8.314e-3 / T)  # 1/s
dH1 = -55e6  # J/kmol
dH2 = -71.5e6  # J/kmol

CA0 = 0.3  # kmol/m3
CP0 = 0.0
VR = 0.01  # m3
v = 0.01667  # m3/s
T0 = 283  # K
rho = 100  # kg/m3
Cp = 600  # J/(kg*K)
Tc = 330  # K
UA = 666.67  # J/(s*K)
# UA = 315.


def Q_gen(CA, CP, T):
    return (-dH1 * k1(T) * CA - dH2 * k2(T) * CP) * VR


def Q_rem(T, T0):
    return (T - T0) * v * rho * Cp + UA * (T - Tc)


fig = plt.figure(1, layout='tight')
ax_QT = fig.add_subplot(projection='3d')

Ts = np.linspace(250, 750)
sv = np.linspace(0.001, 2)
sv, Ts = np.meshgrid(sv, Ts, indexing='ij')
Qgs = -dH1 * VR * k1(Ts) * CA0 / (1 + k1(Ts) / sv) - \
      dH2 * k2(Ts) * (sv * VR * CP0 + (k1(Ts) * CA0 * sv * VR / (k1(Ts) + sv))) / (k2(Ts) + sv)
Tcs = (sv * VR * rho * Cp * (Ts - T0) - Qgs) / UA + Ts
contour = ax_QT.contour(sv, Tcs, Ts, levels=500)

ax_QT.set(xlabel='sv', ylabel='Tc', zlabel='T')
plt.show()
