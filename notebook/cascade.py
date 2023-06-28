import os

import numpy as np
import matplotlib.pyplot as plt
from ode_solve import adaptive_RK

"""cascade reaction A→P→S"""

k1 = lambda T: 8.976446823e5 * np.exp(-41.4235531 / 8.314e-3 / T)  # 1/s
k2 = lambda T: 4.845990778e10 * np.exp(-112.9733266 / 8.314e-3 / T)  # 1/s
dH1 = -55e6  # J/kmol
dH2 = -71.5e6  # J/kmol

CA0 = 0.3  # kmol/m3
CP0 = 0.0
CS0 = 0.0
VR = 0.01  # m3
v = 0.01667  # m3/s
T0 = 250  # K
rho = 100  # kg/m3
Cp = 600  # J/(kg*K)
Tc = 330.
# UA = 600.  # J/(s*K)
UA = 315.


def Q_gen(CA, CP, T):
    return (-dH1 * k1(T) * CA - dH2 * k2(T) * CP) * VR


def Q_rem(T, T0):
    return (T - T0) * v * rho * Cp + UA * (T - Tc)


def dydt(y, T0):
    CA, CP, CS, T = y
    kCA = (CA0 - CA) * v / VR - k1(T) * CA
    kCP = (CP0 - CP) * v / VR - k2(T) * CP + k1(T) * CA
    kCS = (CS0 - CS) * v / VR + k2(T) * CP
    kT = (Q_gen(CA, CP, T) - Q_rem(T, T0)) / (VR * rho * Cp)
    return np.array([kCA, kCP, kCS, kT])


t_max = 1500
dt = 1.0

"""initial values"""
t = np.zeros(1)
CA = np.array([0.2956])
CP = np.array([0.0044])
CS = np.array([0.])
T = np.array([286.2806])
T0_arr = np.array([250.])
Q_gen_transient = np.zeros(1)
Q_rem_transient = np.zeros(1)

while t[-1] < t_max:
    y0 = np.array([CA[-1], CP[-1], CS[-1], T[-1]])
    y, dt, dt_new = adaptive_RK(dydt, y0, args=(T0_arr[-1],), dt=dt, tol=1e-8)
    CA_new, CP_new, CS_new, T_new = y

    t = np.append(t, t[-1] + dt)
    CA = np.append(CA, CA_new)
    CP = np.append(CP, CP_new)
    CS = np.append(CS, CS_new)
    T = np.append(T, T_new)
    Q_gen_transient = np.append(Q_gen_transient, Q_gen(CA_new, CP_new, T_new))
    Q_rem_transient = np.append(Q_rem_transient, Q_rem(T_new, T0_arr[-1]))
    T0_arr = np.append(T0_arr, T0_arr[-1] + 2 / (60 / dt))  # T0 increases 2 K every 60 seconds

    print(f'\rt = {t[-1]:.2f} s, CA = {CA[-1]:.3f} kmol/m3, CP = {CP[-1]:.3f} kmol/m3, '
          f'T = {T[-1]:.2f} K, dt = {dt:.2e} s', end='')
    dt = dt_new

"""plots"""

text_color = '#333333'
plt.rcParams.update({
    'text.color': text_color,
    'axes.labelcolor': text_color,
    'axes.edgecolor': text_color,
    'xtick.color': text_color,
    'ytick.color': text_color,
})

fig_steady = plt.figure(1, figsize=(6.4, 4.8), layout='tight')
fig_transient = plt.figure(2, figsize=(6.4, 4.8), layout='tight')
fig_Qt = plt.figure(3, figsize=(6.4, 3.2), layout='tight')
ax_QT = fig_steady.add_subplot()
ax_Ct = fig_transient.add_subplot(211)
ax_Tt = fig_transient.add_subplot(212)
ax_Qt = fig_Qt.add_subplot()

Tx = np.linspace(250, 1000, 100)
Q_gen_steady = -dH1 * VR * k1(Tx) * CA0 / (1 + k1(Tx) * VR / v) - \
               dH2 * VR * k2(Tx) * (v * CP0 + (k1(Tx) * CA0 * v * VR / (k1(Tx) * VR + v))) / (k2(Tx) * VR + v)
# UA = 600
UA = 315
Q_rem_steady = Q_rem(Tx, T0)

ax_QT.plot(Tx, Q_gen_steady, color='#FF731D', label='$Q_{gen}$')
ax_QT.plot(Tx, Q_rem_steady, color='#5F9DF7', label='$Q_{rem}$')
ax_QT.plot(T, Q_gen_transient, ':', color='darkcyan', label='$Q_{gen}$ (transient)')
ax_QT.legend(loc='lower right')
ax_QT.set(xlabel='T (K)', ylabel='Q (J/s)',
          xlim=(250, 950), ylim=(0, 1.7e6))
ax_QT.ticklabel_format(axis='y', scilimits=[-3, 3])

ax_Tt.plot(t, T, color='#F45050', label='$T_r$')
ax_Tt.plot(t, T0_arr, color='#0079FF', label='$T_0$')
ax_Tt.legend(loc='upper left')
ax_Tt.set(xlabel='t (s)', ylabel='T (K)',
          xlim=(1345, 1355))

ax_Ct.plot(t, CA, color='#FD8A8A', label='$c_A$')
ax_Ct.plot(t, CP, color='#9EA1D4', label='$c_P$')
ax_Ct.plot(t, CS, color='#A8D1D1', label='$c_S$')
ax_Ct.legend(loc='upper left')
ax_Ct.set(xlabel='t (s)', ylabel='Concentration (kmol/m$^3$)',
          xlim=(1345, 1355))

# xlim, ylim = (1355, 1356), (4e4, 4e6)
# xlim, ylim = (2323.8, 2324.4), (2e5, 1e8)
xlim, ylim = (1349, 1350), (5e4, 1e9)
ax_Qt.semilogy(t, Q_gen_transient, color='#FF731D', label='$Q_{gen}$')
ax_Qt.semilogy(t, Q_rem_transient, color='#5F9DF7', label='$Q_{rem}$')
ax_Qt.legend(loc='upper left')
ax_Qt.set(xlabel='t (s)', ylabel='Q (J/s)',
          xlim=xlim, ylim=ylim)
ax_Qt.ticklabel_format(axis='x', useOffset=False)

# # DIR = 'cascade/UA_600'
# DIR = 'cascade/UA_315'
# os.makedirs(DIR, exist_ok=True)
# for i, name in enumerate(['QT', 'C-t', 'Q-t'], start=1):
#     plt.figure(i)
#     plt.savefig(f'{DIR}/{name}.jpg', dpi=600, bbox_inches='tight')

plt.show()
