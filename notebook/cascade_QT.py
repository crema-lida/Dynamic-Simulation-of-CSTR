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
UA = 600.  # J/(s*K)
# UA = 500.


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


fig = plt.figure(layout='tight')
ax = fig.add_subplot()

t_max = 1400
dt = 1.0

"""initial values"""
t = np.zeros(1)
CA = np.array([0.2956])
CP = np.array([0.0044])
CS = np.array([0.])
T = np.array([286.2806])
T0_arr = np.array([290.])
Q = np.zeros(1)

Tx = np.linspace(250, 1000, 100)
Q_gen_steady = -dH1 * VR * k1(Tx) * CA0 / (1 + k1(Tx) * VR / v) - \
               dH2 * VR * k2(Tx) * (v * CP0 + (k1(Tx) * CA0 * v * VR / (k1(Tx) * VR + v))) / (k2(Tx) * VR + v)
ax.plot(Tx, Q_gen_steady, color='#FF731D', label='$Q_{gen}$')
ln_Q_rem, = ax.plot(Tx, Q_rem(Tx, T0_arr[-1]), color='#5F9DF7', label='$Q_{rem}$')
ln_Q_gen_transient, = ax.plot(T, Q, ':', color='darkcyan', label='$Q_{gen}$ (transient)')
dot_Q, = ax.plot(0, '-d', linewidth=0.6, color='darkcyan',
                 markersize=8, fillstyle='left',
                 markerfacecolor='darkcyan', markerfacecoloralt='lightcyan',
                 markeredgecolor='darkcyan', markeredgewidth='0.5')
ax.legend(loc='lower right')
ax.set(xlabel='T (K)', ylabel='Q (J/s)',
       xlim=(250, 950), ylim=(0, 1.25e6))
ax.ticklabel_format(axis='y', scilimits=[-3, 3])

plt.show(block=False)
frame = 0
t_interval = 0

DIR = 'animation/UA_600'
# DIR = 'animation/UA_500'
os.makedirs(DIR, exist_ok=True)

while t[-1] < t_max:
    y0 = np.array([CA[-1], CP[-1], CS[-1], T[-1]])
    y, dt, dt_new = adaptive_RK(dydt, y0, args=(T0_arr[-1],), dt=dt, tol=1e-10)
    CA_new, CP_new, CS_new, T_new = y

    t = np.append(t, t[-1] + dt)
    CA = np.append(CA, CA_new)
    CP = np.append(CP, CP_new)
    CS = np.append(CS, CS_new)
    T = np.append(T, T_new)
    Q = np.append(Q, Q_gen(CA_new, CP_new, T_new))
    T0_arr = np.append(T0_arr, T0_arr[-1] + 2 / (60 / dt))  # T0 increases 2 K every 60 seconds

    if 150 < t[-1] < 1225:
        if t_interval <= 0:
            ln_Q_rem.set_data(Tx, Q_rem(Tx, T0_arr[-1]))
            ln_Q_gen_transient.set_data(T, Q)
            dot_Q.set_data([T[-1], T[-1]], [-Q[-1], Q[-1]])
            fig.canvas.draw_idle()
            fig.canvas.flush_events()
            if t[-1] < 157:
                t_interval = 0.02
            elif t[-1] < 1115:
                t_interval = 60
            else:
                t_interval = 0.02
            plt.savefig(f'{DIR}/frame{frame}.jpg', dpi=300)
            frame += 1
        else:
            t_interval -= dt

    dt = dt_new
    print(f'\rt = {t[-1]:.2f} s, CA = {CA[-1]:.3f} kmol/m3, CP = {CP[-1]:.3f} kmol/m3, '
          f'T = {T[-1]:.2f} K, dt = {dt:.2e} s, frame = {frame}', end='')

plt.show(block=True)
