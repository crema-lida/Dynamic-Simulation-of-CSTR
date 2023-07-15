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


plt.rcParams.update({
    'grid.color': 'lightgray',
})

fig = plt.figure(figsize=(10, 6))
ax3d = fig.add_subplot(121, projection='3d', position=(0.01, 0.15, 0.5, 0.7))
ax3d.view_init(20, -50)
ax2d = fig.add_subplot(122, position=(0.6, 0.15, 0.35, 0.7))

Ts = np.linspace(250, 750)
sv = np.linspace(0.001, 2)
sv, Ts = np.meshgrid(sv, Ts, indexing='ij')
Qgs = lambda T, sv: -dH1 * VR * k1(T) * CA0 / (1 + k1(T) / sv) - \
                dH2 * k2(T) * (sv * VR * CP0 + (k1(T) * CA0 * sv * VR / (k1(T) + sv))) / (k2(T) + sv)
Tcs = lambda T, sv: (sv * VR * rho * Cp * (T - T0) - Qgs(T, sv)) / UA + T

ax3d.plot_surface(sv, Tcs(Ts, sv), Ts, cmap='coolwarm')
ax3d.contourf(sv, Tcs(Ts, sv), Ts, zdir='y', offset=800, cmap='coolwarm', zorder=0)
ax3d.set(xlabel='Space Velocity', ylabel='Coolant Temperature (K)', zlabel='Reactor Temperature (K)',
         xticks=np.linspace(0, 2, 5),
         zlim=(250, 750))
ax3d.set_ylim(top=800)
ax3d.set_box_aspect((1, 1.2, 1))

ax2d.contourf(sv, Ts, Tcs(Ts, sv), cmap='coolwarm')
levels = ax3d.get_yticks()[1:-1]
contour = ax2d.contour(sv, Ts, Tcs(Ts, sv),
                       levels=levels,
                       colors='k',
                       linewidths=0.6)
ax2d.clabel(contour, fmt='%.0f', fontsize=11)
ax2d.set(title='$T_{reactor}/T_{coolant}$ Behavior Changes\nas a Function of Space Velocity',
         xlabel='Space Velocity (s$^{-1}$)',
         ylabel='Reactor Temperature (K)',
         ylim=(250, 750))

plt.show(block=False)

ln, = ax3d.plot([2, 2], [300, 300], [250, 800], zorder=20, color='tab:orange')
dots = [ax3d.plot(2, 300, 0, zorder=20, color='tab:red', marker='h')[0] for _ in range(5)]
ln2, = ax2d.plot([2, 2], [250, 800], color='tab:orange')
dots2 = [ax2d.plot(2, 0, color='tab:red', marker='h')[0] for _ in range(5)]

residual = lambda T: Tcs(T, sv) - 300
sv = 2
while sv > 0:
    for i, x0 in enumerate([700, 550, 450, 390, 300]):
        sol = opt.root(residual, x0)
        Ts = sol.x[0]
        if sv < 1.297 and i < 2:
            Ts = 0
        elif sv < 1.417 and 2 <= i < 4:
            Ts = 0
        dots[i].set_data([sv], [300])
        dots[i].set_3d_properties([Ts])
        dots2[i].set_data([sv], [Ts])
        ln.set_data([sv, sv], [300, 300])
        ln.set_3d_properties([250, 800])
        ln2.set_data([sv, sv], [250, 800])

    sv -= 0.04
    fig.canvas.draw_idle()
    fig.canvas.flush_events()
    plt.pause(0.01)
    plt.savefig(f'animation/frame{frame}.jpg', dpi=300)

plt.show(block=True)
