import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

# cascade reaction A→P→S
k1 = lambda T: 8.976446823e5 * np.exp(-41.4235531 / 8.314e-3 / T)  # 1/s
k2 = lambda T: 4.845990778e10 * np.exp(-112.9733266 / 8.314e-3 / T)

VR = 0.01  # m3
CA0 = 0.3  # kmol/m3
CP0 = 0.0
v = 0.01667  # m3/s
T0 = 250  # K
dH1 = -55e6  # J/kmol
dH2 = -71.5e6
rho = 100  # kg/m3
Cp = 600  # J/(kg*K)
Tc = 333.15
UA = 666.67  # J/(s*K)


def kCA(CA, CP, T):
    return (CA0 - CA) * v / VR - k1(T) * CA


def kCP(CA, CP, T):
    return (CP0 - CP) * v / VR - k2(T) * CP + k1(T) * CA


def kT(CA, CP, T):
    return (-dH1 * k1(T) * CA - dH2 * k2(T) * CP) / (rho * Cp) + \
        (T0 - T) * v / VR + UA * (Tc - T) / (VR * rho * Cp)


"""methods"""


def RK4(CAi, CPi, Ti, dt):
    kCA1 = kCA(CAi, CPi, Ti)
    kCP1 = kCP(CAi, CPi, Ti)
    kT1 = kT(CAi, CPi, Ti)
    args = CAi + dt / 2 * kCA1, CPi + dt / 2 * kCP1, Ti + dt / 2 * kT1
    kCA2 = kCA(*args)
    kCP2 = kCP(*args)
    kT2 = kT(*args)
    args = CAi + dt / 2 * kCA2, CPi + dt / 2 * kCP2, Ti + dt / 2 * kT2
    kCA3 = kCA(*args)
    kCP3 = kCP(*args)
    kT3 = kT(*args)
    args = CAi + dt * kCA3, CPi + dt * kCP3, Ti + dt * kT3
    kCA4 = kCA(*args)
    kCP4 = kCP(*args)
    kT4 = kT(*args)

    CA_new = CAi + dt / 6 * (kCA1 + 2 * kCA2 + 2 * kCA3 + kCA4)
    CP_new = CPi + dt / 6 * (kCP1 + 2 * kCP2 + 2 * kCP3 + kCP4)
    T_new = Ti + dt / 6 * (kT1 + 2 * kT2 + 2 * kT3 + kT4)
    return CA_new, CP_new, T_new


def BS23(CAi, CPi, Ti, dt, tol=1e-6):
    kCA1 = kCA(CAi, CPi, Ti)
    kCP1 = kCP(CAi, CPi, Ti)
    kT1 = kT(CAi, CPi, Ti)
    args = CAi + 0.5 * dt * kCA1, CPi + 0.5 * dt * kCP1, Ti + 0.5 * dt * kT1
    kCA2 = kCA(*args)
    kCP2 = kCP(*args)
    kT2 = kT(*args)
    args = CAi + 0.75 * dt * kCA2, CPi + 0.75 * dt * kCP2, Ti + 0.75 * dt * kT2
    kCA3 = kCA(*args)
    kCP3 = kCP(*args)
    kT3 = kT(*args)

    CA_new = CAi + dt * (2 * kCA1 + 3 * kCA2 + 4 * kCA3) / 9
    CP_new = CPi + dt * (2 * kCP1 + 3 * kCP2 + 4 * kCP3) / 9
    T_new = Ti + dt * (2 * kT1 + 3 * kT2 + 4 * kT3) / 9

    kCA4 = kCA(CA_new, CP_new, T_new)
    kCP4 = kCP(CA_new, CP_new, T_new)
    kT4 = kT(CA_new, CP_new, T_new)

    err = np.array([
        dt * (-5 * kCA1 / 72 + kCA2 / 12 + kCA3 / 9 - kCA4 / 8),
        dt * (-5 * kCP1 / 72 + kCP2 / 12 + kCP3 / 9 - kCP4 / 8),
        dt * (-5 * kT1 / 72 + kT2 / 12 + kT3 / 9 - kT4 / 8),
    ])

    E = np.linalg.norm(err, np.inf)
    maxerr = tol * (1 + np.linalg.norm([CAi, CPi, Ti], np.inf))
    q = min(0.9 * (maxerr / E) ** 0.25, 5.0)
    if E > maxerr:
        return BS23(CAi, CPi, Ti, dt * q)
    else:
        return CA_new, CP_new, T_new, dt * q


def DIRK(CAi, CPi, Ti, dt):
    args = (CAi, CPi, Ti)
    K0 = np.array([kCA(*args), kCP(*args), kT(*args)])

    def residual_K1(K1):
        args = []
        for i, y in enumerate((CAi, CPi, Ti)):
            args.append(y + 0.25 * dt * K1[i])
        return [kCA(*args) - K1[0], kCP(*args) - K1[1], kT(*args) - K1[2]]

    def residual_K2(K2, K1):
        args = []
        for i, y in enumerate((CAi, CPi, Ti)):
            args.append(y + 0.5 * dt * K1[i] + 0.25 * dt * K2[i])
        return [kCA(*args) - K2[0], kCP(*args) - K2[1], kT(*args) - K2[2]]

    sol = opt.root(residual_K1, K0, method='hybr')
    K1 = sol.x
    if not sol.success:
        return DIRK(CAi, CPi, Ti, 0.1 * dt)
    sol = opt.root(residual_K2, K1, args=(sol.x), method='hybr')
    K2 = sol.x
    if not sol.success:
        return DIRK(CAi, CPi, Ti, 0.1 * dt)
    CA_new = CAi + dt * (0.5 * K1[0] + 0.5 * K2[0])
    CP_new = CPi + dt * (0.5 * K1[1] + 0.5 * K2[1])
    T_new = Ti + dt * (0.5 * K1[2] + 0.5 * K2[2])
    return CA_new, CP_new, T_new


t_max = 3000
dt = 1.0

t = np.zeros(1)
CA = np.array([0.2956])
CP = np.array([0.0044])
T = np.array([286.2806])
T0_arr = np.array([250.])
Q = np.zeros(1)

cooldown = 0

while t[-1] < t_max:
    if CA[-1] < 0 or CP[-1] < 0 or np.isnan([CA[-1]]):
        # raise Exception('Solution Diverged')
        dt *= 0.1
        t = np.delete(t, np.s_[-100:])
        CA = np.delete(CA, np.s_[-100:])
        CP = np.delete(CP, np.s_[-100:])
        T = np.delete(T, np.s_[-100:])
        Q = np.delete(Q, np.s_[-100:])
        T0_arr = np.delete(T0_arr, np.s_[-100:])
        cooldown += 200
    if cooldown > 0:
        cooldown -= 1
        if cooldown % 200 == 0:
            dt *= 10

    CAi, CPi, Ti, T0 = CA[-1], CP[-1], T[-1], T0_arr[-1]

    """RK4"""

    # CA_new, CP_new, T_new = RK4(CAi, CPi, Ti, dt)

    """Bogacki-Shampine"""

    # CA_new, CP_new, T_new, dt = BS23(CAi, CPi, Ti, dt)

    """diagonally implicit RK"""

    CA_new, CP_new, T_new = DIRK(CAi, CPi, Ti, dt)

    t = np.append(t, t[-1] + dt)
    CA = np.append(CA, CA_new)
    CP = np.append(CP, CP_new)
    T = np.append(T, T_new)
    Q = np.append(Q, -dH1 * VR * k1(T_new) * CA_new - dH2 * VR * k2(T_new) * CP_new)
    T0_arr = np.append(T0_arr, T0 + 200 / (6000 / dt))

    print(f't = {t[-1]:.2f} s, CA = {CA[-1]:.3f} kmol/m3, CP = {CP[-1]:.3f} kmol/m3, T = {T[-1]:.2f} K, dt = {dt:.2e}')

fig = plt.figure(figsize=(8, 8))
axQ = fig.add_subplot(211)
axT = fig.add_subplot(212)

tau = VR / v

Tx = np.linspace(290, 800, 100)
Q_gen = -dH1 * VR * k1(Tx) * CA0 / (1 + k1(Tx) * VR / v) - \
        dH2 * VR * k2(Tx) * k1(Tx) * CA0 * tau / (k1(Tx) * tau + 1) / (k2(Tx) * tau + 1)
Q_rem = v * rho * Cp * (Tx - T0) + UA * (Tx - Tc)
axQ.plot(Tx, Q_gen, color='tomato')
axQ.plot(Tx, Q_rem, color='skyblue')
axQ.plot(T, Q, ':')
axQ.set(xlabel='T (K)', ylabel='Q (J/s)')
axQ.set_ylim(0, 1e6)

axT.plot(t, T, color='red', label='Tr')
axT.plot(t, T0_arr, color='blue', label='T0')
axC = axT.twinx()
axC.plot(t, CA, color='purple', label='$c_A$')
axC.plot(t, CP, color='green', label='$c_P$')
axT.legend(loc='upper left')
axC.legend(loc='lower left')
axT.set(xlabel='t (s)', ylabel='T (K)',
        xlim=(0, t_max))
axC.set(ylabel='Concentration (kmol/m3)')

plt.show()
