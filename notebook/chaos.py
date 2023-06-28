import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Process, Pipe
from numba import njit
import os
from ode_solve import adaptive_RK

"""two parallel first-order reactions A→C and B→D"""

Da1 = 10.
Da2 = 2.75
beta1 = 0.04
beta2 = 0.0149
alpha = 260.
gamma1 = 25
gamma2 = 25


@njit
def dydt(y):
    T, Y, Z = y
    kY = 1 - Y - Y * Da1 * np.exp(gamma1 * T / (1 + T))
    kZ = 1 - Z - Z * Da2 * np.exp(gamma2 * T / (1 + T))
    kT = alpha * (beta1 * Y * Da1 * np.exp(gamma1 * T / (1 + T)) +
                  beta2 * Z * Da2 * np.exp(gamma2 * T / (1 + T)) - T)
    return np.array([kT, kY, kZ])


def update(conn):
    while True:
        paths, t, dt = conn.recv()
        for i, xyz in enumerate(paths):
            while t[i, -1] - t[i, -2] < 0.001:
                xyz[-1], dt[i], dt_new = adaptive_RK(dydt, xyz[-2], dt=dt[i], tol=1e-8)
                t[i, -1] += dt[i]
                dt[i] = dt_new
        conn.send((paths[:, -1], t[:, -1], dt))


if __name__ == '__main__':
    plt.rcParams.update({
        'grid.color': 'lightgray',
    })

    xlim, ylim, zlim = [0, 0.4], [0, 0.04], [0, 0.09]
    fig = plt.figure(figsize=(8, 6.8), layout='tight')
    ax = fig.add_subplot(projection='3d')
    ax.set(xlabel='Dimensionless Temperature',
           ylabel='\n\nDimensionless\nConcentration (A)',
           zlabel='\nDimensionless\nConcentration (B)',
           xlim=xlim, ylim=ylim, zlim=zlim,
           xticks=np.linspace(*xlim, 5),
           yticks=np.linspace(*ylim, 5),
           zticks=np.linspace(*zlim, 4))
    ax.xaxis.set_pane_color("whitesmoke")
    ax.yaxis.set_pane_color("whitesmoke")
    ax.zaxis.set_pane_color("whitesmoke")
    ax.view_init(30, -150)
    # plt.show(block=False)

    # initial values of T, Y, Z
    N = 3  # 2304
    Ny = 48
    Nz = N // Ny
    paths = np.zeros((N, 1, 3))
    paths[0, 0] = [0.01, 0.01, 0.01]
    paths[1, 0] = [0.01, 0.01, 0.01 + 1e-5]
    paths[2, 0] = [0.01, 0.01, 0.01 - 1e-5]
    # paths[:, 0] = np.array([np.full(Ny * Nz, 0.),
    #                         np.linspace(*ylim, Ny).reshape((1, -1)).repeat(Nz, axis=0).reshape(-1),
    #                         np.linspace(*zlim, Nz).repeat(Ny)]).T
    np.random.shuffle(paths)
    lines = []
    scatters = []
    colors = ['#8dd3c7', '#80b1d3', '#fb8072']
    for i, xyz in enumerate(paths):
        ln, = ax.plot(*xyz.T, linewidth=1 if i < 2 else 0.8, label=f'line {i}', color=colors[i])
        lines.append(ln)
    scatter = ax.scatter(*paths[:, -1].T, color=colors)
    # scatter = ax.scatter(*paths[:, -1].T, color=colors[1], s=1)

    N_processes = 1
    n = N // N_processes
    pipes = []
    processes = []
    for i in range(N_processes):
        conn1, conn2 = Pipe()
        pipes.append(conn1)
        processes.append(Process(target=update, args=(conn2,), daemon=True))
        processes[i].start()

    t = np.zeros((N, 1))
    dt = np.ones(N)
    frame = 0
    t_max = 2.0

    os.makedirs('animation', exist_ok=True)

    while t[0, -1] < t_max:
        paths = np.append(paths, np.zeros((N, 1, 3)), axis=1)
        t = np.append(t, t[:, -1:], axis=1)

        for i, conn in enumerate(pipes):
            s = np.s_[i * n: i * n + n]
            conn.send((paths[s, -2:], t[s, -2:], dt[s]))
        for i, conn in enumerate(pipes):
            s = np.s_[i * n: i * n + n]
            paths[s, -1], t[s, -1], dt[s] = conn.recv()

        print(f'\r t = {t[0, -1]:.2f}', end='')
        for i, xyz in enumerate(paths):
            T, Y, Z = xyz.T
            lines[i].set_data(T, Y)
            lines[i].set_3d_properties(Z)
        scatter._offsets3d = paths[:, -1].T

        # fig.canvas.draw_idle()
        # fig.canvas.flush_events()
        plt.savefig(f'animation/lines_{frame}.jpg', dpi=200)
        frame += 1

    plt.show(block=True)
