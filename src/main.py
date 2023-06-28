import json
import os

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox, CheckButtons
import math
import time

from CSTR import CSTR

text_color = '#333333'
bg_color = 'white'
ln_C_colors = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
               '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
               '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f',
               '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
theme_color = '#79B4B7'
silver = '#D8D8D8'
title = 'Dynamic Simulation of CSTR'
reactions_dir = '../reactions'
scripts_dir = '../scripts'

plt.rcParams.update({
    'backend': 'TkAgg',
    'toolbar': 'None',
    'text.color': text_color,
    'axes.labelcolor': text_color,
    'axes.edgecolor': text_color,
    'xtick.color': text_color,
    'ytick.color': text_color,
})
# a window for controlling variables
fig_setvar = plt.figure(2, figsize=(4, 6), facecolor=bg_color)
plt.get_current_fig_manager().window.wm_geometry('+20+20')
ax_check_btn = fig_setvar.add_axes([0.1, 0.1, 0.8, 0.8])  # [left, bottom, width, height]
# main figure window
fig = plt.figure(1, figsize=(11, 6.8), facecolor=bg_color)
fig.text(0.5, 0.945, title, ha='center', va='bottom', fontweight='light', fontsize=16)
ax_Tt = fig.add_axes([0.06, 0.64, 0.57, 0.28], facecolor=(0, 0, 0, 0), zorder=10)
ax_QT = fig.add_axes([0.06, 0.08, 0.30, 0.36], facecolor=bg_color)
ax_Ct = ax_Tt.twinx()
ax_T0_slider = fig.add_axes([0.22, 0.54, 0.32, 0.009])
ax_Tc_slider = fig.add_axes([0.22, 0.50, 0.32, 0.01])
ax_v_slider = fig.add_axes([0.39, 0.07, 0.008, 0.33])
ax_UA_slider = fig.add_axes([0.44, 0.07, 0.008, 0.33])
ax_radio_btn = fig.add_axes([0.49, 0.08, 0.17, 0.36], facecolor=bg_color)
ax_tolerance_input = fig.add_axes([0.92, 0.16, 0.06, 0.04])
ax_speed_input = fig.add_axes([0.92, 0.1, 0.06, 0.04])
ax_reset_btn = fig.add_axes([0.92, 0.04, 0.06, 0.04])

for ax in (ax_radio_btn, ax_tolerance_input, ax_speed_input, ax_reset_btn, ax_check_btn):
    plt.setp(ax.spines.values(), color=silver)
ax_radio_btn.text(0.54, 0.44, 'Reactions', transform=fig.transFigure, ha='center', va='center',
                  backgroundcolor=bg_color)

# transient curves
ln_Tt, = ax_Tt.plot(0, '--', linewidth=0.8, color='red')
ln_Ct = []
ax_Tt.set(xlabel='t (s)', ylabel='T (K)')
ax_Ct.set(ylabel='Concentration (kmol/m$^3$)')
axT_legend = ax_Tt.legend([ln_Tt], [f'Reactor Temperature: - K'],
                          loc='upper left')
axC_legend = None

# Q-T curves
ln_Q_gen, = ax_QT.plot(0, label='$Q_{gen}$', color='#FF731D')
ln_Q_rem, = ax_QT.plot(0, label='$Q_{rem}$', color='#5F9DF7')
ln_QT, = ax_QT.plot(0, ':', color='darkcyan')
dot_Q, = ax_QT.plot(0, '-d', linewidth=0.6, color='darkcyan',
                    markersize=8, fillstyle='left',
                    markerfacecolor='darkcyan', markerfacecoloralt='lightcyan',
                    markeredgecolor='darkcyan', markeredgewidth='0.5')
ax_QT.set(xlabel='T (K)', ylabel='Q (J/s)')
ax_QT.ticklabel_format(axis='y', scilimits=[-3, 3])
ax_QT.legend(loc='upper left')

# Tr-v-Tc/T0 graph
ax_TvT = fig.add_axes([0.74, 0.36, 0.24, 0.48], facecolor=bg_color)

state = {
    't': NDArray,
    'C': NDArray,
    'T': NDArray,
    'Q': NDArray,
    'T_lim': [270, 800],
    'Q_lim': [0, 0],
    'prompt_start': 0.,
    'reaction': '',
    't_max': 10,
    'tol': 1e-10,
    'contour': 'Tc',
    'params': {},
    'scripts': {},
}
reactor: CSTR
reaction_list = [fname.replace('.json', '') for fname in sorted(os.listdir(reactions_dir))]
prompt = fig.text(0.5, 0.01, '', ha='center', va='bottom', fontweight='light', fontsize=9, color='red')


def draw_contour():
    ax_TvT.clear()
    z = reactor.get_2d_data(state['contour'])
    ax_TvT.contourf(*reactor.mgrid, z,
                    cmap='RdBu_r')
    levels = ax_TvT.get_yticks()[1:-1]
    if len(levels) > 5:
        levels = levels[::2]
    contour = ax_TvT.contour(*reactor.mgrid, z,
                             levels=levels,
                             colors=(text_color,),
                             linewidths=0.5)
    ax_TvT.clabel(contour, fmt='%.0f', fontsize=9)
    ax_TvT.set(title='$T_{reactor}/T_{coolant}$ Behavior Changes\nas a Function of Space Velocity',
               xlabel='Space Velocity (s$^{-1}$)',
               ylabel='Reactor Temperature (K)')


def update_T0(val):
    reactor.T0 = val
    T = np.linspace(*state['T_lim'], 10)
    ln_Q_rem.set_data(T, reactor.Q_rem(T))
    if reactor.UA > 0:
        draw_contour()
    if T0_slider.drag_active:
        update_data()


def update_Tc(val):
    reactor.Tc = val
    T = np.linspace(*state['T_lim'], 10)
    ln_Q_rem.set_data(T, reactor.Q_rem(T))
    if Tc_slider.drag_active:
        update_data()


def update_UA(val):
    reactor.UA = val
    T = np.linspace(*state['T_lim'], 10)
    ln_Q_rem.set_data(T, reactor.Q_rem(T))
    if reactor.UA > 0:
        draw_contour()
    if UA_slider.drag_active:
        update_data()


def update_v(val):
    reactor.v = val
    T = np.linspace(*state['T_lim'])
    Q_gen = reactor.Q_gen_steady(T, reactor.v / reactor.VR)
    ln_Q_gen.set_data(T, Q_gen)
    ln_Q_rem.set_data(T, reactor.Q_rem(T))
    state['Q_lim'] = [np.min(Q_gen), np.max(Q_gen) * 1.5]
    if v_slider.drag_active:
        update_data()


def select_reaction(label):
    if label == state['reaction']:
        return
    with open(f'{reactions_dir}/{label}.json') as f:
        try:
            state['params'] = json.load(f)
        except json.decoder.JSONDecodeError as err:
            prompt.set_text(err)
            state['prompt_start'] = time.time()
    initialize()
    state['reaction'] = label

    idx = reaction_list.index(label)
    edgecolor = [silver] * len(reaction_list)
    edgecolor[idx] = theme_color
    facecolor = [bg_color] * len(reaction_list)
    facecolor[idx] = theme_color
    radio_btn.set_radio_props({'edgecolor': edgecolor, 'facecolor': facecolor})


def reset():
    if 'C' in reactor.initial:
        state['C'] = np.array(reactor.initial['C'], float).reshape((1, -1))
    else:
        state['C'] = np.zeros_like(reactor.C0).reshape((1, -1))
    if 'T' in reactor.initial:
        state['T'] = np.array([reactor.initial['T']], float)
    else:
        state['T'] = np.array([reactor.T0], float)
    state['t'] = np.zeros(1)
    state['Q'] = np.zeros(1)
    reactor.T0 = T0_slider.valinit
    reactor.Tc = Tc_slider.valinit
    reactor.UA = UA_slider.valinit
    reactor.v = v_slider.valinit
    T0_slider.reset()
    Tc_slider.reset()
    UA_slider.reset()
    v_slider.reset()

    if f'scripts.json' not in os.listdir(scripts_dir):
        prompt.set_text('scripts.json not found.')
        state['prompt_start'] = time.time()
        return
    with open(f'{scripts_dir}/scripts.json') as f:
        state['scripts'] = json.load(f)
    ax_check_btn.clear()
    check_btn.__init__(
        ax=ax_check_btn,
        labels=list(state['scripts'].keys()),
        frame_props={'s': 36, 'linewidth': 2, 'edgecolors': theme_color},
        check_props={'linewidth': 2.8, 'color': theme_color},
    )
    ax_check_btn.text(0.2, 0.9, 'Scripts', transform=fig_setvar.transFigure, ha='center', va='center',
                      backgroundcolor=bg_color)
    fig_setvar.canvas.draw_idle()


# widgets for changing parameters

style_slider = {
    'valinit': 0,
    'initcolor': None,
    'color': theme_color,
    'track_color': silver,
    'handle_style': {'size': 14},
}

T0_slider = Slider(
    ax=ax_T0_slider,
    label='Feed Temperature (K)',
    valmin=273,
    valmax=600,
    valfmt='%.0f',
    **style_slider
)
T0_slider.on_changed(update_T0)

Tc_slider = Slider(
    ax=ax_Tc_slider,
    label='Coolant Temperature (K)',
    valmin=273,
    valmax=600,
    valfmt='%.0f',
    **style_slider
)
Tc_slider.on_changed(update_Tc)

UA_slider = Slider(
    ax=ax_UA_slider,
    label='UA\n' + r'$\mathrm{(J/s/K)}$',
    valmin=0,
    valmax=1e4,
    orientation='vertical',
    **style_slider
)
UA_slider.on_changed(update_UA)

v_slider = Slider(
    ax=ax_v_slider,
    label='VFR\n' + r'$\mathrm{(m^3/s)}$',
    valmin=0,
    valmax=1,
    orientation='vertical',
    valfmt='%.3f',
    **style_slider
)
v_slider.on_changed(update_v)

radio_btn = RadioButtons(
    ax=ax_radio_btn,
    labels=reaction_list,
    activecolor=theme_color,
    radio_props={'edgecolor': silver, 's': 36},
)
radio_btn.on_clicked(select_reaction)

tolerance_input = TextBox(
    ax=ax_tolerance_input,
    label='Tolerance: '.ljust(40),
    initial=str(state['tol']),
    textalignment='center',
    color=bg_color,
    hovercolor=silver,
)
tolerance_input.on_submit(lambda value: state.update({'tol': eval(value)}))

speed_input = TextBox(
    ax=ax_speed_input,
    label='Seconds Per Frame: '.ljust(33),
    initial=str(state['t_max']),
    textalignment='center',
    color=bg_color,
    hovercolor=silver,
)
speed_input.on_submit(lambda value: state.update({'t_max': eval(value)}))

reset_btn = Button(
    ax=ax_reset_btn,
    label='Reset',
    color=bg_color,
    hovercolor=silver,
)
reset_btn.on_clicked(lambda event: reset())

check_btn = CheckButtons(
    ax=ax_check_btn,
    labels=[],
)

window_closed = False


def on_close():
    global window_closed
    window_closed = True


fig.canvas.mpl_connect('close_event', lambda event: on_close())


def run_scripts(t):
    globs = {
        "t": t,
        "sin": math.sin,
        "cos": math.cos,
        "pi": math.pi,
    }
    scripts = {'T0': [], 'Tc': [], 'v': [], 'UA': []}
    for is_active, (varname, script) in zip(check_btn.get_status(), state['scripts'].values()):
        if is_active:
            scripts[varname].append(script)

    for (varname, script), slider in zip(scripts.items(), (T0_slider, Tc_slider, v_slider, UA_slider)):
        if len(script) > 0:
            slider.active = False
        else:
            slider.active = True
        for expr in script:
            globs.update({"T0": reactor.T0, "Tc": reactor.Tc, "v": reactor.v, "UA": reactor.UA})
            try:
                setattr(reactor, varname, max(eval(expr, globs), 0))
            except Exception as err:
                prompt.set_text(f'Error in executing script: {err}')
                state['prompt_start'] = time.time()
        if (val := getattr(reactor, varname)) != slider.val:
            slider.set_val(val)


def update_data():
    t, C, T, Q = state['t'], state['C'], state['T'], state['Q']
    run_scripts(t[-1])
    C_new, T_new, t_new = reactor.step(C[-1], T[-1], state['t_max'], state['tol'], method='RK45')
    if np.any(np.isnan(C_new[-1])):
        reset()
        prompt.set_text('Warning: Program has been reset because nan values were encountered. '
                        'Try to reduce tolerance or slow down your mouse.')
        state['prompt_start'] = time.time()
        return
    Q_new = reactor.Q_gen(C_new, T_new)

    t = np.append(t, t[-1] + t_new)
    C = np.append(C, C_new, axis=0)
    T = np.append(T, T_new)
    Q = np.append(Q, Q_new)

    idx = np.argmax(t[-1] - t < state['t_max'] * 100)
    t = np.delete(t, np.s_[:idx])
    C = np.delete(C, np.s_[:idx], axis=0)
    T = np.delete(T, np.s_[:idx])
    Q = np.delete(Q, np.s_[:idx])

    state.update({'t': t, 'C': C, 'T': T, 'Q': Q})

    [line.set_data(t, Ci) for line, Ci in zip(ln_Ct, C.T)]
    ln_Tt.set_data(t, T)
    ln_QT.set_data(T, Q)
    dot_Q.set_data([T[-1], T[-1]], [-Q[-1], Q[-1]])

    ax_Tt.set_xlim(t[0], t[-1])
    ax_Ct.set_xlim(t[0], t[-1])
    axC_ymax = max(np.max(reactor.C0), np.max(C))
    ax_Ct.set_ylim(-axC_ymax * 0.05, axC_ymax * 1.05)
    T_max = max([np.nanmax(T), state['T_lim'][1]])
    ax_Tt.set_ylim(top=T_max + (T_max - 270) * 0.05)
    ax_QT.set_xlim(right=T_max)
    if state['Q_lim'][0] != state['Q_lim'][1]:
        axQ_ymax = max(state['Q_lim'][1], min(np.max(Q), 1.5 * state['Q_lim'][1]))
        ax_QT.set_ylim(state['Q_lim'][0], axQ_ymax)

    axT_legend.texts[0].set_text(f'Reactor Temperature: {T[-1]:.2f} K')
    if prompt.get_text() != '' and time.time() - state['prompt_start'] > 10.:
        prompt.set_text('')

    fig.canvas.draw_idle()
    fig.canvas.flush_events()


def initialize():
    import random
    global reactor, ln_Ct, axC_legend

    reactor = CSTR(state['params'][:-1], **state['params'][-1])

    colors = iter(random.sample(ln_C_colors, k=len(reactor.components)))
    [line.remove() for line in ln_Ct]
    ln_Ct = [ax_Ct.plot(0, label=comp, color=next(colors), linewidth=1.2)[0] for comp in reactor.components]
    if axC_legend is not None:
        axC_legend.remove()
    axC_legend = ax_Ct.legend(loc='lower left')

    T0_slider.valinit = reactor.T0
    Tc_slider.valinit = reactor.Tc
    UA_slider.valinit = reactor.UA
    UA_max = 1e4 if reactor.UA < 1e4 else reactor.UA * 2
    UA_slider.valmax = UA_max
    UA_slider.ax.set_ylim(0, UA_max)
    v_slider.valinit = reactor.v
    v_slider.valmax = min(10 * reactor.v, 2 * reactor.VR)
    v_slider.ax.set_ylim(0, v_slider.valmax)

    T_max = reactor.solve_steady_state()
    draw_contour()

    state['T_lim'][1] = T_max
    ax_Tt.set_ylim(*state['T_lim'])
    ax_QT.set_xlim(*state['T_lim'])
    for slider in T0_slider, Tc_slider:
        slider.valmax = T_max
        slider.ax.set_xlim(right=T_max)
    T = np.linspace(*state['T_lim'])
    Q_gen = reactor.Q_gen_steady(T, reactor.v / reactor.VR)
    ln_Q_gen.set_data(T, Q_gen)
    ln_Q_rem.set_data(T, reactor.Q_rem(T))
    state['Q_lim'] = [np.min(Q_gen), np.max(Q_gen) * 1.5]
    reset()


if __name__ == '__main__':
    fig_setvar.canvas.manager.set_window_title('User Defined Variable Behavior')
    fig.canvas.manager.set_window_title('Initializing...')
    plt.show(block=False)
    fig.canvas.draw_idle()
    fig.canvas.flush_events()
    radio_btn.set_active(0)
    fig.canvas.manager.set_window_title(title)

    while not window_closed:
        update_data()
