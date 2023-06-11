# Dynamic Simulation of CSTR
This app allows you to control different operating parameters (e.g., feed temperature and volume flow rate) and study the dynamic behavior of a CSTR. Furthermore, you are able to define your own power-law reactions and scripts to change input parameters.

## Transient and Steady-State Curves
<img width="720" alt="main_window" src="https://github.com/real-Crema/Dynamic-Simulation-of-CSTR/assets/100750226/baea4c3c-66ed-4d5a-a533-f32ab3191ced">

The change of reactor temperature and concentrations of components with time is obtained by solving a system of ordinary differential equations:

$$\frac{(c_{j0}-c_j)v}{V_R}+r_j=\frac{dc_j}{dt}$$

$$\sum\frac{\Delta H_i\ r_i}{\rho C_P}+\frac{v(T_0-T)}{V_R}+\frac{UA(T_c-T)}{V_R\rho C_p}=\frac{dT}{dt}$$

where $r_j$ in the material balance equation denotes the net generation rate of component $j$, and $r_i$ in the heat balance equation is solely the consumption rate of component $i$ in each reaction.

$$r_i=k_{0i}\ \mathrm{exp}(\frac{E_a}{RT})\prod c_k^n$$

These equations are solved using explicit Runge-Kutta method of order 5(4). The error is controlled assuming accuracy of the fourth-order method, but steps are taken using the fifth-order accurate formula.

The steady-state heat generation/removal rate ($Q_{gen},\ Q_{rem}$) as a function of reactor temperature is obtained by solving steady-state material balance equations when $dc_j/dt=0$.

## $T_{reactor}/T_{coolant}$ Behavior Changes as a Function of Space Velocity ($v/V_R$)

The steady-state heat balance equation ($dT/dt=0$) is

$$\sum\Delta H_i\ r_i V_R=v \rho C_p(T-T_0)+UA(T-T_c)$$

which simply means

$$Q_{gen}=Q_{rem}$$

We can rearrange this equation to find the required coolant temperature $T_c$ for any given steady-state reactor temperature $T$ and volume flow rate $v$.

$$T_c=T+\frac{v \rho C_p (T-T_0)-Q_{gen}}{UA}$$

The contour plot below demonstrates the relationship between $T_{reactor}$ (y values) and $T_{coolant}$ (z values) at different space velocities (x values) in a cascade reaction $A\longrightarrow P \longrightarrow S$.

<img width="398" alt="TvT_contour" src="https://github.com/real-Crema/Dynamic-Simulation-of-CSTR/assets/100750226/de0d5d99-cab9-41b0-b51c-f2d21c3f41ee">

If we draw a straight line at a certain space velocity, like the yellow one, we find 5 intersections at $T_{coolant}$=300 K, which correspond to 5 possible steady-states. Similarly, the red line has 3 intersections at $T_{coolant}$=400 K. As we decrease $v/V_R$, we moves from multiple steady-state to monotonic behavior.
## Define custom reactions

```
[
    {
        "A": [-1, 1],  // "reactant": [stoichiometry, exponent]
        "P": [1, 0],   // "product": [stoichiometry, exponent]
        "k0": 1e13,    // 1/s
        "Ea": 100,     // kJ/mol
        "dH": -2e7     // J/kmol
    },
    {
        "C0": {
            "A": 5,    // kmol/m3
            "P": 0     // kmol/m3
        },
        "VR": 10,      // m3
        "v": 1e-2,     // m3/s
        "T0": 300,     // K
        "rho": 850,    // kg/m3
        "Cp": 2200,    // J/(kg*K)
        "Tc": 300,     // K
        "UA": 1000     // J/(s*K)
    }
]
```
