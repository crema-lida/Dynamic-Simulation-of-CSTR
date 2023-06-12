# Dynamic Simulation of CSTR
This app allows you to control different operating parameters (e.g., feed temperature and volume flow rate) and study the dynamic behavior of a CSTR. Furthermore, you are able to define your own power-law reactions and scripts to change input parameters.

## Transient and steady-state curves
<img width="720" alt="cascade" src="https://github.com/real-Crema/Dynamic-Simulation-of-CSTR/assets/100750226/f9548e07-1e59-4ddb-8fa4-465019202390">

The change of reactor temperature and concentrations of components with time is obtained by solving a system of ordinary differential equations:

$$\frac{(c_{j0}-c_j)v}{V_R}+r_j=\frac{dc_j}{dt}$$

$$\sum\frac{\Delta H_i\ r_i}{\rho C_P}+\frac{v(T_0-T)}{V_R}+\frac{UA(T_c-T)}{V_R\rho C_p}=\frac{dT}{dt}$$

where $r_j$ in the material balance equation denotes the net generation rate of component $j$, and $r_i$ in the heat balance equation is solely the consumption rate of component $i$ in each reaction.

$$r_i=-k_{0i}\ \mathrm{exp}(-\frac{E_a}{RT})\prod c_k^n$$

These equations are solved using explicit Runge-Kutta method of order 5(4). The error is controlled assuming accuracy of the fourth-order method, but steps are taken using the fifth-order accurate formula.

The steady-state heat generation/removal rate ($Q_{gen},\ Q_{rem}$) as a function of reactor temperature is obtained by solving steady-state material balance equations when $dc_j/dt=0$.

## $T_{reactor}/T_{coolant}$ Behavior changes as a function of space velocity ($v/V_R$)

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

Reactions can be defined by a `.json` file under `reactions` directory. An example file for the first-order reaction $A\longrightarrow P$ is shown below.

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

We can define multiple reactions in one file, so that we can simulate a complex system composed of many simple power-law reactions. Take the van de Vusse reaction for example:

$$A\stackrel{k_1}{\longrightarrow} B \stackrel{k_2}{\longrightarrow} C$$

$$2A\stackrel{k_3}{\longrightarrow} D$$

The `.json` file for the van de Vusse reaction looks like this:

```
[
    {
        "A": [-1, 1],
        "B": [1, 0],
        "k0": 3.575e8,
        "Ea": 81.1305,
        "dH": 4.2e6
    },
    {
        "B": [-1, 1],
        "C": [1, 0],
        "k0": 3.575e8,
        "Ea": 81.1305,
        "dH": -11e6
    },
    {
        "A": [-2, 2],
        "D": [1, 0],
        "k0": 2.512e6,
        "Ea": 71.16784,
        "dH": -41.85e6
    },
    {
        "C0": {
            "A": 5.1,
            "B": 0,
            "C": 0,
            "D": 0
        },
        "VR": 0.01,
        "v": 1.667e-3,
        "T0": 387.05,
        "rho": 934.2,
        "Cp": 3010,
        "Tc": 350,
        "UA": 240.8,
        "initial": {
            "C": [2.2291, 1.0417, 0.9140, 0.9152],
            "T": 352.741
        }
    }
]
```

We use 3 `{}` to define those reactions and then use the last `{}` to specify other necessary operating parameters. Be aware that enthalpy of reaction ( `"dH"` ) and reaction rate is calculated based on the first component written in each `{}`. In this case, we also set initial values for concentrations and reactor temperature by specifying `"initial"` (optional). If specified, values in `initial["C"]` will be assigned to each component at t=0; if not specified, the initial values of `"C"` and `"T"` will be identical to `"C0"` and `"T0"`.

## Define custom scripts 

The control over operating variables ($T_0,\ T_c,\ v,\ UA$) can be automated via user-defined scripts in `scripts/scripts.json`. Each entry contains the script name and a piece of python code.

```
{
  "T0 increases linearly": ["T0", "T0 + 0.5"],
  "T0 decreases linearly": ["T0", "T0 - 0.5"],
  "T0 fluctuates in a sine wave": ["T0", "350 + 20 * sin( 0.02 * pi * t)"],
  "Tc increases linearly": ["Tc", "Tc + 0.5"],
  "Tc decreases linearly": ["Tc", "Tc - 0.5"],
  "Tc fluctuates in a sine wave": ["Tc", "320 + 20 * sin( 0.02 * pi * t)"],
  "v increases linearly": ["v", "v + 0.0005"],
  "v decreases linearly": ["v", "v - 0.0005"],
  "UA increases linearly": ["UA", "UA + 20"],
  "UA decreases linearly": ["UA", "UA - 20"]
}
```

Each brackets `[]` consists of two strings: the variable to control and the value to be assigned to it. All scripts are executed once in each frame. That is to say, if you set Seconds Per Frame to 10, and select the first script, the program will do `T0 = T0 + 0.5` at each frame, so that `T0` will increase 0.5 K every 10 seconds.

Once a script is created or modified, save the file and click the Reset button, and it will automatically appear on this window:

<img width="427" alt="control-variables" src="https://github.com/real-Crema/Dynamic-Simulation-of-CSTR/assets/100750226/a7b4fa02-f0cd-4748-b2ed-409079536a61">

