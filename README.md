# Dynamic Simulation of CSTR
This app allows you to control different operating parameters (e.g., feed temperature and volume flow rate) and study the dynamic behavior of a CSTR. Furthermore, you are able to define other power-law reactions and scripts to change input parameters.

## Main Window
<img width="818" alt="main_window" src="https://github.com/real-Crema/Dynamic-Simulation-of-CSTR/assets/100750226/baea4c3c-66ed-4d5a-a533-f32ab3191ced">

## How To Read Contour Plot



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
