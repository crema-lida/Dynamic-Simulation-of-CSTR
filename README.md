## Main Window

<img src="C:\Users\Crema\AppData\Roaming\Typora\typora-user-images\image-20230611215103123.png" alt="image-20230611215103123" style="zoom: 50%;" />

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