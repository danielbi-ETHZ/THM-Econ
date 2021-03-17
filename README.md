# THM$ Code for HT-ATES

A thermo-hydro-mechanical-economic (THM$) approach to understand
high-temperature aquifer thermal energy storage (HT-ATES).
While the approach is analytical, 
it uses Python 3.7.7. to solve implicit equations, interpolate, pass variables,
plot, etc.
This approach was first laid out in Birdsell et al. (2021).

## Layout of the directory

The root directory has a number of functions that contain the main analytical
reservoir and economic relationships and commonly-used, global variables.
These should be imported with:

```
import economic_functions as EF
import reservoir_functions as RF
from global_vars import *
```

The *./scripts* directory has example problem(s) that show how these functions
are called.


## Variables

This code relies extensively on global variables, so it is very important to
keep track of names. The global variables are set to be equal to the values used
in Table 1 of Birdsell et al. (2020), with the exceptions of: (a) the
permeability,
(b) the aquifer thickness, (c) the depth, and (d) the faulting regime (i.e. alphaII) which can all change within
individual scripts, based on the plots that are generated within the scripts
functions. The global variables are stored in ```global_vars.py```

## Note on the calculation of optimal well spacing.
Since the optimal well spacing is defined by one of two intrinsic functions, it is solved
with *scipy.optimize.fsolve()*. L_star_fun gives the well spacing according to the reservoir constraints (i.e. L_res), and L_star_fun2 gives the well spacing according to the economic constraints (i.e. L_econ).
Both must be calculated, since we do not know which regime we are in a priori. The optimal well spacing is the smaller of the well spacing values from the reservoir and economic constraints.
From what Daniel Birdsell can tell,
the *fsolve()* function cannot easily use an imported function, or a function
with multiple arguments. Therefore, the
*L_star_fun* and the *L_star_fun2* are defined within each script that they are
used. This involves some copy and pasting and some values that must be calculated
before *fsolve()* is used. Therefore it is important to define and/or calculate
*capital_cost_internal*, *permeability_*, and *b_*. An example is provided below:
```
permeability_ = 1.e-13;
b_ = 20.
reservoir_depth = 575.
well_costs_LD = EF.Well_cost(adjust = True,reservoir_depth = reservoir_depth)
capital_cost_internal = 4.0*well_costs_LD[i]
Lstars_res = scipy.optimize.fsolve(L_star_fun,100)
Lstars_LCOH = scipy.optimize.fsolve(L_star_fun2,100)
mdot_res = RF.m_max_ResVol(L = Lstars_res, b = b_)
mdot_LCOH = RF.mdot_minLCOH(well_cost = well_costs_LD,
                            transmissivity = T, Lstar = Lstars_LCOH)
if mdot_res<mdot_LCOH:
    mdot = mdot_res
    Lstars = Lstars_res
else:
    mdot = mdot_LCOH
    Lstars = Lstars_LCOH
```

## Three Constraints

The THM$ approach balances three constraints:

1. Constraint 1 defines a maximum flow rate so that the reservoir's thermal
capacity is not over-utilized. It is hydro-thermal in nature and is calculated
with: ```RF.m_max_ResVol(L, b)```.
2. Constraint 2 defines a maximum flow rate so that the reservoir pressure does
not lead to hydraulic fracturing (HF). It is hydro-mechanical in nature and is
calculated with: ```RF.m_max_HF(L , b , k, reservoir_depth = None)```.
3. Constraint 3 defines the maximum flow rate such that the levelized cost of
heat (LCOH) is minimized. It is realted to economics and hydro-thermal reservoir
 engineering to a lesser degree.
 It is calculated with: ```RF.mdot_minLCOH(well_cost, transmissivity, Lstar)```

## Contact

To get in contact: https://geg.ethz.ch/daniel-birdsell/


## References

Birdsell, D.T., Adams, B.M., and Saar, M.O., 2021. Minimum Transmissivity and Optimal Well Spacing and Flow Rate for High-Temperature Aquifer Thermal Energy Storage. _Applied Energy_.

<!-- To speed up computations, it is possible to pre-calculate Lstar as a function
of reservoir_depth and permeability. In this way, the optimal well spacing
does not need to re-calculated many times. This can be done by running
*pre_calc_lstars.py*, which outputs a python file with variables which can
be imported. -->
