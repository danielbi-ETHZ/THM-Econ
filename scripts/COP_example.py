################################################################################
######### © 2020 ETH Zurich, Institute of Geophysics, Daniel T. Birdsell #######
################################################################################

### Generate Figure 2
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..') ## point to python functions in main directory for import
from global_vars import *
from plotting_vars import *
import reservoir_functions as RF
import economic_functions as EF

################################################################################
######## PART 1: Initilize functions, parameters, figure #######################
################################################################################
### Define the well spacing functions within each script so they can be accessed
###   by fsolve.
def L_star_fun(L_star):
    """
    This calcualtes the zero of the optimal well spacing, L_star.
    """
    import numpy as np
    return L_star**2.0*np.log(L_star/D) - \
           2.0*np.pi*rho_w*permeability_/viscosity * \
           ((alphaII*rho_r-rho_w)*gravity*reservoir_depth) * \
           Cpw*t_inj/((rho_w*Cpw*porosity+rho_r*Cpr*(1.0-porosity)))

def L_star_fun2(L_star):
    """
    This calcualtes the optimal well spacing, L_star, if the reservoir
    constraints imply a flow rate that is higher than the flow rate that
    would minimize the LCOH.

    define capital_costs and CRF externally
    """
    import numpy as np

    return (capital_cost_internal*CRF*rho_w*rho_w * np.pi * permeability_*b_ /
           (2.0*dollars_per_kWhth*joule_to_kWh*t_inj*viscosity *
            np.log(L_star/D)))**0.5 - \
           (rho_w*Cpw*porosity + rho_r*Cpr*(1.0-porosity)) * \
           (L_star**2 * b_)/ (Cpw*t_inj)

### Select reservoir_depth and associated dT.
reservoir_depth = 575;
dT = 35; # dT is T_CV - T_DH. Requires finding T_CV a-priori with heat_loss_fn

### Set Parameters - DO NOT CHANGE
permeability_ = 1.e-13;
b_ = 20. # reservoir thickness
transmissivity = permeability_*b_; #transmissivity is permeability * thickness

### Deine plotting colors (from plotting_vars.py)
colorI =  Yellow
colorII = Blue
colorIII = Red
COP2_color = (0.8,0.8,0.8)

### Plotting parameters
transmissivity_vec = [transmissivity/100.,transmissivity/10.,transmissivity,
                      transmissivity*10.,transmissivity*100.]
transmissivity_solid = transmissivity; # Solid line for transmissivity
xlim_min = -0.; xlim_max = 2.1
ylim_min = 0.; ylim_max = 4.5

### Initialize a figure
plt.figure()
plt.xlabel(r'log$_{10}$($\.{m}$ [kg/s])')
plt.ylabel(r'log$_{10}$(COP)')
plt.xlim([xlim_min,xlim_max])
plt.ylim([ylim_min,ylim_max])

################################################################################
######## PART 2: Calculations ##################################################
################################################################################
### Find optimal well spacing for base-case scenario.
###   This involves evaluating the reservoir and economic constraints.
###   In this case, the flow rate and well spacing are essentially the same,
###   but this is only because the base case was chosen so all three Constraints
###   are equal.
well_costs_LD = EF.Well_cost(adjust = True,reservoir_depth = reservoir_depth)
capital_cost_internal = 4.0*well_costs_LD
import scipy.optimize
Lstar_res = scipy.optimize.fsolve(L_star_fun,100) #L_res
Lstar_LCOH = scipy.optimize.fsolve(L_star_fun2,100) #L_econ
mdot_res = RF.m_max_ResVol(L = Lstar_res, b = b_)
mdot_LCOH = RF.mdot_minLCOH(well_cost = well_costs_LD,
                            transmissivity = transmissivity, Lstar = Lstar_LCOH)
if mdot_res < mdot_LCOH:
    Lstar = Lstar_res
else:
    Lstar = Lstar_LCOH

### Set the well spacing for all transmissivity to Lstar from base case
L = Lstar;

### Calculate and plot the COP as a function of flow rate.
mdot = np.linspace(.1,200,51)
for Tcount in transmissivity_vec:
    T = Tcount
    COP = np.zeros(len(mdot))
    for i in range(0,len(COP)):
        COP[i] = RF.COP_fn_mdot(L = Lstar, transmissivity = T, mdot = mdot[i],
                                dT = dT)
    if T == transmissivity_solid:
        plt.plot(np.log10(mdot),np.log10(COP),'k-',linewidth=normal_width)
    else:
        plt.plot(np.log10(mdot),np.log10(COP),'k:',linewidth=normal_width)


### Calculate and plot Constraint 1.
bvec = [10, 20, 40] ## reservoir thicknesses.
b_solid = b_; # this is the solid line in the plot
for bcount in bvec:
    b = bcount
    m_dot_max = RF.m_max_ResVol(L = Lstar, b = b)
    if b == b_solid:
        plt.plot([np.log10(m_dot_max),np.log10(m_dot_max)],[-2,10],'-',
                  linewidth=base_case_width,color = colorI)
    else:
        plt.plot([np.log10(m_dot_max),np.log10(m_dot_max)],[-2,10],'--',
                  linewidth=normal_width,color=colorI)

### Calculate and plot Constraint 2.
COP_min = RF.COP_HF(L=L, dT = dT, reservoir_depth = reservoir_depth)
plt.plot([np.log10(.1),np.log10(200)],[np.log10(COP_min),np.log10(COP_min)],
          color=colorII)

### Calculate and plot Constraint 3.
COP_LCOH = RF.COP_fn_mdot(L = Lstar, transmissivity= transmissivity_solid,
                          mdot = mdot_LCOH, dT = dT)
plt.plot(np.log10(mdot_LCOH),np.log10(COP_LCOH),marker='*',markersize = '10',
         color=colorIII)

### Plot COP = 2 as a boundary to avoid (grey)
COP2 = 2.0
plt.plot([xlim_min,xlim_max],[np.log10(COP2),np.log10(COP2)],color = COP2_color,
          linewidth=normal_width)

### Plot mdot = 3.3 kg/s (grey)
mdot_line = 3.3
plt.plot([np.log10(mdot_line),np.log10(mdot_line)],[ylim_min,ylim_max],
          color = COP2_color,linewidth=normal_width)

plt.savefig(FIGURE_DIRECTORY+'/COP_example.png')
plt.show()
