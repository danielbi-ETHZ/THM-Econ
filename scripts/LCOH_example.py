################################################################################
######### © 2020 ETH Zurich, Institute of Geophysics, Daniel T. Birdsell #######
################################################################################

### Generate Figure 6
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
import sys
sys.path.append('..') ## point to python functions in main directory for import
import economic_functions as EF
import reservoir_functions as RF
from global_vars import *
from plotting_vars import *
import scipy.optimize

################################################################################
######## PART 1: Initilize functions, parameters, figure #######################
################################################################################

### To change sign on depth for plotting purposes
import matplotlib.ticker as ticker

@ticker.FuncFormatter
def major_formatter(x, pos):
    label = str(int(-x)) if x < 0 else str(int(x))
    return label

### Define the well spacing functions within each script so they can be accessed
###   by fsolve.
def L_star_fun(L_star):
    """
    This calcualtes the zero of the optimal well spacing, L_star.
    """
    import numpy as np
    # import pudb; pudb.set_trace()
    return L_star**2.0*np.log(L_star/D) - \
           2.0*np.pi*rho_w*T/viscosity * \
           ((alphaII*rho_r-rho_w)*gravity*reservoir_depth) * \
           Cpw*t_inj/(b_*(rho_w*Cpw*porosity+rho_r*Cpr*(1.0-porosity)))

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

################################################################################
###### Part 2: Calculations ####################################################
################################################################################

### In the loops below, we calculate LCOH as a function of depth, thickness, and
###    faulting regime (Fig6a.png)
###    The values that are required to calculate LCOH (annualized capital cost,
###    annual operating cost, and heat recovered) are presented as a function of
###    depth in Fig6b.png
alpha_vec = [1.0,0.8] # Faulting regime
### First, loop over the stress state
for alphaII in alpha_vec:
    b_vec = [10,20,40] # Reservoir thickness
    permeability_ = 1.0e-13;
    ### Then, loop over the reservoir thickness
    for b_count in range(0,len(b_vec)):
        b_ = b_vec[b_count]
        T = permeability_*b_; #transmissivity = permeability * thickness
        depthmin = 50; depthmax = 2667; #depth range for plots
        depths = np.linspace(depthmin,depthmax,1001)
        Lstars = np.zeros(len(depths)) #the optimal well spacing
        Lstars_res = np.zeros(len(depths)) # spacing from reservoir constraints
        Lstars_LCOH = np.zeros(len(depths)) # spacing from economic constraints
        Rth = np.zeros(len(depths)) # thermal radius
        well_costs_LD = np.zeros(len(depths)) #Large diameter
        mdot = np.zeros(len(depths)) #optimal flow rate
        mdot_res = np.zeros(len(depths)) #flow rate from reservoir constraints
        mdot_LCOH = np.zeros(len(depths)) #flow rate from economic constraints
        dP_inj = np.zeros(len(depths)) #injection pressure
        Tcv = np.zeros(len(depths)) #control volume temp
        T_G = np.zeros(len(depths)) #geothermal temp
        heat_loss = np.zeros(len(depths)) #heat loss [J]
        heat_injection = np.zeros(len(depths)) #heat injected [J]
        final_recovery = np.zeros(len(depths)) #Joules recovered
        efficiency = np.zeros(len(depths)) #thermal efficiency
        Ed = np.zeros(len(depths)) #annual thermal energy produced [kWh]
        LCOH = np.zeros(len(depths)) #LCOH [$/kWh]
        annual_capital_cost = np.zeros(len(depths)) #annualized cap cost [$]
        annual_operating_cost = np.zeros(len(depths)) #annual operating cost [$]
        ### Finally, loop over the depth
        for i in range(0,len(depths)):
            ### Calculate optimal flow rate and well spacing at each depth
            reservoir_depth = depths[i]
            well_costs_LD[i] = EF.Well_cost(adjust = True,reservoir_depth = reservoir_depth)
            capital_cost_internal = 4.0*well_costs_LD[i]
            Lstars_res[i] = scipy.optimize.fsolve(L_star_fun,100)
            Lstars_LCOH[i] = scipy.optimize.fsolve(L_star_fun2,100)
            mdot_res[i] = RF.m_max_ResVol(L = Lstars_res[i], b = b_)
            mdot_LCOH[i] = RF.mdot_minLCOH(well_cost = well_costs_LD[i],
                                        transmissivity = T, Lstar = Lstars_LCOH[i])
            if mdot_res[i]<mdot_LCOH[i]:
                mdot[i] = mdot_res[i]
                Lstars[i] = Lstars_res[i]
            else:
                mdot[i] = mdot_LCOH[i]
                Lstars[i] = Lstars_LCOH[i]

            ### Change in pressure at injection well
            dP_inj[i] = RF.dP_inj_fn(mdot = mdot[i], transmissivity = permeability_*b_, L = Lstars[i])

            ### Thermal radius
            Rth[i] = RF.thermal_radius(mdot = mdot[i], t_inj = t_inj, b = b_)

            ### Calculate the heat injection, heat injection, and heat recovery
            [Tcv[i], T_G[i], heat_loss[i], final_recovery[i], efficiency[i]] = \
                            RF.heat_loss_fn(T_WH = T_WH, T_DH = T_DH, depth = depths[i],
                                         Rth = Rth[i], b = b_)
            heat_injection[i] = mdot[i]*Cpw*t_inj*(T_WH-T_DH)

            Ed[i] = final_recovery[i]*joule_to_kWh  #the yearly energy produced as a fun of depth [kWh]

            ### Calculate LCOH, operating cost, capital cost
            [LCOH[i], annual_capital_cost[i], annual_operating_cost[i]] = \
                              EF.simplified_LCOH_fn(well_cost = well_costs_LD[i], Ed = Ed[i],
                              mdot = mdot[i], dP_inj = dP_inj[i], return_all=True)


        ### Plot Figure 6(a) for each time through loop:
        plt.figure(10)
        if alphaII == 1:
            if b_count ==0:
                plt.plot(LCOH,-depths,'-',color = colors[-3], linewidth=normal_width)
            elif b_count == 1:
                plt.plot(LCOH,-depths,'-',color = colors[-2], linewidth=base_case_width)
            elif b_count == 2:
                plt.plot(LCOH,-depths, '-', color = colors[-1], linewidth=normal_width)
        elif alphaII == 0.8:
            if b_count ==0:
                plt.plot(LCOH,-depths,':',color = colors[-3],linewidth=normal_width)
            elif b_count == 1:
                plt.plot(LCOH,-depths,':',color = colors[-2], linewidth=normal_width)
            elif b_count == 2:
                plt.plot(LCOH,-depths, ':', color = colors[-1], linewidth=normal_width)
                leg = [
                       r'$b=$'+str(b_vec[0])+ r' m',
                       r'$b=$'+str(b_vec[1])+ r' m',
                       r'$b=$'+str(b_vec[2])+ r' m',
                       ]
                plt.legend(leg)
                plt.xlim([0,.15])
                plt.grid('on')
                plt.xlabel(r'LCOH [\$/kW$_{th}h$]',fontsize=label_font_size)
                plt.ylabel('Reservoir Depth [m]',fontsize=label_font_size)
                plt.xticks(fontsize=tick_font_size)
                plt.yticks(fontsize=tick_font_size)
                plt.text(0.0019,-125,'(a)',color='k',fontsize=normal_font_size)
                plt.savefig(FIGURE_DIRECTORY+'/LCOH_example_a.png')

        ### Plot Figure 6(b) using base case data
        if alphaII == 1.0 and b_count == 1:
            fig, ax = plt.subplots()
            ax.plot(Ed/1.e7,-depths,'k:',linewidth=normal_width)
            ax.plot(annual_operating_cost/1.e6,-depths,'k--',linewidth=base_case_width)
            ax.plot(annual_capital_cost/1.e6,-depths,'k-',linewidth=normal_width)
            ax.plot([0,3],[dstar_LCOH_min_alpha1[1],dstar_LCOH_min_alpha1[1]])
            ax.plot([0,3],[dstar_kmin_alpha1[1],dstar_kmin_alpha1[1]])
            ax.plot([0,3],[-570,-570])
            leg = [
                       r'$Q$',
                       r'$C_{op}$',
                       r'$C_{cap} \cdot CRF$',
                       r'$d^*_{LCOH,min}$',
                       r'$d^*_{kmin}$',
                       r'$d^*_{constraints}$'
                       ]
            plt.legend(leg,loc='upper right')
            plt.grid('on')
            ax.yaxis.set_major_formatter(major_formatter)
            plt.xlabel('Annual Costs [Million USD] or Heat Recovered [$10^{10}$ W$_{th}$-h]',fontsize=label_font_size);
            plt.ylabel('Reservoir Depth [m]',fontsize=label_font_size)
            plt.xticks(fontsize=tick_font_size)
            plt.yticks(fontsize=tick_font_size)
            plt.text(0.08,-2700,'(b)',color='k',fontsize=normal_font_size)
            plt.savefig(FIGURE_DIRECTORY+'/LCOH_example_b.png')

plt.show()
