################################################################################
######### © 2020 ETH Zurich, Institute of Geophysics, Daniel T. Birdsell #######
################################################################################

def COP_fn_mdot(L, transmissivity, mdot,dT):
    """
    Calculates COP as a function of mass flow rate.

    L is the well spacing [m]
    transmissivity is the transmissivity [m^3]
    mdot is the mass flow rate [kg/s]
    dT is the difference in temperature between T_CV and T_DH, which depends
        on depth
    """
    import global_vars as GV
    import numpy as np

    # import pudb; pudb.set_trace()
    COP = (GV.rho_w**2*GV.Cpw*dT*np.pi*transmissivity) / \
          (mdot*GV.viscosity*np.log(L/GV.D))
    return COP

def m_max_ResVol(L, b):
    """
    Defines the maximum mass flow rate [kg/s] based on the amount of heat the
    reservoir volume can hold. The reservoir volume is approximated as L^2*b.
    This is Constraint I in Birdsell paper.

    L is the well spacing [m].
    b is the aquifer thickness [m].
    """
    import global_vars as GV

    mdot = (GV.rho_w*GV.Cpw*GV.porosity + GV.rho_r*GV.Cpr*(1.0-GV.porosity)) * \
           (GV.alphaI**2 * L**2 * b) / (GV.Cpw*GV.t_inj);
    return mdot

def COP_HF(L, dT, reservoir_depth = None):
    """
    Defines the maximum mass flow rate [kg/s] to avoid hydraulic fracuting. This
    is Constraint II in Birdsell paper.

    L is the distance between the wells [m].
    reservoir_depth is the depth of the reservoir [m]. Note - if it is changing
        within the script, it should be provided a value. Otherwise it will
        default to the value in the global_vars.py
    dT is the difference in temperature between T_CV and T_DH, which depends
        on depth
    """
    import global_vars as GV

    if reservoir_depth is None:
        print('Warning - you are using the default reservoir_depth')
        reservoir_depth = GV.reservoir_depth;
    # import pudb; pudb.set_trace()
    COP_II = (GV.rho_w*GV.Cpw*dT)/2.0/(GV.alphaII*GV.rho_r - GV.rho_w)/ \
               GV.gravity/reservoir_depth
    return COP_II

def mdot_minLCOH(well_cost, transmissivity, Lstar):
    """
    Calculates the flow rate that gives the minimum LCOH, from d(LCOH)/d(mdot)

    well_cost is the construction cost of one well [$]. The capital cost are 4
        times the cost of a single well because there are two wells and the
        capital cost is assumed to double.
    transmissivity is the reservoir transmissivity [m3]
    Lstar is the optimal well spacing
    """
    import global_vars as GV
    import numpy as np

    capital_cost = 4.0*well_cost
    cap_cost = capital_cost
    CRF = (GV.r*(1.0+GV.r)**GV.lifetime) / ((1.0+GV.r)**GV.lifetime-1.0)
    mdot_max_val = np.sqrt(capital_cost*CRF*GV.rho_w**2 *np.pi*transmissivity/ \
                           (2.0*GV.joule_to_kWh * GV.dollars_per_kWhth* \
                            GV.t_inj*GV.viscosity*np.log(Lstar/GV.D)))
    return mdot_max_val

def m_max_COP(L, b, k, dT):
    """
    Defines the maximum mass flow rate (kg/s) based on the COP>1 constraint.
    This is not

    L is the distance between the wells.
    b is the aquifer thickness
    k is the permeability
    dT is the difference in temperature between T_CV and T_DH, which depends
        on depth
    """
    import global_vars as GV
    import numpy as np

    # import pudb; pudb.set_trace()
    mdot = (GV.rho_w*GV.Cpw*dT)*(GV.rho_w*k*b*np.pi)/(2.0*GV.viscosity * \
            np.log(L/GV.D))
    return mdot

def m_max_HF(L , b , k, reservoir_depth = None):
    """
    defines the maximum mass flow rate (kg/s) to avoid hydraulic fracuting. this
    assumes that pore pressure should be less than 80% of lithostatic stress.

    L is half the distance between the wells.
    b is the aquifer thickness
    k is the permeability
    """
    import global_vars as GV
    import numpy as np

    if reservoir_depth is None:
        reservoir_depth = GV.reservoir_depth
        print('Warning: m_max_HF used GV definition of reservoir_volume')

    max_dP_allowed = ((GV.alphaII*GV.rho_r) - GV.rho_w) * GV.gravity*reservoir_depth; #Pa

    mdot =(((GV.alphaII*GV.rho_r) - GV.rho_w) * GV.gravity*reservoir_depth) * \
               (2.0*np.pi*GV.rho_w*k*b) / (GV.viscosity*np.log(L/GV.D))
    return mdot

def thermal_radius(mdot , b, t_inj = None ):
    """
    Calculates the thermal radius

    mdot is the flow rate [kg/s]
    t_inj is the duration of the injection [s]
    b is the reservoir thickness
    """
    import global_vars as GV
    import numpy as np

    if t_inj is None:
        t_inj = GV.t_inj
        print("Warning t_inj in RF.thermal_radius takes default")
    R_th = np.sqrt(GV.Cpw*mdot*t_inj/((GV.Cpw*GV.porosity*GV.rho_w + \
                  GV.Cpr*(1.0-GV.porosity)*GV.rho_r)*np.pi*b))
    return R_th

def dP_inj_fn(mdot, transmissivity , L ):
    """
    Calculates the change in pressure for a single well from Darcy's equation.
    See Schaetlze (1980).

    mdot is flow rate [kg/s]
    transmissivity is reservoir transmissivity [m3]
    L is the well spacing [m]
    """
    import global_vars as GV
    import numpy as np

    dP_inj = mdot*GV.viscosity*np.log(L/GV.D)/(2.0*np.pi*GV.rho_w*transmissivity)
    return dP_inj


def heat_loss_fn(depth, Rth, b, T_WH = None, T_DH = None, T_surf = None,
               geothermal_gradient = None,Lstar = None):
    """
    depth = reservoir depth [m]
    Rth = thermal radius [m]
    b = reservoir thickness [m]
    T_WH is the waste heat temperature
    T_DH is the district heating rejection temperature
    T_surf is the average surface temperature
    geothermal_gradient [C/m] is the geothermal gradient
    Lstar is well spacing (only used if heat loss area is based on well spacing).
    """
    import global_vars as GV
    import numpy as np

    if T_WH is None:
        T_WH = GV.T_WH
    if T_DH is None:
        T_DH = GV.T_DH
    if T_surf is None:
        T_surf = GV.T_surf
    if geothermal_gradient is None:
        geothermal_gradient = GV.geothermal_gradient
    t_rest = GV.t_inj #assumes that the resting period is 1/4 year.
    T_G = T_surf + geothermal_gradient*depth

    thermal_radius = True
    if thermal_radius:
        interfacial_area = 2.0*np.pi*Rth**2
        reservoir_volume = np.pi*Rth**2*b
    else: ## based on reservoir volume, probably wont use
        interfacial_area = 2.0*(alphaI*Lstar)**2
        reservoir_volume = Lstar**2*b
    C1 = (GV.rho_w*GV.Cpw*GV.porosity + GV.rho_r*GV.Cpr*(1.0-GV.porosity))*reservoir_volume #J/C
    C2 = GV.lambda_eff*interfacial_area/GV.DeltaL
    Tcv = (T_WH-T_G)*np.exp(-C2/C1*t_rest)+T_G #temperature after rest
    heat_loss = (T_WH - Tcv)*C1 #joules
    final_recovery = (Tcv - T_DH)*C1 #joules
    efficiency = (Tcv-T_DH)/(T_WH-T_DH) # Eout/Ein. Eout = mCpf*(Tcv-TDH). Ein = m*Cpf*(Twh-Tdh)
    return [Tcv, T_G, heat_loss, final_recovery, efficiency]
