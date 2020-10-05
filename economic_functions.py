################################################################################
########### © 2020 ETH Zurich, Daniel T. Birdsell ##############################
################################################################################

def Well_cost(reservoir_depth, adjust = None):
    """
    Gives the cost of drilling and construction of one well in 2010 USD.
    Taken from GETEM, for large-diameter (31.1 cm) well.
    See also Ben Adams SI Tables S1 and S15.

    reservoir_depth is in [m]
    adjust = True accounts for inflation in O&G industry, adjusting from USD in
        year 2010 to USD in year 2019.
    """
    import global_vars as GV

    if adjust is None:
        adjust = True
    if adjust:
        PPIOG_2010 = 2.123; #
        PPIOG_2019 = 2.195; #
        PPI_adjustment = PPIOG_2019/PPIOG_2010
    else:
        PPI_adjustment = 1.0 # 2002 dollars

    reservoir_depth = reservoir_depth/0.3048 # convert meters to feet
    C_well = 0.033*reservoir_depth**2 + 350*reservoir_depth + 290000
    C_well = C_well*PPI_adjustment
    return C_well

def LCOH_fn(well_cost = None, Ed = None, mdot = None, dP_inj = None):
    """
    Calculates the LCOH based on the well cost as capital cost
    and the pumping cost for operating cost, and the thermal energy recovered.

    Results match simplified_LCOH_fn, but simplified_LCOH_fn is preferred
    because it matches the notation in Birdsell et al. paper.

    well_cost is the cost of one well [$]
    Ed is the annual heat recovered [kWh]
    mdot is flow rate [kg/s]
    dP_inj is the injection pressure

    returns:
    LCOH [$/kWh]
    cap_cost of project [$]
    pump_cost is total project pumping cost, to present day $
    benefit_USD - is the total project amount of revenue,
        adjusted to present-day dollars, if the heat has a value of
        dollars_per_kWhth (usually $0.10/kWhth)
    """
    import global_vars as GV
    import numpy as np
    capital_cost = 4.0*well_cost
    ### The work into pumping is 2*mdot*dP_total/rho_w.
    ### The 2 is because of loading and unloading. dP_total= 2*mdot
    joules_per_second_pumping = (2.0*mdot*dP_inj/GV.rho_w) ## This is the J/s during pumping
    total_pumping_time = 2.0*GV.t_inj #seconds of pumping
    joules = joules_per_second_pumping*total_pumping_time
    annual_pumping_cost = joules*GV.joule_to_kWh*GV.dollars_per_kWhth

    cap_cost = 0 #Calculates the value with discount rate. initialize 0
    pump_cost = 0 #Calculates total pumping cost with discount rate. initialize to 0
    benefit = 0 # initialize (kWhth)
    benefit_USD = 0 # initizlize (USD)
    ### Loop through the lifetime
    for t in range(0,GV.lifetime+1):
        if t == 0:
            cap_cost = capital_cost/(1.+GV.r)**t + cap_cost
        else:
            pump_cost = annual_pumping_cost/(1.+GV.r)**t + pump_cost
            benefit = Ed/(1.+GV.r)**t + benefit
            benefit_USD = (Ed/(1.+GV.r)**t)*GV.dollars_per_kWhth + benefit_USD
    LCOH = (cap_cost+pump_cost)/benefit #$/kWhth
    return [LCOH, cap_cost, pump_cost, benefit_USD]

def simplified_LCOH_fn(well_cost, Ed, mdot, dP_inj, return_all = False):
    """
    Calculates the LCOH based on the well cost as capital cost
    and the pumping cost for operating cost, and the thermal energy recovered.
    The difference between this and LCOH_fn, is this one uses CRF and Gives
    annualized values of capital cost, pumping cost, etc.

    simplified_LCOH_fn() is equivalent to LCOH_fn(), but is preferred because it
    matches the notation used in Birdsell et al. paper.

    well_cost is the cost of one well [$]
    Ed is the annual heat recovered [kWh]
    mdot is flow rate [kg/s]
    dP_inj is the injection pressure
    return_all - if it is False, it only returns the LCOH

    returns:
    LCOH [$/kWh]
    """
    import global_vars as GV
    import numpy as np

    capital_cost = 4.0*well_cost
    ### The work into pumping is 2*mdot*dP_total/rho_w.
    ### The 2 is because of loading and unloading. dP_total= 2*mdot
    joules_per_second_pumping = (2.0*mdot*dP_inj/GV.rho_w) ## This is the J/s during pumping
    total_pumping_time = 2*GV.t_inj #seconds of pumping
    joules = joules_per_second_pumping*total_pumping_time
    annual_pumping_cost = joules*GV.joule_to_kWh*GV.dollars_per_kWhth

    cap_cost = capital_cost
    CRF = (GV.r*(1+GV.r)**GV.lifetime) / ((1+GV.r)**GV.lifetime-1)
    LCOH = cap_cost*CRF/Ed + annual_pumping_cost/Ed
    if return_all:
        return [LCOH, cap_cost*CRF, annual_pumping_cost]
    else:
        return LCOH
