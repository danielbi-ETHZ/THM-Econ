################################################################################
######### © 2020 ETH Zurich, Institute of Geophysics, Daniel T. Birdsell #######
################################################################################

### Define global parameters using units of J, kg, C, s, m
### The notation follows Birdsell et al. fairly closely.
### Do not change these variable names when calling global_vars from a script.
rho_w = 1000; #Water or other fluid density [kg/m3]
Cpw = 4186; #Heat capacity of the water or other fluid [J/kg/C]
T_WH = 90.; T_DH = 45.; #Waste heat and district heating return temperature [C]
geothermal_gradient = 30./1000. #Geothermal Gradient [C/m]
T_surf = 10; #Surface temperature [C]
D = 0.311-0.05; #The inner diameter of the well [m]. The large diameter from
                #GETEM is 31.1 cm or 0.311 m. 5 cm is subtracted for casing.
viscosity = 0.0005; #Fluid viscosity [Pa-s]
rho_r = 2500; #rock density [kg/m3]
Cpr = 850; #Rock heat capacity [J/kg/C]
porosity = 0.15; # Porosity [-]
t_inj = 365.25*24*60*60/4. #The length of each stage [s].
gravity = 9.81; #gravity [m2/s]
reservoir_depth = 575; #the depth of the reservoir [m]
alphaI = 1.0 #This scales the reservoir size. Constraint I in Birdsell et al.
alphaII = 1.0 # This is the percent of lithostatic that leads to HF.
lambda_r = 3.0 #Thermal conductivity of rock [W/m-C]
lambda_w = 0.6; #Thermal conductivity of water or other fluid [W/m-C]
lambda_eff = porosity*lambda_w+(1.0-porosity)*lambda_r #effective thermal
                                                       #conductivity [W/m-C]
DeltaL = 5.0 ; # Distance for temperature gradient in conductive heat loss [m]
lifetime = 25 ; #years of ATES operation [yr]
r = 0.03 #time value of money (discount rate) [-]
CRF = (r*(1+r)**lifetime) / ((1+r)**lifetime-1); #capital recovery factor [-]
dollars_per_kWhth = 0.10 #10 cents per kWh for electricity
joule_to_kWh = 2.77778e-7 #this is C1 in the paper. Multiply with Joules to
                          #convert to kW*h.

## A vector of permeabilities that are used to calculate LCOH at
pre_calc_perm = [
          1.0e-12,
          9.75e-13,9.5e-13,9.25e-13,9.0e-13,8.75e-13,8.5e-13,8.25e-13,8.0e-13,
          7.75e-13,7.5e-13,7.25e-13,7.0e-13,6.75e-13,6.5e-13,6.25e-13,6.0e-13,
          5.75e-13,5.5e-13,5.25e-13,5.0e-13,4.75e-13,4.5e-13,4.25e-13,4.0e-13,
          3.75e-13,3.5e-13,3.25e-13,3.0e-13,2.75e-13,2.5e-13,2.25e-13,2.0e-13,
          1.75e-13,1.6e-13,1.5e-13,1.25e-13,1.1e-13,1.0e-13,
          9.75e-14,9.5e-14,9.25e-14,9.0e-14,8.75e-14,8.5e-14,8.25e-14,8.0e-14,
          7.75e-14,7.5e-14,7.25e-14,7.0e-14,6.75e-14,6.5e-14,6.25e-14,6.0e-14,
          5.75e-14,5.5e-14,5.25e-14,5.0e-14,4.75e-14,4.5e-14,4.25e-14,4.0e-14,
          3.75e-14,3.5e-14,3.25e-14,3.0e-14,2.75e-14,2.5e-14,2.25e-14,2.0e-14,
          1.9e-14,1.8e-14,1.7e-14,1.6e-14,1.5e-14,1.4e-14,1.3e-14,1.2e-14,
          1.1e-14,1.0e-14,
          9.75e-15,9.5e-15,9.25e-15,9.0e-15,8.75e-15,8.5e-15,8.25e-15,8.0e-15,
          7.75e-15,7.5e-15,7.25e-15,7.0e-15,6.75e-15,6.5e-15,6.25e-15,6.0e-15,
          5.75e-15,5.5e-15,5.25e-15,5.0e-15,4.75e-15,4.5e-15,4.25e-15,4.0e-15,
          3.75e-15,3.5e-15,3.25e-15,3.0e-15,2.75e-15,2.5e-15,2.25e-15,2.0e-15,
          1.75e-15,1.5e-15,1.25e-15,1.0e-15
              ]

## A vector of permeabilities that are used to calculate LCOH at
pre_calc_perm_long = [
          9.75e-10,9.5e-10,9.25e-10,9.0e-10,8.75e-10,8.5e-10,8.25e-10,8.0e-10,
          7.75e-10,7.5e-10,7.25e-10,7.0e-10,6.75e-10,6.5e-10,6.25e-10,6.0e-10,
          5.75e-10,5.5e-10,5.25e-10,5.0e-10,4.75e-10,4.5e-10,4.25e-10,4.0e-10,
          3.75e-10,3.5e-10,3.25e-10,3.0e-10,2.75e-10,2.5e-10,2.25e-10,2.0e-10,
          1.75e-10,1.5e-10,1.25e-10,1.0e-10,
          9.75e-11,9.5e-11,9.25e-11,9.0e-11,8.75e-11,8.5e-11,8.25e-11,8.0e-11,
          7.75e-11,7.5e-11,7.25e-11,7.0e-11,6.75e-11,6.5e-11,6.25e-11,6.0e-11,
          5.75e-11,5.5e-11,5.25e-11,5.0e-11,4.75e-11,4.5e-11,4.25e-11,4.0e-11,
          3.75e-11,3.5e-11,3.25e-11,3.0e-11,2.75e-11,2.5e-11,2.25e-11,2.0e-11,
          1.75e-11,1.5e-11,1.25e-11,1.0e-11,
          9.75e-12,9.5e-12,9.25e-12,9.0e-12,8.75e-12,8.5e-12,8.25e-12,8.0e-12,
          7.75e-12,7.5e-12,7.25e-12,7.0e-12,6.75e-12,6.5e-12,6.25e-12,6.0e-12,
          5.75e-12,5.5e-12,5.25e-12,5.0e-12,4.75e-12,4.5e-12,4.25e-12,4.0e-12,
          3.75e-12,3.5e-12,3.25e-12,3.0e-12,2.75e-12,2.5e-12,2.25e-12,2.0e-12,
          1.75e-12,1.5e-12,1.25e-12,1.0e-12,
          9.75e-13,9.5e-13,9.25e-13,9.0e-13,8.75e-13,8.5e-13,8.25e-13,8.0e-13,
          7.75e-13,7.5e-13,7.25e-13,7.0e-13,6.75e-13,6.5e-13,6.25e-13,6.0e-13,
          5.75e-13,5.5e-13,5.25e-13,5.0e-13,4.75e-13,4.5e-13,4.25e-13,4.0e-13,
          3.75e-13,3.5e-13,3.25e-13,3.0e-13,2.75e-13,2.5e-13,2.25e-13,2.0e-13,
          1.9e-13,1.8e-13,1.7e-13,1.6e-13,1.5e-13,1.4e-13,1.3e-13,1.2e-13,
          1.1e-13,1.0e-13,
          9.75e-14,9.5e-14,9.25e-14,9.0e-14,8.75e-14,8.5e-14,8.25e-14,8.0e-14,
          7.75e-14,7.5e-14,7.25e-14,7.0e-14,6.75e-14,6.5e-14,6.25e-14,6.0e-14,
          5.75e-14,5.5e-14,5.25e-14,5.0e-14,4.75e-14,4.5e-14,4.25e-14,4.0e-14,
          3.75e-14,3.5e-14,3.25e-14,3.0e-14,2.75e-14,2.5e-14,2.25e-14,2.0e-14,
          1.75e-14,1.5e-14,1.25e-14,1.0e-14,
          9.75e-15,9.5e-15,9.25e-15,9.0e-15,8.75e-15,8.5e-15,8.25e-15,8.0e-15,
          7.75e-15,7.5e-15,7.25e-15,7.0e-15,6.75e-15,6.5e-15,6.25e-15,6.0e-15,
          5.75e-15,5.5e-15,5.25e-15,5.0e-15,4.75e-15,4.5e-15,4.25e-15,4.0e-15,
          3.75e-15,3.5e-15,3.25e-15,3.0e-15,2.75e-15,2.5e-15,2.25e-15,2.0e-15,
          1.75e-15,1.5e-15,1.25e-15,1.0e-15
              ]

### dT = T_CV - T_DH, which depends on depth.
###    Only important for COP caluclations in first figure
# dT = 33;  # reservoir_depth = 200
# dT = 35; # reservoir_depth = 575
# dT = 36.2; #reservoir_depth = 800
# dT = 39.5; # reservoir_depth = 1500
