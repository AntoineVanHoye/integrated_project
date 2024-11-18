import Numpy as np 


# Variables

M_cruise = 0.85
Max_Thrust = 75*2 # ? [kN] (assumption made, bombardier global 8000)
D = 75000 * 2 # [N] because of the max_thrust
alt_design = 12500 # [m] (41010,499 ft)
S = 40 # Front Surface [m^2] 
c_d = 0.031 # Drag coefficent
V_stall = 1 # [m/s]

# Fct that gives me the density of the air in fct of the altitude
def air_density(altitude):
    
    # Constants for ISA atmosphere up to 11,000 meters (Troposphere)
    sea_level_density = 1.225  # kg/m^3 at sea level
    sea_level_temperature = 288.15  # Kelvin
    lapse_rate = -0.0065  # Temperature lapse rate in K/m
    gas_constant = 287.05  # Specific gas constant for dry air in J/(kg·K)
    gravity = 9.80665  # Gravitational acceleration in m/s^2
    
    # Base height for the troposphere model
    tropopause_height = 11000  # meters
    
    if altitude <= tropopause_height:
        # Troposphere model (up to 11,000 m)
        temperature = sea_level_temperature + lapse_rate * altitude
        density = sea_level_density * (temperature / sea_level_temperature) ** ((-gravity) / (lapse_rate * gas_constant))
    else:
        # Stratosphere model (above 11,000 m) - simplified constant temperature
        temperature = sea_level_temperature + lapse_rate * tropopause_height
        exp_factor = (-gravity * (altitude - tropopause_height)) / (gas_constant * temperature)
        density = sea_level_density * (temperature / sea_level_temperature) ** ((-gravity) / (lapse_rate * gas_constant))
        density *= math.exp(exp_factor)
    
    return density



# ---------------
# Placard diagram 
# ---------------

# Design cruise Mach (M_c), Mach you obtained at maximum thrust

M_c = 1.06 * M_cruise 

# Speed at design altitude (max trust so T = D)

#Find the rho at alt_design
rho = air_density(alt_design)
v = np.sqrt((2 * D)/(rho * S *c_d)) # [m/s]

# Service ceiling (compute the rho at stall speed, then deduce the altitude)

rho_stall = (2*D)/((V_stall**2) * S * c_d )

alt_max = air_density(rho_stall)



