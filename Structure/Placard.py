import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

import matplotlib.pyplot as plt

# -----------------------------------
# Questions 
# -----------------------------------
# 1. Formula for the thrust ?
# 2. Last zone with minimum speed ? I don't really understand it ...
# 3. Reexplain the dive speed and which vlaue do we have to take ?
# 4. Turbulent zone and its application on the dive velocity ? 



# Constants
R = 287  # [m^2/s^2K]
gamma = 1.4

# Assumed values
M_cruise = 0.9  # [-] Cruise Mach number (not at max thrust)
factor = 1.06 # Expl of a factor between the design cruise mach number and the mach number at maximum thrust (DEPEND ON THE ENGINE, etc), so will be changed  
M_C = M_cruise * factor # Mach obtained at maximum engine thrust (slide 7 cours structure PI)
S_wing = 141   # Wing area [m^2]
C_L = 0.33  # Coefficient of lift
C_D = 0.022 # Coefficient of drag
Weight = 471511.49122  # [N]
T_max = 160200 # Thrust of the two engines [N] (based on the one of the DassaultFalcon 10X)


speed = np.linspace(0, 500, 1000) # True airspeed values up to 500m/s for the x-axis 


# Function to calculate air density using the ISA model
def air_density(altitude):
    # Up to 11 km (Troposphere)
    if altitude <= 11000:
        T = 288.15 - 0.0065 * altitude  # Temperature [K]
        p = 101325 * (T / 288.15) ** 5.2561  # Pressure [Pa]
    else:
        # Simplification for stratosphere, constant T [K] above 11 km
        T = 216.65  # Constant temperature [K]
        p = 22632 * np.exp(-9.81 * (altitude - 11000) / (287.05 * T))
    rho = p / (287.05 * T)  # Air density [kg/m^3]
    return rho, T

# Function to calculate the true airspeed at a given altitude
def true_airspeed_at_altitude(altitude):
    T = air_density(altitude)[1]
    a = np.sqrt(gamma * R * T)  # [m/s] Speed of sound
    v = M_C * a  # [m/s] Aircraft velocity
    return v








"""
# ------------------------
# Cruise speed (V_C)
# ------------------------
"""

# Define the different altitude zones 
#altitude_1 = np.linspace(12500, 14000, 1500) (Defined after)  # Max alt (assume 14000 [m], but not actually true) -> Design alt (12500[m])
altitude_2 = np.linspace(11000, 12500, 1500)  # Design alt (12500[m]) -> Stratosphere (11000[m])
altitude_3 = np.linspace(0, 11000, 11000)     # Stratosphere (11000[m]) -> Stall limit (at 0[m] for now, not true !)
altitude_4 = np.linspace(0, 0, 0)             # Stall limit (at 0[m] for now, not true !) -> Ground at 0[m]




# 1.Lift and thrust limit : 
# -------------------------    
 
   
# The lift limit is :
def lift_at_altitude(altitude):
    rho, T = air_density(altitude)
    a = np.sqrt(gamma * R * T)  # [m/s] Speed of sound
    v = M_C * a  # [m/s] Aircraft velocity
    L = 0.5 * rho * v**2 * S_wing * C_L
    return L

# Define the difference between lift and weight
def lift_minus_weight(altitude):
    return lift_at_altitude(altitude) - Weight

# Find the root of the equation lift - weight = 0
result_1 = root_scalar(lift_minus_weight, bracket=[0, 20000], method='brentq')

# Extract the altitude of intersection and print it
if result_1.converged:
    intersection_altitude_lift = result_1.root
    print(f"Altitude of intersection for lift: {intersection_altitude_lift:.2f} meters")
else:
    print("Root-finding did not converge.")
# My altitude obtained with the lift method is the variable ' intersection_altitude_lift'


# The thrust limit : NB : I think I need a relationship between the produce thrust and the density of the air because here my method isn't correct because I didn't consider the fact that the thrust is dependant of the altitude (and so the density fo the air) 
"""
def drag_at_altitude(altitude):
    rho, T = air_density(altitude)
    a = np.sqrt(gamma * R * T)  # [m/s] Speed of sound
    v = M_C * a  # [m/s] Aircraft velocity
    D = 0.5 * rho * v**2 * S_wing * C_D
    return D
# Define the difference between drag and thrust
def drag_minus_thrust(altitude):
    return drag_at_altitude(altitude) - T_max

# Find the root of the equation drag - thrust = 0
result_2 = root_scalar(drag_minus_thrust, bracket=[0, 20000], method='brentq')

# Extract the altitude of intersection and print it
if result_2.converged:
    intersection_altitude_drag = result_2.root
    print(f"Altitude of intersection for drag: {intersection_altitude_drag:.2f} meters")
else:
    print("Root-finding did not converge.")
# My altitude obtained with the drag method is the variable ' intersection_altitude_drag'
"""

# Having no result for the thrust, I temporally take the value of the lift method
L_anf_T_limit_alt = intersection_altitude_lift

altitude_1 = np.linspace(12500, L_anf_T_limit_alt, int(L_anf_T_limit_alt-12500)) # 'Dynamic definition of the altitude range
L_anf_T_limit_alt_line = [L_anf_T_limit_alt] * len(speed) # To trace the horizontal line on the graph

# Trace the vertical line
v_true_below_LT_limit = np.linspace(12500, L_anf_T_limit_alt, int(L_anf_T_limit_alt-12500))
for i in range(len(altitude_1)) : 
    
    # Compute a value for the true airspeed with the air_density = f(altitude_1)
    v_true_below_LT_limit[i] = true_airspeed_at_altitude(altitude_1[i]) 



# 2.Design altitude 12 500 [m] :
# ------------------------------
design_alt = 12500
design_alt_line = [design_alt] * len(speed) # to trace the horizontal line

# Trace the vertical line
v_true_below_design = np.linspace(11000, 12500, 1500)
for i in range(len(altitude_2)) : 
    
    # Compute a value for the true airspeed with the air_density = f(altitude_1)
    v_true_below_design[i] = true_airspeed_at_altitude(altitude_2[i]) 



# 3.Stratospheric limit 11 000 [m] :
# ----------------------------------
stratospheric = 11000
stratospheric_line = [stratospheric] * len(speed) # to trace a line
#print('Air density at 11000 [m]',air_density(11000)[0]) # Air density at 11000 [m]
cste = 0.5*(air_density(11000)[0])*(true_airspeed_at_altitude(11000))**2 # Drag has to be kept cst => slide 10
#print('The value of the constant is', cste)
 


       
# 3.(pt2) Below stratospheric limit < 11 000 [m] : 
# ------------------------------------------------
v_true_below_strat = np.linspace(0, 11000, 11000) # same size than altitude_3

for i in range(len(altitude_3)) : 
    
    # Compute a value for the true airspeed with the air_density = f(altitude_3)
    v_true_below_strat[i] = np.sqrt(2*cste/air_density(altitude_3[i])[0]) 




# 4. Minimum speed to not stall (close to the ground) :
# -----------------------------------------------------
# Verify because the teaching assistant wasn't sure and we still have some speed left 
# But on the graph in the course we clearly see that the true airspeed becomes constant











"""
# ------------------------
# Dive speed (V_D)
# ------------------------
"""
# Design dive Mach is defined as the minimum between :
# 1.25*M_C
# Mach actually obtained after a 20-second dive at 7.5° followed by a 1.5-g pullout --> MD ~ 1.07 MC
# I don't know which one to choose ...

# What about the turbulent zone on the slides ? How do I know where the zone is ... look on the internet ?




# Plotting
plt.figure(figsize=(16, 9))
# Plotting of the horizontal lines
plt.plot(speed, L_anf_T_limit_alt_line , color='red', label='Lift anf thrust limit', linestyle='-')
plt.plot(speed, stratospheric_line, color='blue', label='Stratospheric limit', linestyle='--')
plt.plot(speed, design_alt_line, color='cyan', label='Design altitude', linestyle='--')
# Plotting V_C curve
plt.plot(v_true_below_LT_limit, altitude_1, color='blue')
plt.plot(v_true_below_design, altitude_2, color='blue')
plt.plot(v_true_below_strat,altitude_3, color='blue')
plt.title('Placard diagram', fontsize=14)
plt.xlabel('True airspeed [m/s]', fontsize=12)
plt.ylabel('Altitude [m]', fontsize=12)
#plt.grid(True)
plt.legend()
plt.show()
