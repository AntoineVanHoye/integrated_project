import numpy as np
import matplotlib.pyplot as plt
import math
import initmass as init

#geometry
x_cg = 10.4                 # coordinate in x of the center of gravity (estimate). Y coordinate will for now assumed to be on the chord
y_cg = 0                    # coordinate in y of the center of gravity (estimate)
Vstall = 52.57              # stalling velocity with flaps(m/s)
Cl0 = 0.244                 # lift coefficient at 0 degree angle
a = 1.754                   # Cl derivative in alpha
S = 204.41                  # Surface of the wing m^2

thetamin_lat_clear = 5      #
x_nt = 17.633               # nose-tail lenght (m)
x_fus = 16.8                # lenght of fuselage

#takeoff general data
theta_LOF = 15                  # maximum pitch angle margin to take before takeoff
S = 1829                        # maximum airplane runway (less than 1800 m idealy)
rhoair_grnd_20 = 1.395          # air density on the ground kg/m3 at 20° Fahrenheit
rhoair_grnd_85 = 1.146          # air density on the ground kg/m3 at 85° Fahrenheit
g = 9.81

'''
C_r = 4.622
C_t = 1.448
#C_R=  (S)*(C_r-C_t)/9

M = 2/3*(C_r + C_t-(C_r*C_t)/(C_r+C_t)) #for the standard straight tapered wing
#M = 1/2 * (C_r+C_t)
b = 20
H = b/8*(C_r**2+2*C_r*C_t+3*C_t**2)/(C_r**2+C_r*C_t+C_t**2)

print(M,H)
'''

#Determination of key variables

#Theta =         # Maximum pitch angle
#Phi =           # Maximum roll angle

#Gamma =         # Dihedral angle (discuss with Antoine)
#Hg =            # Wing height
#t =             # Distance between landing gears


V_LOF = 1.1 * Vstall        # Lift-off speed(emprirical formula)

Cl_20 = 2*(init.W0*1000000/453592*g)/(rhoair_grnd_20*V_LOF**2*S)    # Cl at lift-off for 20° / 85°
Cl_85 = 2*(init.W0*1000000/453592*g)/(rhoair_grnd_85*V_LOF**2*S)    


AlphaLOF_20 = (Cl_20-Cl0)/(a)*180/np.pi      # Angle of attack at lift-off for 20° / 85°
AlphaLOF_85 = (Cl_85-Cl0)/(a)*180/np.pi
print(f"At 20 degree F :{AlphaLOF_20}",f"At 85 degree F :{AlphaLOF_85}")

#ThetaTD =       # Touch down angle

#approx : we neglect the diameter of the wheel for now
H = (x_nt-x_cg-np.tan(np.pi/2-np.radians(theta_LOF))*y_cg)/(np.tan(np.radians(theta_LOF))+np.tan(np.pi/2-np.radians(theta_LOF)))
h = H * np.tan(np.radians(theta_LOF))
print(H,h)

#lm =            # Distance between plane gravity center and aft landing gear
#Zcg =           # Plane gravity center height


 #affichage figure

# Load the points from a Selig-format file
with open("sc2018.txt", "r") as file:
    # Skip the first line (airfoil name)
    airfoil_name = file.readline().strip()
    
    # Load coordinates from the remaining lines
    points = np.loadtxt(file)

# Extract x and y coordinates
x, y = points[:, 0]*x_fus, points[:, 1]*x_fus

transla = np.ones((len(x),1)) *16.185
x_t, y_t = points[:, 0]*1.448 + transla , points[:, 1]*1.448
# Plot the airfoil
plt.figure(figsize=(8, 4))
plt.plot(x, y, marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(x_t, y_t, marker=None, markersize=2, linestyle='-', color='red')
plt.title(f"Airfoil Shape: {airfoil_name}")
plt.xlabel("x-coordinate")
plt.ylabel("y-coordinate")
plt.axis("equal")  # Ensure equal scaling on both axes
plt.grid(True)
plt.show()