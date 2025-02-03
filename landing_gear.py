import numpy as np
import matplotlib as plt
import math

import initmass as init

#geometry
x_cg = 10.4                 # cordinate in x of the center of gravity (estimate). Y coordinate will for now assumed to be on the chord
Vstall = 52.57              # stalling velocity with flaps(m/s)
Cl0 = 0.244                 # lift coefficient at 0 degree angle
a = 1.754                   # Cl derivative in alpha
S = 204.41                  # Surface of the wing m^2

thetamin_lat_clear = 5      #
x_nt = 17.633                   # nose-tail lenght (m)

#takeoff general data
theta_LOF = 15                  # minimum pitch angle margin to take before takeoff
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

Cl_20 = 2*(init.W0*1000000/453592*g)/(rhoair_grnd_20*V_LOF**2*S)
Cl_85 = 2*(init.W0*1000000/453592*g)/(rhoair_grnd_85*V_LOF**2*S)


AlphaLOF_20 = (Cl_20-Cl0)/(a)*180/np.pi      # Angle of attack at lift-off
AlphaLOF_85 = (Cl_85-Cl0)/(a)*180/np.pi
print(f"At 20 degree F :{AlphaLOF_20}",f"At 85 degree F :{AlphaLOF_85}")



#ThetaTD =       # Touch down angle





#lm =            # Distance between plane gravity center and aft landing gear
#Zcg =           # Plane gravity center height