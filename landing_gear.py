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

x_fus = 16.8                # lenght of fuselage

#takeoff general data
theta_LOF = 15                  # maximum pitch angle margin to take before takeoff
theta_TOA =                     # turnover angle
S = 1829                        # maximum airplane runway (less than 1800 m idealy)
rhoair_grnd_20 = 1.395          # air density on the ground kg/m3 at 20° Fahrenheit
rhoair_grnd_85 = 1.146          # air density on the ground kg/m3 at 85° Fahrenheit
g = 9.81
w_naft = 0.08                       # weight fraction on the forward lg at takeoff
w_maft = 1 - w_naft                    # weight fraction on the rear lg at takeoff
w_nfwd = 0.15                   # weight fraction on the forward lg for the most forward position of cg
w_mfwd = 1 - w_nfwd             # weight fraction on the rear lg for the most forward position of cg
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

transla = np.ones((len(x))) *16.185
transla_b = np.ones((len(x))) *7.794
x_t, y_t = points[:, 0]*1.448 + transla , points[:, 1]*1.448 #draw the rear wing limit
x_b, y_b = points[:, 0]*9.006+ transla_b , points[:, 1]*9.006

H = (x_t[-1]-x_cg-np.tan(np.pi/2-np.radians(theta_LOF))*(y_cg + y_t[-1]))/(np.tan(np.radians(theta_LOF))+np.tan(np.pi/2-np.radians(theta_LOF)))
h = H * np.tan(np.radians(theta_LOF))


x_w, y_w = x_cg + h , y_cg - H
x_n = x_w -(x_w-x_cg)/w_naft
y_n = y_w
x_cgfwd = x_w - w_nfwd *(x_w-x_n)
y_cgfwd = y_cg


print(f'Coordinates of forward landing gear (m): {x_w,y_w}')

# Plot the airfoil
plt.figure(figsize=(8, 4))
plt.plot(x, y, marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(x_t, y_t, marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(x_b, y_b, marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(x_cg,y_cg,marker='o', markersize=7,color= 'red',label ='cg.aft')
plt.plot(x_w ,y_w ,marker='o', markersize=7,color= 'black')

plt.plot(x_n,y_n ,marker='o', markersize=7,color= 'black')
plt.plot(x_cgfwd,y_cgfwd ,marker='x', markersize=7,color= 'red',label = 'cg.fwd')


#Draw lines
line = [[x_t[-1],x_w],[y_t[-1],y_w]]
line_cg = [[x_cg,x_w],[y_cg,y_w]]
plt.plot(line[0], line[1], marker=None, markersize=2, linestyle='--', linewidth=1, color='black')
plt.plot(line_cg[0], line_cg[1], marker=None, markersize=2, linestyle='--', linewidth=1, color='black')

line = [[x_t[-1],x_w],[y_t[-1],y_w]]
line_cg = [[x_cg,x_w],[y_cg,y_w]]
plt.plot(line[0], line[1], marker=None, markersize=2, linestyle='--', linewidth=1, color='black')
plt.plot(line_cg[0], line_cg[1], marker=None, markersize=2, linestyle='--', linewidth=1, color='black')


plt.title(f"Longitudinal disposition of the landing gear")
plt.xlabel("x-[m]")
plt.ylabel("y-[m]")
plt.axis("equal")  # Ensure equal scaling on both axes
plt.grid(True)



#Uper view of landing gear repartition

#plot fuselage lines

plt.figure(figsize=(8, 4))
z_cg = 0
z = np.zeros((len(x)))
z_b = np.ones((len(x))) * 4.5
z_t = np.ones((len(x))) * 14.5
z_n = 0
plt.plot(z,x, marker=None, markersize=2, linestyle='-', color='blue')

plt.plot(z_t ,x_t,  marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(-z_t ,x_t, marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(z_b ,x_b, marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(-z_b ,x_b, marker=None, markersize=2, linestyle='-', color='blue')

line = [[z[len(x)//2],z_b[len(x)//2]],[x[len(x)//2],x_b[len(x)//2]]]
line_minus = [[-z[len(x)//2],-z_b[len(x)//2]],[x[len(x)//2],x_b[len(x)//2]]]
plt.plot(line[0] ,line[1],  marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(line_minus[0] ,line_minus[1],  marker=None, markersize=2, linestyle='-', color='blue')

line = [[z_b[len(x)//2],z_t[len(x)//2]],[x_b[len(x)//2],x_t[len(x)//2]]]
line_minus = [[-z_b[len(x)//2],-z_t[len(x)//2]],[x_b[len(x)//2],x_t[len(x)//2]]]
plt.plot(line[0] ,line[1],  marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(line_minus[0] ,line_minus[1],  marker=None, markersize=2, linestyle='-', color='blue')

line = [[z[0],z_b[0]],[x[0],x_b[0]]]
line_minus = [[-z[0],-z_b[0]],[x[0],x_b[0]]]
plt.plot(line[0] ,line[1],  marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(line_minus[0] ,line_minus[1],  marker=None, markersize=2, linestyle='-', color='blue')

line = [[z_t[0],z_b[len(x)//2]],[x_t[0],x_b[len(x)//2]+4.622]]
line_minus = [[-z_t[0],-z_b[len(x)//2]],[x_t[0],x_b[len(x)//2]+4.622]]
plt.plot(line[0] ,line[1],  marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(line_minus[0] ,line_minus[1],  marker=None, markersize=2, linestyle='-', color='blue')

# Draw cg and landing gear
plt.plot(z_cg,x_cg,marker='o', markersize=7,color= 'red',label ='cg.aft')
plt.plot(z_n,x_n,marker='o', markersize=7,color= 'black',label ='cg.aft')


plt.title(f"Latteral disposition of the landing gear")
plt.xlabel("y-[m]")
plt.ylabel("x-[m]")
plt.axis("equal")  # Ensure equal scaling on both axes
plt.grid(True)


plt.show()