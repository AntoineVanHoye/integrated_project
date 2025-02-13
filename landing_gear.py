import numpy as np
import matplotlib.pyplot as plt
import math
import initmass as init
import weight as whe
import pandas as pd
import matplotlib.patches as patches


#geometry
x_cg = 8.980632020973697        # coordinate in x of the center of gravity (estimate). Y coordinate will for now assumed to be on the chord
y_cg = 0                        # coordinate in y of the center of gravity (estimate)
Vstall = 52.57* 3.2808399       # stalling velocity with flaps(m/s)
Cl0 = 0.244                     # lift coefficient at 0 degree angle
a = 1.754                       # Cl derivative in alpha
S = 204.41                      # Surface of the wing m^2

x_fus = 16.8                    # lenght of fuselage

#takeoff general data
theta_LOF = 15                  # maximum pitch angle margin to take before takeoff
theta_margin = 0
theta_TOA = 45                  # turnover angle
S = 1829                        # maximum airplane runway (less than 1800 m idealy)
rhoair_grnd_20 = 1.395          # air density on the ground kg/m3 at 20° Fahrenheit
rhoair_grnd_85 = 1.146          # air density on the ground kg/m3 at 85° Fahrenheit
g = 9.81
w_naft = 0.08                   # weight fraction on the forward lg at takeoff
w_maft = 1 - w_naft             # weight fraction on the rear lg at takeoff
w_nfwd = 0.15                   # weight fraction on the forward lg for the most forward position of cg
w_mfwd = 1 - w_nfwd             # weight fraction on the rear lg for the most forward position of cg

#Determination of key variables

#Theta =         # Maximum pitch angle
#Phi =           # Maximum roll angle

#Gamma =         # Dihedral angle (discuss with Antoine)
#Hg =            # Wing height
#t =             # Distance between landing gears


V_LOF = 1.1 * Vstall        # Lift-off speed(emprirical formula)

Cl_20 = 2*(init.W0*453592/1000000*g)/(rhoair_grnd_20*V_LOF**2*S)    # Cl at lift-off for 20° / 85°
Cl_85 = 2*(init.W0*453592/1000000*g)/(rhoair_grnd_85*V_LOF**2*S)    


AlphaLOF_20 = (Cl_20-Cl0)/(a)*180/np.pi      # Angle of attack at lift-off for 20° / 85°
AlphaLOF_85 = (Cl_85-Cl0)/(a)*180/np.pi
print(f"At 20 degree F :{AlphaLOF_20}",f"At 85 degree F :{AlphaLOF_85}")

#ThetaTD =       # Touch down angle

#approx : we neglect the diameter of the wheel for now

#lm =            # Distance between plane gravity center and aft landing gear
#Zcg =           # Plane gravity center height


 #affichage figure
    
    # Load coordinates from the remaining lines
points = pd.read_csv(r"Airfoils\NACA45118_XYZ.csv", sep=",", header=None).to_numpy()
angle = np.radians(2)                                   # attack angle

# Matrice de rotation (attack angle)
rotation_matrix = np.array([
    [np.cos(angle), -np.sin(angle)],
    [np.sin(angle), np.cos(angle)]
])
points = np.dot(points[:,:2]*x_fus, rotation_matrix)    #proportionnal to the fuselage lenght
# Extract x and y coordinates

x, y = points[:, 0], points[:, 1]



points2 = pd.read_csv(r"Airfoils\NACA45118_XYZ.csv", sep=",", header=None).to_numpy()

angle = np.radians(2)                                   # attack angle

# Matrice de rotation (attack angle)
rotation_matrix = np.array([
    [np.cos(angle), -np.sin(angle)],
    [np.sin(angle), np.cos(angle)]
])
points2 = np.dot(points2[:,:2]*10.3733339696600000000, rotation_matrix)    #proportionnal to the fuselage lenght

transla_b = np.ones((len(points2[:,0]))) *6.42666603034
x_b, y_b = points2[:, 0]+ transla_b , points2[:, 1] #draw root wing limit



points3 = pd.read_csv(r"Airfoils\NASASC(2)0010.csv", sep=",", header=None).to_numpy()
angle = np.radians(1)                                   # attack angle

# Matrice de rotation (attack angle)
rotation_matrix = np.array([
    [np.cos(angle), -np.sin(angle)],
    [np.sin(angle), np.cos(angle)]
])
points3 = np.dot(points3[:,:2]*1.497, rotation_matrix)    #proportionnal to the fuselage lenght


transla = np.ones((len(points3[:,0]))) *10.56880165407
x_t, y_t = points3[:, 0] + transla , points3[:, 1] #draw the rear wing limit


H = (x[-1]-x_cg-np.tan(np.pi/2-np.radians(theta_LOF))*(-y_cg + y[-1]))/(np.tan(np.radians(theta_LOF+theta_margin))+np.tan(np.pi/2-np.radians(theta_LOF)))
h = H * np.tan(np.radians(theta_LOF+theta_margin))


x_w, y_w = x_cg + h , y_cg - H
x_n = x_w -(x_w-x_cg)/w_naft
y_n = y_w
print(x_n/16.8)
x_cgfwd = x_w - w_nfwd *(x_w-x_n)
y_cgfwd = y_cg


print(f'Coordinates of forward landing gear (m): {x_w,y_w}')





#Uper view of landing gear repartition

#plot fuselage lines

plt.figure(figsize=(8, 4))
z_cg = 0
z = np.zeros((len(x)))
z_b = np.ones((len(x_b))) * 4.5
z_t = np.ones((len(x_t))) * 14.5
z_n = 0
plt.plot(z*3.28084,x*3.28084, marker=None, markersize=2, linestyle='-', color='blue')

plt.plot(z_t*3.28084 ,x_t*3.28084,  marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(-z_t*3.28084 ,x_t*3.28084, marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(z_b*3.28084 ,x_b*3.28084, marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(-z_b*3.28084 ,x_b*3.28084, marker=None, markersize=2, linestyle='-', color='blue')

line = [[z[len(x)//2],z_b[len(x)//2]],[x[len(x)//2],x_b[len(x)//2]]]
line_minus = [[-z[len(x)//2],-z_b[len(x)//2]],[x[len(x)//2],x_b[len(x)//2]]]
plt.plot(np.array(line[0])*3.28084 ,np.array(line[1])*3.28084,  marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(np.array(line_minus[0])*3.28084 ,np.array(line_minus[1])*3.28084,  marker=None, markersize=2, linestyle='-', color='blue')

line = [[z_b[len(x)//2],z_t[len(x)//2]],[x_b[len(x)//2],x_t[len(x)//2]]]
line_minus = [[-z_b[len(x)//2],-z_t[len(x)//2]],[x_b[len(x)//2],x_t[len(x)//2]]]
plt.plot(np.array(line[0])*3.28084 ,np.array(line[1])*3.28084,  marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(np.array(line_minus[0])*3.28084 ,np.array(line_minus[1])*3.28084,  marker=None, markersize=2, linestyle='-', color='blue')

line = [[z[0],z_b[0]],[x[0],x_b[0]]]
line_minus = [[-z[0],-z_b[0]],[x[0],x_b[0]]]
plt.plot(np.array(line[0])*3.28084 ,np.array(line[1])*3.28084,  marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(np.array(line_minus[0])*3.28084 ,np.array(line_minus[1])*3.28084,  marker=None, markersize=2, linestyle='-', color='blue')

line = [[z_t[0],z_b[len(x)//2]],[x_t[0],x_b[len(x)//2]+4.622]]
line_minus = [[-z_t[0],-z_b[len(x)//2]],[x_t[0],x_b[len(x)//2]+4.622]]
plt.plot(np.array(line[0])*3.28084 ,np.array(line[1])*3.28084,  marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(np.array(line_minus[0])*3.28084 ,np.array(line_minus[1])*3.28084,  marker=None, markersize=2, linestyle='-', color='blue')

# Draw cg and landing gear
z_n = 0
z_w =z_n + (x_w-x_n)*np.tan(np.arcsin((y_cg-y_w)/(x_cg-x_n)*1/(np.tan(np.radians(theta_TOA)))))
plt.plot(z_cg*3.28084,x_cg*3.28084,marker='o', markersize=7,color= 'red',label ='cg.aft')
plt.plot(z_n*3.28084,x_n*3.28084,marker='x', markersize=7,color= 'black',label ='landing gear')
plt.plot(z_w*3.28084,x_w*3.28084,marker='x', markersize=7,color= 'black')
plt.plot(-z_w*3.28084,x_w*3.28084,marker='x', markersize=7,color= 'black')
plt.plot(z_cg*3.28084,x_cgfwd*3.28084 ,marker='x', markersize=7, linestyle='',color= 'red',label = 'cg.fwd(max)')

ax = plt.gca()
ax.annotate('', xy=(-z_w*3.28084, x_w*3.28084), xytext=(z_w*3.28084, x_w*3.28084),arrowprops=dict(arrowstyle='<|-|>', color='gray', linewidth=2))
plt.text((0)/2*3.28084, (x_w)*3.28084, 'D_lg', fontsize=12, ha='right', va='bottom', color='black')

plt.legend()
plt.title(f"Latteral disposition of the landing gear")
plt.xlabel("y-[inch]")
plt.ylabel("x-[inch]")
plt.axis("equal")  # Ensure equal scaling on both axes
plt.grid(True)
plt.savefig(r'Plots\Latteral_disposition_landing_gear.pdf',format='pdf')

# Weight determination (statistical)

W_l = (init.W0)
#method 1
Kfl = 0.11
Krf = 0.15
Ksl = 1
Klg = Ksl + Krf + Kfl

w_g1 = 0.046*Klg
w_nose = w_nfwd * w_g1
w_rear = w_mfwd * w_g1
print(w_g1,w_nose,w_rear)


#method 2 (statist)
n= 1.17
Kg = 20.45

#table values
Kcg = 0.718

w_g2 = Kg * Kcg*(W_l/1000)**n*1/W_l
w_nose = w_nfwd * w_g1
w_rear = w_mfwd * w_g1
print(w_g2,w_nose,w_rear)

# wheel

M_sl_m = w_maft* W_l/2* 1.07*1.25                           # maximum static load main gear
M_sl_n = w_nfwd* W_l* 1.07*1.25                             # maximum static load nose gear
M_bl_n = (M_sl_n + (10*H*W_l)/(32.2*(x_w-x_n)))* 1.07*1.25  # maximum braking load nose gear

N_w_main = 2
N_w_nose = 2

M_w_m = M_sl_m / N_w_main
M_w_n = M_sl_n / N_w_nose
M_wb_n = M_bl_n / N_w_nose
M_wb_n2 = M_wb_n /1.3
print(f'Total weight: {W_l}, Static load main weel: {M_w_m}, Static load nose wheel: {M_w_n}, Braking load nose gear: {M_wb_n} (all lbs), Braking load nose gear (for diameter): {M_wb_n2} (all lbs)')

# main wheel :
D_max_m = 39.80     #inch
W_max_m = 14        #inch
D_wheel_m = 16      #inch
Rr_m = 16.5         # rolling radius main gear

# nose wheel :
D_max_n = 29.40     #inch
W_max_n = 7.85      #inch
D_wheel_n = 16      #inch
Rr_n = 12.7         # rolling radius nose gear

# braking

KE_b = 1/2*(W_l/32.2)*(Vstall)**2

KE_b_w = KE_b /4

print(f'Braking kinetic energy per wheel : {KE_b_w}')


S = 3625.94344325

N = np.int64((30*W_l*Rr_m/12)/(4*np.pi*((8/12)**3-(2/12)**3)*S*4*32.2)+1)
print(f'n: {N}')
# damping (oleo pneumatic)

V_v = 10            # ft/s safety limit for vertical landing speed to absorb
nu =0.75            # efficiency of the shock absorber
nu_t = 0.47         # efficiency of the tire
Ngear = 3          # Gear load factor


KE_v = 1/2*(W_l/32.2)*(V_v)**2 # kinetic energy to absorb

S = ((V_v**2)/(2*32.2*nu*Ngear)-nu_t/nu*(D_max_m/2-Rr_m)/12)*12 +1 #stroke of the shock absorber (inches)

print(f'Kinetic energy to absorb: {KE_v}, Stroke of the shock absorber: {S}')

#supposing minimal lenght of 8 inches of the shock absorber 
S = 8
L = S * 2.5
print(f'Length of the shock absorber: {L}')
#D_m_d = 1.3*((4*M_sl_m)/(1800*np.pi))**1/2      #diameter damping main gear
D_m_d = 0.04*(M_sl_m)**(1/2)
#D_n_d = 1.3*((4*M_sl_n)/(1800*np.pi))**1/2      #diameter damping nose gear
D_n_d = 0.04*(M_bl_n)**(1/2)
print(f'Diameter damping main gear: {D_m_d}, Diameter damping nose gear: {D_n_d}')


#CAD
M =  0.453592 * np.array(whe.get_weight())
#print (M)
Vol_aile = 3.98238578552851391
Vol_aile2 = 3.98921089024740887
Vol_fus = 12.04511198465683556
densite_aile = M[2]/ (2*Vol_aile)
densite_aile2 = M[2]/ (2*Vol_aile2)


#print(densite_aile)
#print(densite_aile2)

# Plot the airfoil
plt.figure(figsize=(8, 4))
plt.plot(x*3.28084, y*3.28084, marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(x_t*3.28084, y_t*3.28084, marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(x_b*3.28084, y_b*3.28084, marker=None, markersize=2, linestyle='-', color='blue')
plt.plot(x_cg*3.28084,y_cg*3.28084,marker='o', markersize=7, linestyle='',color= 'red',label ='cg.aft')
plt.plot(x_w *3.28084,y_w*3.28084 ,marker='x', markersize=7, linestyle='',color= 'black',label ='landing gear')


plt.plot(x_n*3.28084,y_n*3.28084 ,marker='x', markersize=7,color= 'black')
plt.plot(x_cgfwd*3.28084,y_cgfwd*3.28084 ,marker='x', markersize=7, linestyle='',color= 'red',label = 'cg.fwd(max)')


#Draw lines
line = [[x[-1],x_w],[y[-1],y_w]]
line_cg = [[x_cg,x_w],[y_cg,y_w]]
plt.plot(np.array(line[0])*3.28084, np.array(line[1])*3.28084, marker=None, markersize=2, linestyle='--', linewidth=1, color='black')
plt.plot(np.array(line_cg[0])*3.28084, np.array(line_cg[1])*3.28084, marker=None, markersize=2, linestyle='--', linewidth=1, color='black')


circle_w = plt.Circle((x_w*3.28084, (y_w+(Rr_m)*0.0254)*3.28084), (Rr_m)*0.0254*3.28084, color='black', fill=False)
circle_n = plt.Circle((x_n*3.28084, (y_n+(Rr_n)*0.0254)*3.28084), (Rr_n)*0.0254*3.28084, color='black', fill=False)
circle_w_roll = plt.Circle((x_w*3.28084, (y_w+(Rr_m)*0.0254)*3.28084), (D_max_m/2)*0.0254*3.28084, color='black',linestyle='--', fill=False)
circle_n_roll = plt.Circle((x_n*3.28084, (y_n+(Rr_n)*0.0254)*3.28084), (D_max_n/2)*0.0254*3.28084, color='black',linestyle='--', fill=False)
arc = patches.Arc((x_w*3.28084,y_w*3.28084), 5*3.28084, 5*3.28084, angle=0, theta1=0, theta2=theta_LOF, color='gray', linewidth=2)
arc2 = patches.Arc((x_w*3.28084,y_w*3.28084), 5*3.28084, 5*3.28084, angle=0, theta1=90, theta2=theta_LOF+90, color='gray', linewidth=2)
ax = plt.gca()
ax.add_patch(circle_w)
ax.add_patch(circle_n)
ax.add_patch(circle_w_roll)
ax.add_patch(circle_n_roll)
ax.add_patch(arc)
ax.add_patch(arc2)
ax.annotate('', xy=(x_w*3.28084, y_w*3.28084), xytext=(x_n*3.28084, y_n*3.28084),arrowprops=dict(arrowstyle='<|-|>', color='gray', linewidth=2))
plt.text((x_w+x_n)/2*3.28084, (y_w+y_n)/2*3.28084, 't', fontsize=12, ha='center', va='bottom', color='black')

ax.annotate('', xy=(x_cg*3.28084, y_cg*3.28084), xytext=(x_cg*3.28084, y_w*3.28084),arrowprops=dict(arrowstyle='<|-|>', color='gray', linewidth=2))
plt.text((x_cg+x_cg)/2*3.28084, (y_cg+y_w)/2*3.28084, 'Zcg', fontsize=12, ha='right', va='bottom', color='black')

plt.text((x_w + 2.5)*3.28084, (y_w +0.3)*3.28084, f'θ', fontsize=12, verticalalignment='center')
plt.text((x_w -0.5)*3.28084, (y_w +3)*3.28084, f'βcg', fontsize=12, verticalalignment='center')
plt.legend()
plt.title(f"Longitudinal disposition of the landing gear")
plt.xlabel("x-[inch]")
plt.ylabel("z-[inch]")
plt.axis("equal")  # Ensure equal scaling on both axes
plt.grid(True)
plt.savefig(r'Plots\Longitudinal_disposition_landing_gear.pdf',format='pdf')



print(f't:{(x_w-x_n)*3.28084}',f'H:{(y_w-y_cg)*3.28084}',f'{((x_w-x_cg)**2+(y_w-y_cg)**2)**0.5*3.28084}',f'{2*z_w*3.28084}')

plt.show()