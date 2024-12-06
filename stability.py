import numpy as np
from matplotlib import pyplot as plt

#C_m : pitching moment coefficient

data_fus = np.loadtxt("fus_data.txt", delimiter=",")
data_wings = np.loadtxt("wings_data.txt", delimiter=",")
data_tail = np.loadtxt("tail_data.txt", delimiter=",")

AoA_fus = data_fus[:,0]
Cm_fus = data_fus[:,4]
AoA_wings = data_wings[:,0]
AoA_tail = data_tail[:,0]
CL_tail = data_tail[:,1]
Cm_wings = data_wings[:,4]
CL_wings = data_wings[:,1]

a = (CL_wings[68]- CL_wings[67])/(AoA_wings[68]-AoA_wings[67])
a1 = (CL_tail[68]- CL_tail[67])/(AoA_tail[68]-AoA_tail[67])

a1_over_a = a1/a


#at cruise : AoA of the fuselage : O° and 6° for the wings (with a wash of 3°)
AoA_root = 6
AoA_tip = 3
AoA_wing = (AoA_root - AoA_tip)/2
Cm_fus_cruise = -0.1072 #at 0° AoA
Cm_wings_cruise =  -0.1448 #at AoA_wing
CL_wing = 0.675
CL_fus = 0.059
CL_T = 7

#calculate the forces 
S_fus = 110.67
S_wing = 54.22
rho = 0.288
Mach = 0.9
R = 287                 #[m^2/s^2K]
gamma = 1.4
T = 216.5
speed_sound = np.sqrt(gamma * R * T) #speed of sound 
V = speed_sound * Mach
wing_AR = 6.021

L_wing = CL_wing * 1/2 * rho* S_wing * V**2
L_fus = CL_fus * 1/2 * rho* S_fus* V**2
L_T = 5

x_AC_wing = 2
x_CG_plane = 11
x_AC_fus = 5
x_AC_tail = 7
y_AC_wing = 3
y_AC_tail = 4
MAC_wing = 8
S_T = 2

#tail volume ratio effectivness 

V_T = S_T * (x_AC_tail - x_CG_plane)/(S_wing * MAC_wing)

def plot_CL(AoA, CL): 
    plt.plot(AoA, CL, color = 'b')
    plt.xlabel("Angle of attack α (°)")
    plt.ylabel("Lift coefficient C_L")
    plt.axhline(0, color='gray', linestyle='--', label="C_L = 0") #to see the Cm where there is no angle of attck
    plt.legend()
    plt.show()

    return

plot_CL(AoA_wings,CL_wings)
plot_CL(AoA_tail, CL_tail)

#to see the trend of Cm with respect to the angle of attack
def plot_Cm(AoA_fus, AoA_wings, Cm_fus, Cm_wings):

    alpha_fus = AoA_fus
    alpha_wings = AoA_wings

    plt.plot(alpha_fus, Cm_fus, color = 'b', label="Fuselage")
    plt.xlabel("Angle of attack α (°)")
    plt.ylabel("Pitching moment coefficient C_m")
    plt.axhline(0, color='gray', linestyle='--', label="C_m = 0") #to see the Cm where there is no angle of attck
    plt.legend()
    plt.show()

    plt.plot(alpha_wings, Cm_wings, color = 'r', label="Wings")
    plt.xlabel("Angle of attack α (°)")
    plt.ylabel("Pitching moment coefficient C_m")
    plt.axhline(0, color='gray', linestyle='--', label="C_m = 0") #to see the Cm where there is no angle of attck
    plt.legend()
    plt.show()

    return

plot_Cm(AoA_fus, AoA_wings,Cm_fus,Cm_wings)

##################################################################
######EQUILIBRIUM IN PITCH
##################################################################

def equilibrium() : 
    Cm_tot = Cm_fus_cruise + 2* Cm_wings_cruise + 2 * CL_wing* (x_AC_wing - x_CG_plane) + CL_fus * (x_AC_fus - x_CG_plane) - 2*CL_T *V_T

    #check the equilibrium at cruise

    if Cm_tot == 0 : 
        print("The airplane is at equilibrium because Cm_tot = 0")
    
    else : 
        print("The airplane is not at equilibrium because Cm_tot = ", Cm_tot)
    
    return 

equilibrium()

##################################################################
######LONGITUDINAL STATIC STABILITY
##################################################################

def downwash():
    It = x_AC_tail - x_AC_wing
    lamb = 0.137
    b = 29
    m = (y_AC_tail - y_AC_wing)/b*2

    deps = (1.75*a)/(np.pi*wing_AR*((2*lamb*It)/b)**1.4*(1 + m))
    
    return deps

eps = downwash()

def long_stat_stab(): #in the pitching plane
    #check the stability
    deps = downwash()
    hn = x_AC_wing + 2*V_T*a1_over_a*(1- deps) #position of the neutral point 
    derivative = hn - x_CG_plane 
    Kn = - derivative #static margin

    if Kn >= 0.05 : 
        print("The static margin has a correct value and is equal to : ", (Kn*100))

    else : 
        print("The static margin has to be changed and is equal to : ", (Kn*100))
    
    return 

long_stat_stab()
##################################################################
######LATERAL STATIC STABILITY
##################################################################

def yaw_stab(): #stable if dCn/dbeta < 0

    return 
