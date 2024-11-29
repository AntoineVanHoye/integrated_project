import numpy as np
from matplotlib import pyplot as plt

#C_m : pitching moment coefficient
"""
data_fus = np.loadtxt("fus_data.txt", usecols=(0, 4))
data_wings = np.loadtxt("wings_data.txt", usecols=(0, 4))

AoA_fus = data_fus[0]
Cm_fus = data_fus[1]
AoA_wings = data_wings[0]
Cm_wings = data_wings[1]
"""
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
a = np.sqrt(gamma * R * T) #speed of sound 
V = a * Mach


L_wing = CL_wing * 1/2 * rho* S_wing * V**2
L_fus = CL_fus * 1/2 * rho* S_fus* V**2
L_T = 5

x_AC_wing = 2
x_CG_plane = 5
x_AC_fus = 5
x_AC_tail = 7
MAC_wing = 8
S_T = 2

#tail volume ratio effectivness 

V_T = S_T * (x_AC_tail - x_CG_plane)/(S_wing * MAC_wing)


"""
#to see the trend of Cm with respect to the angle of attack
def plot_Cm(AoA_fus, AoA_wings, Cm_fus, Cm_wings):

    alpha_fus = AoA_fus
    alpha_wings = AoA_wings

    plt.plot(alpha_fus, Cm_fus, color = 'b', label="Fuselage")
    plt.plot(alpha_wings, Cm_wings, color = 'r', label="Wings")
    plt.xlabel("Angle of attack α (°)")
    plt.ylabel("Pitching moment coefficient C_m")
    plt.axhline(0, color='gray', linestyle='--', label="C_m = 0") #to see the Cm where there is no angle of attck
    plt.legend()
    plt.show()

    return

plot_Cm(AoA_fus, AoA_wings,Cm_fus,Cm_wings)
"""
##################################################################
 ######EQUILIBRIUM
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

def long_stat_stab(): #in the pitching plane
    #check the stability
    hn = x_AC_wing + 2*V_T*(1-deps)*a1_over_a #position of the neutral point 
    derivative = hn - x_CG_plane 
    Kn = - derivative #static margin

    if Kn >= 0.05 : 
        print("The static margin has a correct value and is equal to : ", Kn)

    else : 
        print("The static margin has to be changed and is equal to : ", Kn)
    
    return 

long_stat_stab()

def instability(Xac, Xcg):
    if (Xcg - Xac < 1): 
        print("The center of gravity has a good position with respect to the neutral point") 
    else : 
        print("The center of gravity has not a good position with respect to the neutral point") 

    

 
 ##################################################################
 ######LATERAL STATIC STABILITY
 ##################################################################


def roll_stab(): #stable if dCr/dalpha < 0, phi : roll angle 

    derivative = np.diff(Cr)/np.diff(phi)
    return 

def yaw_stab(): #stable if dCn/dbeta < 0

    return 
