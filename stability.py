import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import trapz
from wings import fusGeometry as fusgeom
from wings import wingGeometry as winggeom


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

#a = (CL_wings[68]- CL_wings[67])/(AoA_wings[68]-AoA_wings[67])*np.pi/180
#a1 = (CL_tail[68]- CL_tail[67])/(AoA_tail[68]-AoA_tail[67])
a = 1.754
a1 = 1.5

a1_over_a = a1/a


#at cruise : AoA of the fuselage : O° and 6° for the wings (with a wash of 3°)
AoA_root = 6
AoA_tip = 3
AoA_wing = (AoA_root - AoA_tip)/2
Cm_fus_cruise = -0.1072 #at 0° AoA
Cm_wings_cruise =  -0.1448 #at AoA_wing
CL_wing = 0.574 #total wing lift coeff in cruise 
CL_fus = 0.056
CL_T = CL_wing*15/100

b = 29 #span
S_fus = 110.67
S_wing = 66.44
surf_tot = 182.56

width_fus = 9

Cm0_airfoil_fus = -0.1158 
Cm0_airfoil_wing = -0.1533

stop_wing = b/2 - width_fus
n = 10
"""
y_values_wing = np.linspace(0,stop_wing,n)
y_values_fus = np.linspace(stop_wing,b/2,n)
"""
_, _, _, _, _, _, _, _, y_wing, leading_wing, trailing_wing, quarter_wing = winggeom() 
_, _, _, _, _, _, _, y_fus, leading_fus, trailing_fus, quarter_fus = fusgeom() 
c_fus = trailing_fus - leading_fus
c_wing = trailing_wing - leading_wing

MAC_tot = 9.949 
Cm0_wing = (2/(surf_tot*MAC_tot)) * trapz(Cm0_airfoil_wing*c_wing**2, y_wing)
Cm0_fus = (2/(surf_tot*MAC_tot)) * trapz(Cm0_airfoil_fus*c_fus**2, y_fus)
Cm0_tot = Cm0_wing + Cm0_fus


print("The blended wing body has a Cm0_tot of",Cm0_tot)

"""
#calculate the forces 








L_wing = CL_wing * 1/2 * rho* S_wing * V**2
L_fus = CL_fus * 1/2 * rho* S_fus* V**2
L_T = 5
"""
rho = 0.288
Mach = 0.9
R = 287                 #[m^2/s^2K]
gamma = 1.4
T = 216.5
speed_sound = np.sqrt(gamma * R * T) #speed of sound 
V = speed_sound * Mach
AR_tot = 1.5
wing_AR = 6.021
#values with respect to the airplane nose 
x_AC_wing = 11.91
x_CG_plane = 10.4
x_AC_fus = 1.33
x_AC_tail = 16.5
z_AC_wing = 0
z_AC_tail = 1.51
z_AC_tot = 0
MAC_tail = 1.483
S_T = 5.422
MAC_wing = 4.141
MAC_fus = 13.305
fus_length = 16.8
fus_width = 9
x_AC_tot = 8.789
MAC_tot = 9.949

Cl_tot = 0.244
#tail volume ratio effectivness 
hor_tail_surf = 3.182
V_T = hor_tail_surf * (x_AC_tail - x_CG_plane)/(surf_tot* MAC_tot)

def plot_CL(AoA, CL): 
    plt.plot(AoA, CL, color = 'b')
    plt.xlabel("Angle of attack α (°)")
    plt.ylabel("Lift coefficient C_L")
    plt.axhline(0, color='gray', linestyle='--', label="C_L = 0") #to see the Cm where there is no angle of attck
    plt.legend()
    plt.show()

    return

#plot_CL(AoA_wings,CL_wings)
#plot_CL(AoA_tail, CL_tail)

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

#plot_Cm(AoA_fus, AoA_wings,Cm_fus,Cm_wings)

##################################################################
######EQUILIBRIUM IN PITCH
##################################################################
def required_CLT(): 
    CLT = (Cm0_tot + Cl_tot* (x_CG_plane - x_AC_tot)/MAC_tot)/V_T
    L_T = CLT * 1/2 * rho * V**2 * hor_tail_surf

    print("The required lift coefficient of the tail to reach the equilibirum in pitch is CLT =", CLT, "[-] and the required lift equals ", L_T, "[N]")
    return CLT

required_CLT()

def equilibrium() : 
    
    Cm_tot = Cm0_tot + Cl_tot* (x_CG_plane - x_AC_tot)/MAC_tot - CL_T *V_T
    
    #check the equilibrium at cruise

    if Cm_tot > -1e-3 and Cm_tot < 1e-3 : 
        print("The airplane is at equilibrium because Cm_tot = 0")
    
    else : 
        print("The airplane is not at equilibrium because Cm_tot = ", Cm_tot)
    
    return 

equilibrium()

##################################################################
######LONGITUDINAL STATIC STABILITY
##################################################################

def downwash():
    It = x_AC_tail - x_AC_tot
    lamb = 0.086190476 #global taper ratio 
    m = (z_AC_tail - z_AC_tot)/b*2 #because (z_AC_tail - z_AC_wing) = mb/2

    deps = (1.75*a)/(np.pi*AR_tot*(((2*lamb*It)/b)**1.4)*(1 + np.abs(m)))
    
    return deps

eps = downwash()
print("deps/dalpha is equal to", eps)


def long_stat_stab(): #in the pitching plane
    #check the stability
    #neutral point : position of the cg in order to have the derivative equals 0
    deps = downwash()
    
    hn = x_AC_tot + V_T*a1_over_a*(1- deps) #- 0.5*fus_width**2 * fus_length/(S_wing*a*MAC_wing)#position of the neutral point  

    derivative = (x_CG_plane - hn)/MAC_tot
    Kn = - derivative #static margin

    if Kn >= 0.05 and Kn < 0.15 : 
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
