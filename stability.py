import numpy as np

#C_m : pitching moment coefficient

data_fus = np.loadtxt("fus_data.txt", usecols=(0, 4))
data_wings = np.loadtxt("wings_data.txt", usecols=(0, 4))

AoA_fus = data_fus[0]
Cm_fus = data_fus[1]
AoA_wings = data_wings[0]
Cm_wings = data_wings[1]

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
    plt.grid()
    plt.show()

    return

plot_Cm(AoA_fus, AoA_wings,Cm_fus,Cm_wings)

##################################################################
 ######LONGITUDINAL STATIC STABILITY
 ##################################################################

def long_stat_stab(): #in the pitching plane
    #goal : verify if dCm/dalpha <0

    Cm_tot = Cm_fus + Cm_wings

    derivative = np.diff(Cm_tot) / np.diff(AoA_fus)

    stable = np.all(derivative < 0)

    if stable: 
        print("The stability condition dcm/dalpha is respected. The airplane is stable in pitch.")
    else : 
        print("The condition is not respected. The airplane is not stable in pitch.")

    return 

long_stat_stab()


    
def static_margin(Xac, Xcg, MAC): 
    Kn = (Xac - Xcg)/MAC

    if Kn >= 0.05 : 
        print("The static margin has a correct value and is equal to : ", Kn)

    else : 
        print("The static margin has to be changed and is equal to : ", Kn)
 
 ##################################################################
 ######LATERAL STATIC STABILITY
 ##################################################################


def roll_stab(): #stable if dCr/dalpha < 0

    derivative = np.diff(Cr)/np.diff(phi)
    return 

def yaw_stab(): #stable if dCn/dbeta < 0

    return 