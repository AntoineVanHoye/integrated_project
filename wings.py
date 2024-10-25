import numpy as np
from matplotlib import pyplot as plt


#---Guesses---#
b = 31.6992         #Span (guest) 104[ft] from Global8000
AR = 10              #Aspect ratio (guest)
S = 0               #total area
weight = 471511.49122 #[N] = 106000lb (guess from weight code)
alti = 12500        #[m]
rho = 0.288         #[kg/m^3]
M = 0.85            #[-] Mach number
R = 287             #[m^2/s^2K]
gamma = 1.4
e = 0.85            #Ostxald's efficiency factor
CD0 = 0.02          # Zero lift drag coeff
sweep_quater = 45          #[°] sweep angle
Lambda = 0.3        # [-] taper ratio

#---Commande---#
polar_Cl_Cd = True


#---Code---#


CL = np.linspace(0, 1.0, 100)

CD = CD0 + CL*CL / (np.pi * AR * 0.85) # (Zero Lift drag + Induce drag)

CL_CD = CL / CD
# Trouver l'indice du maximum de CL/CD
max_index = np.argmax(CL_CD)

# Extraire le CL correspondant au maximum de CL/CD
CL_max = CL[max_index]
CL_CD_max = CL_CD[max_index]

# Afficher les résultats
print(f"Le maximum de CL/CD est {CL_CD_max:.2f} pour CL = {CL_max:.2f}")

if polar_Cl_Cd:
    plt.plot(CL, CL/CD)
    plt.scatter(CL_max, CL_CD_max, marker="x", color="r")
    plt.xlabel('$CL$')
    plt.ylabel('$CL/CD$')
    plt.show()

Cl_max = CL_max -  (CL_max*0.1)


T = 288.15 - (6.5e-3*alti)

a = np.sqrt(gamma * R * T) #[m^2] speed of sound
v = M*a #[m/s]


#cl_cd_ratio = 18

S = 2*weight/(rho*(v**2)*Cl_max)

c_quater = S/b
c_root = c_quater/(2*(1-Lambda))
c_tip = Lambda*c_root


#S = b**2/AR


print("Cl max used:", Cl_max, "[-]")
print("Total area:", S, "[m^2]")
print("Corde:", c_quater, "[m]")

sweep_leeding = np.tan(sweep_quater) + (4/AR) * ((1-Lambda)/(1+Lambda) * (0.25 - 0))
sweep_trailing = np.tan(sweep_quater) + (4/AR) * ((1-Lambda)/(1+Lambda) * (0.25 - 1))

x = np.linspace(0, b/2, 100)

#chord_len =  

#help of prof.
# a (CL_alpha) 
#CL = CL_alpha * (alpha - alpha_L0)

