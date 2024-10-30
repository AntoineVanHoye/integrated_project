import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as integrate

#---Guesses---#
b = 31.6992         #Span (guest) 104[ft] from Global8000
AR = 10              #Aspect ratio (guest)
S = 0               #total area
weight = 471511.49122 #[N] = 106000lb (guess from weight code)
alti = 12500        #[m]
rho = 0.288         #[kg/m^3]
M = 0.9            #[-] Mach number
R = 287             #[m^2/s^2K]
gamma = 1.4
e = 0.85            #Ostxald's efficiency factor
CD0 = 0.02          # Zero lift drag coeff
sweep_quarter = 10         #[°] sweep angle
Lambda = 0.6        # [-] taper ratio

#---Commande---#
polar_Cl_Cd = False
wing_plot = False

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
print(f"Le maximum de CL/CD est {CL_CD_max:.2f} pour CL = {CL_max:.2f} \n")

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

c_quarter = S/b
#c_root = c_quarter/(2*(1-Lambda))
#c_tip = Lambda*c_root


#S = b**2/AR


print(f"Cl max used = {Cl_max:.2f} [-]")
print(f"Total area = {S:.2f} [m^2]")
print(f"Corde quarter = {c_quarter:.2f} [m] \n")
#print("Codre root:", c_root, "[m]")
#print("Codre tip:", c_tip, "[m]")

sweep_quarter = sweep_quarter*((2*np.pi)/180)
sweep_leading = np.arctan(np.tan(sweep_quarter) + (4/AR) * (((1-Lambda)/(1+Lambda)) * (0.25 - 0)))
sweep_trailing = np.arctan(np.tan(sweep_quarter) + (4/AR) * (((1-Lambda)/(1+Lambda)) * (0.25 - 1)))


y = np.linspace(0, b/2, 10)

quarter_line = (np.tan(sweep_quarter))*y + (0.25*c_quarter)
leading_edge = (np.tan(sweep_leading))*y 
trailing_edge = (np.tan(sweep_trailing))*y + c_quarter


c = leading_edge - trailing_edge

result_mac, _ = integrate.quad(lambda y: ((np.tan(sweep_trailing))*y -(np.tan(sweep_leading))*y + c_quarter)**2, 0, b/2)
result_yac, _ = integrate.quad(lambda y: ((np.tan(sweep_trailing))*y -(np.tan(sweep_leading))*y + c_quarter)*y, 0, b/2)

MAC = (2/(S))*result_mac
y_ac = (2/(S))*result_yac

print(f"MAC = {MAC:.2f} [m]")
print(f"Y_a_c = {y_ac:.2f} [m] \n")
 

beta = np.sqrt(1-(M**2))

beta_AR = beta*AR

sweep_beta = np.arctan2(np.tan(sweep_quarter), beta) * (180/(2*np.pi))
print(f"beta AR = {beta_AR:.2f}")
print(f"sweep Beta = {sweep_beta:.2f} [deg]\n")

x_ac = 0.23*MAC  # x_ac/MAC => choice done according to slide 26 Conception Aero Design

print(f"x_ac = {x_ac:.2f} [m] \n")

if wing_plot:
    plt.plot(y, leading_edge)
    plt.plot(y, trailing_edge)
    plt.plot(y, quarter_line, color='red',)
    plt.scatter(y_ac, (np.tan(sweep_leading))*y_ac + x_ac, marker="x")

    plt.xlabel('$Y$')
    plt.ylabel('$X$')
    # Fixer l'échelle des axes
    plt.axis('equal')
    plt.show() 




# ------------------- Schuemann wing planform ------------------- #
print("# ------------------- Schuemann wing planform ------------------- # \n")

c0 = 20 #[m]
c1 = 5 #[m]
c2 = 1 #[m]

b = 40


c = [c0, c1, c2]
y = [10, 10]
Lambda = c[-1]/c[0]

S = np.zeros(2)
S_tot = 0
for i in range(len(S)):
    S[i] = y[i] * (c[i] + c[i+1])/2
    S_tot += S[i]

Sw = 0
for i in range(len(S)):
    Sw += (b/S_tot) * ((c[i]*S[i]) + (c[i+1]*S[i]))


Cwre = 0
for i in range(len(S)):
    Cwre += (2/Sw) * (c[i]*S[i])



cr = Cwre / c[0]

MCG = (2*cr/3)*((1 + Lambda + Lambda**2)/(1 + Lambda))

Yac = (b/6) * ((1+ (2*Lambda))/(1+Lambda))

print(f"MCG = {MCG} [m]")
print(f"Yac = {Yac:.2f} [m]")
print(f"Cwre = {Cwre:.2f} [m]")
print(f"Sw = {Sw:.2f} [m^2]")
print(f"cr = {cr:.2f} [-]")

#help of prof.
# a (CL_alpha) 
#CL = CL_alpha * (alpha - alpha_L0)

