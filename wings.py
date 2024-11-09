import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import trapz

#---Guesses---#
span_max = 36              #[m] Span  max span for airport
cabin_width = 8    #[m] /!\ guess
cabin_lenght = 17   #[m] /!\ guess
AR = 5              #Aspect ratio (guess)
S = 0               #total area
weight = 471511.49122 #[N] = 106000lb (guess from weight code)
alti = 12500        #[m]
rho = 0.288         #[kg/m^3]
M = 0.9            #[-] Mach number
R = 287             #[m^2/s^2K]
gamma = 1.4
e = 0.85            #Ostxald's efficiency factor
CD0 = 0.02          # Zero lift drag coeff
#sweep_quarter = 15         #[°] sweep angle
Lambda = 0.6        # [-] taper ratio

#---Commande---#
polar_Cl_Cd = False
wing_plot = True
cl_plot = True


#---Code---#


#cl_alpha = np.genfromtxt('SC2-0410.csv', delimiter=',')   # SC 0410
#cl_alpha[0, 0] =-19.5
#print(cl_alpha)


def guess_CL_max():
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

    return Cl_max

Cl_max = guess_CL_max()

T = 288.15 - (6.5e-3*alti)

a = np.sqrt(gamma * R * T)  #[m^2] speed of sound
v = M*a #[m/s] 

S_total = 2*weight/(rho*(v**2)*Cl_max)

surface_fuslage = (cabin_width * cabin_lenght)/2
surface_wing = S_total - surface_fuslage

beta = np.sqrt(1-(M**2))

print(f"Cl max used = {Cl_max:.2f} [-]")
print(f"Total area = {S_total:.2f} [m^2]\n")


def wing_side():
    """
    The choice of the airfoil must be made (cl_alpha and alpha_0)
    """
    b = span_max - cabin_width - 2
    AR_wing = (b**2)/surface_wing
    
    """
    cl_alpha = 1.621 # SC(2)-0012 M0.85 Re6M
    alpha_l0 = 0
    CD_wing = 999 # /!\
    """
    cl_alpha = (1.7-0)/(4+12) * (180/np.pi) # NACA 64209 M0.85 Re6M
    alpha_l0 = -12*(np.pi/180)
    CD_wing = 0.007 
    
    sweep_quarter = 20 #[°]
    taper_ratio = 0.2

    c_root = (surface_wing/b)* (2/(1+taper_ratio))
    c_tip = taper_ratio*c_root

    sweep_quarter = sweep_quarter*((2*np.pi)/180)
    sweep_leading = np.arctan(np.tan(sweep_quarter) + (4/AR_wing) * (((1-taper_ratio)/(1+taper_ratio)) * (0.25 - 0)))
    sweep_trailing = np.arctan(np.tan(sweep_quarter) + (4/AR_wing) * (((1-taper_ratio)/(1+taper_ratio)) * (0.25 - 1)))

    sweep_beta = np.arctan2(np.tan(sweep_quarter), beta) 

    y = np.linspace(0, b/2, 10)
    quarter_line = np.zeros(len(y))
    leading_edge = np.zeros(len(y))
    trailing_edge = np.zeros(len(y))
    for i in range(len(y)):
        quarter_line[i] = (np.tan(sweep_quarter))*y[i] + (0.25*c_root)
        leading_edge[i] = (np.tan(sweep_leading))*y[i] 
        trailing_edge[i] = (np.tan(sweep_trailing))*y[i] + c_root

    c = trailing_edge -leading_edge
    wing_surface_integrate = trapz(c, y)*2 #numerical integration via method of trapeze
    print(f"Surface of wing = {surface_wing:.2f}")
    print(f"Surface of wing via integration = {wing_surface_integrate:.2f} \n")
    if wing_plot:
        plt.plot(y, leading_edge)
        plt.plot(y, trailing_edge, color='green')
        plt.plot(y, quarter_line, color='red',)

        plt.xlabel('$Y$')
        plt.ylabel('$X$')
        # Fixer l'échelle des axes
        plt.axis('equal')
        plt.show() 


    #twist
    twist_angle = 3 * (np.pi/180) #[°]
    alpha_01 = -0.17
    alpha_L0 = alpha_l0 + (alpha_01 * twist_angle)
    # Lift
    AoA = np.linspace(-10, 10, 20) * ((np.pi)/180)
    CL_w = np.zeros(len(AoA))
     
    k = (beta * cl_alpha)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR_wing)) + np.sqrt((1/((k * np.cos(sweep_beta)))**2) + ((2/(beta * AR_wing))**2) )))/beta

    for i in range(len(AoA)):    
        CL_w[i] = a*(AoA[i] - alpha_L0)
        if AoA[i] <= 0:
            if AoA[i+1] >= 0:
                CL_w0 = (CL_w[i] + a*(AoA[i+1] - alpha_L0))/2

    if cl_plot:
        plt.plot(AoA*(180/(np.pi)), CL_w)
        #plt.scatter(CL_max, CL_CD_max, marker="x", color="r")
        plt.xlabel('$AoA$')
        plt.ylabel('$Cl_w$')
        plt.show()
    

    print(f"CL_0W wing = {CL_w0:.3f}")
    return c_root, c_tip, wing_surface_integrate, CL_w0, CD_wing

c_root_wing, c_tip_wing, surface_wing, CL_wing, CD_wing = wing_side()


def wing_fuselage():
    """
    The choice of the airfoil must be made (cl_alpha and alpha_0)
    """
    b = cabin_width + 2
    
    AR_fuslage = (b**2)/surface_fuslage
    
    """
    cl_alpha = (0.82-0.2)/(1.5+2) # SC(2)-0714 M0.75 Re6M
    alpha_L0 = -3.5 * (np.pi/180) #[rad]
    CD_fuselage = 0.01 
    """
    cl_alpha = ((0.8+0.2)/(5+5)) *(180/np.pi)# Eppler 642 M0 Re1M C_m = -0.05
    alpha_L0 = -4 * (np.pi/180) #[rad]
    CD_fuselage = 0.01 
    
    sweep_quarter = 30 #[°]
    taper_ratio = 0.1
    c_root = cabin_lenght # (surface_wing/b)* (2/(1+taper_ratio))
    
    taper_ratio =  c_root_wing/c_root
    
    sweep_quarter = sweep_quarter*((2*np.pi)/180)
    sweep_leading = np.arctan(np.tan(sweep_quarter) + (4/AR_fuslage) * (((1-taper_ratio)/(1+taper_ratio)) * (0.25 - 0)))
    sweep_trailing = np.arctan(np.tan(sweep_quarter) + (4/AR_fuslage) * (((1-taper_ratio)/(1+taper_ratio)) * (0.25 - 1)))

    sweep_beta = np.arctan2(np.tan(sweep_quarter), beta) 

    y = np.linspace(0, b/2, 10)
    quarter_line = np.zeros(len(y))
    leading_edge = np.zeros(len(y))
    trailing_edge = np.zeros(len(y))
    for i in range(len(y)):
        quarter_line[i] = (np.tan(sweep_quarter))*y[i] + (0.25*c_root)
        leading_edge[i] = (np.tan(sweep_leading))*y[i] 
        trailing_edge[i] = (np.tan(sweep_trailing))*y[i] + c_root

    c = trailing_edge - leading_edge
    fuselage_surface_integrate = trapz(c, y)*2 #numerical integration via method of trapeze
    print(f"Surface of wing = {surface_fuslage:.2f}")
    print(f"Surface of wing via integration = {fuselage_surface_integrate:.2f} \n")

    if wing_plot:
        plt.plot(y, leading_edge)
        plt.plot(y, trailing_edge, color='green')
        plt.plot(y, quarter_line, color='red',)
        
        plt.xlabel('$Y$')
        plt.ylabel('$X$')
        # Fixer l'échelle des axes
        plt.axis('equal')
        plt.show() 


    #Lift
    AoA = np.linspace(-10, 10, 50) * ((np.pi)/180)
    CL_w = np.zeros(len(AoA))
    
     
    k = (beta * cl_alpha)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR_fuslage)) + np.sqrt((1/((k * np.cos(sweep_beta)))**2) + ((2/(beta * AR_fuslage))**2) )))/beta

    for i in range(len(AoA)):    
        CL_w[i] = a*(AoA[i] - alpha_L0)
        if AoA[i] <= 0:
            if AoA[i+1] >= 0:
                CL_w0 = (CL_w[i] + a*(AoA[i+1] - alpha_L0))/2
    
    print(f"CL_w0 fuselage = {CL_w0:.3f}")
    
    if cl_plot:
        plt.plot(AoA*(180/(np.pi)), CL_w)
        #plt.scatter(CL_max, CL_CD_max, marker="x", color="r")
        plt.xlabel('$AoA$')
        plt.ylabel('$Cl_w$')
        plt.title("Lift fuselage")
        plt.show()
    
    return fuselage_surface_integrate, CL_w0, CD_fuselage


surface_fuslage, CL_fuselage, CD_fuselage = wing_fuselage()

surface_total = surface_wing + surface_fuslage

CL = ((CL_fuselage * surface_fuslage) + (CL_wing * surface_wing))/surface_total

print(f"Total lift recomputed = {CL:.3}")

#print(c_tip_wing/cabin_lenght) # 0.027813732207197624
delta = 0.005 #graph slide 61 lecture 6 aerodinimics

CD_induce = ((CL**2)/(np.pi* AR))*(1+delta)

CD = CD_induce + CD_wing + CD_fuselage
#CD_gest = ((Cl_max**2)/(np.pi* AR))*(1+delta) + CD_wing + CD_fuselage
print(f"CD = {CD:.3f}")







"""
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
 



beta_AR = beta*AR

sweep_beta = np.arctan2(np.tan(sweep_quarter), beta) 
print(f"beta AR = {beta_AR:.2f}")
print(f"sweep Beta = {(sweep_beta*(180/(2*np.pi))):.2f} [deg]\n")

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


# maximum thickness at s/2

M_star = 1.15 - (Cl_max/(4*(np.cos(sweep_quarter))**2))
t_over_c = (3/(10*M)) * np.cbrt((1/(M*np.cos(sweep_quarter))) - (M*np.cos(sweep_quarter))) * (1- ((5+ M**2 * (np.cos(sweep_quarter)**2))/(5 + M_star**2))**3.5)**(2/3)

print(f"t_over_c = {t_over_c:.2f} [m] \n")




# ------ Lift ------ #

#sweep_beta = np.arctan(np.tan(sweep_quarter)/beta)

cl_alpha = 1.621 # SC(2)-0012 M0.85 Re6M

AoA = np.linspace(-10, 10, 20) * ((np.pi)/180)
CL_w = np.zeros(len(AoA))
alpha_L0 = 0 

k = (beta * cl_alpha)/(2*np.pi)
a = ((2*np.pi)/((2/(beta*AR)) + np.sqrt((1/((k * np.cos(sweep_beta)))**2) + ((2/(beta * AR))**2) )))/beta

for i in range(len(AoA)):    
    CL_w[i] = a*(AoA[i] - alpha_L0)

if cl_plot:
    plt.plot(AoA*(180/(np.pi)), CL_w)
    #plt.scatter(CL_max, CL_CD_max, marker="x", color="r")
    plt.xlabel('$AoA$')
    plt.ylabel('$Cl_w$')
    plt.show()




#page 345 reference book (SNORRI)
#have to find \tau


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

"""











# ------------ Panknin and Culver twist formulas ------------ #
# read page 324 of book of reference (SNORRI)
# to do that we have to compute the lift