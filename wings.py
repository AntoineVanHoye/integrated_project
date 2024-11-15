import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import trapz

#---Guesses---#
span_max = 36              #[m] Span  max span for airport
cabin_width = 7    #[m] /!\ guess
cabin_lenght = 12   #[m] /!\ guess
AR = 1.5              #Aspect ratio (guess)
AoA_wing = 4              #[°]
weight = 471511.49122 #[N] = 106000lb (guess from weight code)
alti = 12500        #[m]
M = 0.9            #[-] Mach number
R = 287             #[m^2/s^2K]
gamma = 1.4
e = 0.85            #Ostxald's efficiency factor
CD0 = 0.02          # Zero lift drag coeff
#sweep_quarter = 15         #[°] sweep angle
#Lambda = 0.6        # [-] taper ratio


#---Commande---#
polar_Cl_Cd = False
wing_plot = True
cl_plot = False


# --- globale constant --- #
T = 288.15 - (6.5e-3*alti)
a = np.sqrt(gamma * R * T)  #[m^2] speed of sound
v = M*a #[m/s]
rho = 0.288         #[kg/m^3]

#---Code---#
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

surface_total = 2*weight/(rho*(v**2)*Cl_max)

surface_fuselage = (cabin_width * cabin_lenght)
surface_wing = surface_total - surface_fuselage

beta = np.sqrt(1-(M**2))

print(f"Cl max used = {Cl_max:.2f} [-]\n")
print(f"Total area = {surface_total:.2f} [m^2]")
print(f"Surface of fuselage = {surface_fuselage:.2f}")
print(f"Surface of wing = {surface_wing:.2f} \n")

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

    #c = trailing_edge -leading_edge
    #wing_surface_integrate = trapz(c, y)*2 #numerical integration via method of trapeze
    
    if wing_plot:
        plt.plot(y, leading_edge)
        plt.plot(y, trailing_edge, color='green')
        plt.plot(y, quarter_line, color='red',)

        plt.xlabel('$Y$')
        plt.ylabel('$X$')
        # Fixer l'échelle des axes
        plt.axis('equal')
        plt.show() 

    # --- Twist angle --- #
    twist_angle = 3 * (np.pi/180) #[°]
    alpha_01 = -0.17
    alpha_L0 = alpha_l0 + (alpha_01 * twist_angle)

    # --- Lift --- #
    AoA = np.linspace(-10, 10, 20) * ((np.pi)/180)
    CL_w = np.zeros(len(AoA))
    
    k = (beta * cl_alpha)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR_wing)) + np.sqrt((1/((k * np.cos(sweep_beta))))**2 + ((2/(beta * AR_wing))**2) )))/beta

    for i in range(len(AoA)):    
        CL_w[i] = a*(AoA[i] - alpha_L0)
        if AoA[i] <= 0:
            if AoA[i+1] >= 0:
                CL_w0 = (CL_w[i] + a*(AoA[i+1] - alpha_L0))/2

    CL_w0 = a*((AoA_wing * (np.pi/180)) - alpha_L0) # choce of AoA of the wing 

    if cl_plot:
        plt.plot(AoA*(180/(np.pi)), CL_w)
        #plt.scatter(CL_max, CL_CD_max, marker="x", color="r")
        plt.xlabel('$AoA$')
        plt.ylabel('$Cl_w$')
        plt.show()
    
    print(f"AR wing: {AR_wing:.3f} \nCL_w0 wing = {CL_w0:.3f} \n")

    Lift_wing = 0.5*rho*(v**2)*surface_wing*CL_w0
    Drag_wing = 0.5*rho*(v**2)*surface_wing*CD_wing
    return CL_w0,  CD_wing#Lift_wing, Drag_wing


def wing_fuselage():
    """
    The choice of the airfoil must be made (cl_alpha and alpha_0)
    """
    b = cabin_width + 2
    
    AR_fuselage = (b**2)/surface_fuselage
    print(f"AR fuselage: {AR_fuselage}")
    """
    cl_alpha = (0.82-0.2)/(1.5+2) # SC(2)-0714 M0.75 Re6M
    alpha_L0 = -3.5 * (np.pi/180) #[rad]
    CD_fuselage = 0.01 
    """
    # --- airfoil --- #
    cl_alpha = ((0.8+0.2)/(5+5)) *(180/np.pi)# Eppler 642 M0 Re1M C_m = -0.05
    alpha_L0 = -4 * (np.pi/180) #[rad] 
    CD_fuselage = 0.01 
    
    sweep_quarter = 30 #[°]
    taper_ratio = 0.1
    c_root = cabin_lenght # (surface_wing/b)* (2/(1+taper_ratio))
    c_tip = ((2*surface_fuselage)/b) - c_root

    taper_ratio =  c_tip/c_root
    
    sweep_quarter = sweep_quarter*((2*np.pi)/180)
    sweep_leading = np.arctan(np.tan(sweep_quarter) + (4/AR_fuselage) * (((1-taper_ratio)/(1+taper_ratio)) * (0.25 - 0)))
    sweep_trailing = np.arctan(np.tan(sweep_quarter) + (4/AR_fuselage) * (((1-taper_ratio)/(1+taper_ratio)) * (0.25 - 1)))

    sweep_beta = np.arctan2(np.tan(sweep_quarter), beta) 

    y = np.linspace(0, b/2, 10)
    quarter_line = np.zeros(len(y))
    leading_edge = np.zeros(len(y))
    trailing_edge = np.zeros(len(y))
    for i in range(len(y)):
        quarter_line[i] = (np.tan(sweep_quarter))*y[i] + (0.25*c_root)
        leading_edge[i] = (np.tan(sweep_leading))*y[i] 
        trailing_edge[i] = (np.tan(sweep_trailing))*y[i] + c_root

    #c = trailing_edge - leading_edge
    #fuselage_surface_integrate = trapz(c, y)*2 #numerical integration via method of trapez
    if wing_plot:
        plt.plot(y, leading_edge)
        plt.plot(y, trailing_edge, color='green')
        plt.plot(y, quarter_line, color='red',)
        
        plt.xlabel('$Y$')
        plt.ylabel('$X$')
        # Fixer l'échelle des axes
        plt.axis('equal')
        plt.show() 

    # --- Lift --- #
    AoA = np.linspace(-10, 10, 50) * ((np.pi)/180)
    CL_w = np.zeros(len(AoA))
    
    k = (beta * cl_alpha)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR_fuselage)) + np.sqrt( ((1/(k * np.cos(sweep_beta)))**2) + ((2/(beta * AR_fuselage))**2) )))/beta

    for i in range(len(AoA)):    
        CL_w[i] = a*(AoA[i] - alpha_L0)
        if AoA[i] <= 0:
            if AoA[i+1] >= 0:
                CL_w0 = (CL_w[i] + a*(AoA[i+1] - alpha_L0))/2
    
    print(f"CL_w0 fuselage = {CL_w0:.3f}\n")
    
    if cl_plot:
        plt.plot(AoA*(180/(np.pi)), CL_w)
        #plt.scatter(CL_max, CL_CD_max, marker="x", color="r")
        plt.xlabel('$AoA$')
        plt.ylabel('$Cl_w$')
        plt.title("Lift fuselage")
        plt.show()
    
    Lift_fuselage = 0.5*rho*(v**2)*surface_fuselage*CL_w0
    Drag_fuselage = 0.5*rho*(v**2)*surface_fuselage*CD_fuselage
    return CL_w0, CD_fuselage #Lift_fuselage, Drag_fuselage

"""
CL = ((CL_fuselage * (surface_fuselage)) + (CL_wing * (surface_wing)))/(surface_total)

print(f"Total lift recomputed = {CL:.3}")

#print(c_tip_wing/cabin_lenght) # 0.027813732207197624
delta = 0.005 #graph slide 61 lecture 6 aerodinimics

CD_induce = ((CL**2)/(np.pi* AR)) * (1+delta)

CD = CD_induce + (((CD_fuselage * surface_fuselage) + (CD_wing * surface_wing))/surface_total)
#CD_gest = ((Cl_max**2)/(np.pi* AR))*(1+delta) + CD_wing + CD_fuselage
print(f"CD = {CD:.3f}")


L_tot = 0.5*rho*(v**2)*surface_total*CL

L_wing = 0.5*rho*(v**2)*surface_wing*CL_wing
L_fuselage = 0.5*rho*(v**2)*surface_fuselage*CL_fuselage
L_tot_add = ((L_wing*surface_wing) + (L_fuselage * surface_fuselage))/surface_total
CL = L_tot_add/(0.5*rho*(v**2)*surface_total)

print(f"Lift total method 1: {L_tot:.3f}")
print(f"Lift total method 2: {L_tot_add:.3f}")
print(CL)

"""

def get_Lift_and_drag(AR, delta):
    #lift_wing, drag_wing = wing_side()
    #lift_fuselage, drag_fuselage = wing_fuselage()

    Cl_wing, Cd_wing = wing_side()
    Cl_fuselage, Cd_fuselage = wing_fuselage()

    """
    # --- Adimentionnalisation --- #
    Cl_wing = lift_wing/(0.5*rho*(v**2)*surface_total)
    Cd_wing = drag_wing/surface_total

    Cl_fuselage = lift_fuselage/(0.5*rho*(v**2)*surface_total)
    Cd_fuselage = drag_fuselage/(0.5*rho*(v**2)*surface_total)
    """
    # --- total lift computation --- #
    Cl_tot = ((Cl_wing*surface_wing) + (Cl_fuselage*surface_fuselage))/surface_total 

    # --- total drag computation --- #
    Cd_induce = ((Cl_tot**2)/(np.pi* AR)) * (1+delta)
    Cd_tot = Cd_induce + (((Cd_wing*surface_wing) + (Cd_fuselage*surface_fuselage))/surface_total)
    
    return Cl_tot, Cd_tot

delta = 0.005 #graph slide 61 lecture 6 aerodinimics
lift_coef, drag_coef = get_Lift_and_drag(AR, delta)

print(f"\n CL = {lift_coef:.3f} \n CD = {drag_coef:.3f} \n")




# ------------ Panknin and Culver twist formulas ------------ #
# read page 324 of book of reference (SNORRI)
# to do that we have to compute the lift