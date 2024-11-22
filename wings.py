import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import trapz

#---Guesses---#
span_max = 29           #[m] Span  max span for airport
cabin_width = 7         #[m] 
cabin_lenght = 16.8     #[m] 
AR = 1.8                #Aspect ratio (guess)
AoA_wing = 6            #[°]
weight = 471511.49122   #[N] = 106000lb (guess from weight code)
alti = 12500            #[m]
M = 0.9                 #[-] Mach number
R = 287                 #[m^2/s^2K]
gamma = 1.4
e = 0.85                #Ostxald's efficiency factor
CD0 = 0.02              # Zero lift drag coeff
delta = 0.005 #graph slide 61 lecture 6 aerodinimics
#sweep_quarter = 15     #[°] sweep angle
#Lambda = 0.6           # [-] taper ratio


#---Commande---#
polar_Cl_Cd = False
wing_plot = False
cl_plot = False


# --- globale constant --- #
T = 216.5#288.15 - (6.5e-3*alti)
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

surface_fuselage = (cabin_lenght + ((cabin_width/2 +1)/np.tan(30*np.pi/180)))/2 * (cabin_width+2) #(cabin_width * cabin_lenght)
surface_wing = surface_total - surface_fuselage

beta = np.sqrt(1-(M**2))


def wingGeometry():
    b = span_max - cabin_width - 2
    AR_wing = (b**2)/surface_wing
    sweep_leading = 40 #[°]
    c_tip = 0.8 #[m]
    c_root = (2*surface_wing/b) - c_tip
    taper_ratio = c_tip/c_root
    sweep_leading = sweep_leading*((np.pi)/180)
    sweep_quarter = np.arctan(np.tan(sweep_leading) + (4/AR_wing) * (((1-taper_ratio)/(1+taper_ratio)) * (0 - 0.25)))
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

    #geometrical_twist = twist_angle*((taper_ratio * y_over_S)/(1-(1-taper_ratio)*y_over_S ))

    return b, AR_wing, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line

def wingPlot(wing_plot):
    if wing_plot == False:
        return
    _, _, _, _, _, _, _, y, leading_edge, trailing_edge, quarter_line = wingGeometry() 
    plt.plot(y, leading_edge)
    plt.plot(y, trailing_edge, color='green')
    plt.plot(y, quarter_line, color='red',)

    plt.xlabel('$Y$')
    plt.ylabel('$X$')
    # Fixer l'échelle des axes
    plt.axis('equal')
    plt.show() 

    return
wingPlot(wing_plot)

def wingCL():
    # ----- Airfoil ----- #
    cl_alpha = (1.4073-0.269)/(5+5) * (180/np.pi) # SC(2)-1010 M0.85 Re12M
    cl_max = 1.4
    alpha_l0 = -7*(np.pi/180)
    CD_wing = 0.01091 #at 6° aoa  and at 0° aoa = 0.006
    """
    cl_alpha = (0.5765-0)/(5+0) * (180/np.pi) # SC(2)-0012 M0.85 Re12M
    alpha_l0 = 0
    CD_wing = 0.0059 
    
    cl_alpha = (0.8293+0.3558)/(5+5) * (180/np.pi) # SC(2)-0410 M0.85 Re12M
    alpha_l0 = -2*(np.pi/180)
    CD_wing = 0.006

    cl_alpha = (1.7-0)/(4+12) * (180/np.pi) # NACA 64209 M0.85 Re6M
    alpha_l0 = -12*(np.pi/180)
    CD_wing = 0.007 
    """
    b, AR_wing, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, _, _, _, _ = wingGeometry()

    # --- Twist angle --- #
    twist_angle = 3 * (np.pi/180) #[°]
    alpha_01 = -0.17
    eta_a_tip = twist_angle 
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

    CL_w0 = a*((AoA_wing * (np.pi/180)) - alpha_L0) # choose of AoA of the wing 

    CL_max = np.cos(sweep_quarter) * 0.95 * ((cl_max + cl_max)/2)

    if cl_plot:
        plt.plot(AoA*(180/(np.pi)), CL_w)
        #plt.scatter(CL_max, CL_CD_max, marker="x", color="r")
        plt.xlabel('$AoA$')
        plt.ylabel('$Cl_w$')
        plt.show()
    
    return CL_w0,  CD_wing, CL_max, alpha_L0

def fusGeometry():
    b = cabin_width + 2
    AR_fuselage = (b**2)/surface_fuselage
    sweep_leading =  60 #[°]
    c_root = cabin_lenght # (surface_wing/b)* (2/(1+taper_ratio))
    c_tip = ((cabin_width/2 +1)/np.tan(30*np.pi/180)) #((2*surface_fuselage)/b) - c_root
    taper_ratio =  c_tip/c_root
    sweep_leading = sweep_leading*((np.pi)/180)
    sweep_quarter = np.arctan(np.tan(sweep_leading) + (4/AR_fuselage) * (((1-taper_ratio)/(1+taper_ratio)) * (0 - 0.25)))
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
    return b, AR_fuselage, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line

def fusPlot(wing_plot):
    if wing_plot == False:
        return
    _, _, _, _, _, _, _, y, leading_edge, trailing_edge, quarter_line = fusGeometry() 
    plt.plot(y, leading_edge)
    plt.plot(y, trailing_edge, color='green')
    plt.plot(y, quarter_line, color='red',)

    plt.xlabel('$Y$')
    plt.ylabel('$X$')
    # Fixer l'échelle des axes
    plt.axis('equal')
    plt.show() 
    return
fusPlot(wing_plot)

def fuselageCL():
    # --- airfoil --- #
    cl_alpha = ((1.0498+0.2062)/(5+5)) * (180/np.pi) # SC(2) 0518 M0 Re12M C_m = -0.1158
    cl_max = 1.87
    alpha_L0 = -3.5 * (np.pi/180) #[rad] 
    CD_fuselage = 0.00636 
    """
    cl_alpha = (0.82-0.2)/(1.5+2) # SC(2)-0714 M0.75 Re6M
    alpha_L0 = -3.5 * (np.pi/180) #[rad]
    CD_fuselage = 0.01 

    cl_alpha = ((0.8+0.2)/(5+5)) *(180/np.pi) # Eppler 642 M0 Re1M C_m = -0.058
    alpha_L0 = -4 * (np.pi/180) #[rad] 
    CD_fuselage = 0.01 
    """
    b, AR_fuselage, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, _, _, _, _ = fusGeometry()

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
    
    if cl_plot:
        plt.plot(AoA*(180/(np.pi)), CL_w)
        #plt.scatter(CL_max, CL_CD_max, marker="x", color="r")
        plt.xlabel('$AoA$')
        plt.ylabel('$Cl_w$')
        plt.title("Lift fuselage")
        plt.show()
    
    CL_max = np.cos(sweep_quarter) * 0.95 * ((cl_max + cl_max)/2)
    
    return CL_w0, CD_fuselage, CL_max

def getMAC():
    _, _, _, _, _, _, _, y_wing, leading_wing, trailing_wing, quarter_wing = wingGeometry() 
    _, _, _, _, _, _, _, y_fus, leading_fus, trailing_fus, quarter_fus = fusGeometry() 

    # -- fuselage -- #
    c_fus = trailing_fus - leading_fus
    MAC_fus = (2/surface_total) * trapz(c_fus*2, y_fus) #numerical integration via method of trapez
    cy_fus = c_fus*y_fus
    yac_fus = (2/surface_total) * trapz(cy_fus, y_fus)
    xac_fus = MAC_fus*0.2

    # -- wing -- #
    c_wing = trailing_wing - leading_wing
    MAC_wing = (2/surface_total) * trapz(c_wing*2, y_wing) #numerical integration via method of trapez
    cy_wing = c_wing*y_wing
    yac_wing = (2/surface_total) * trapz(cy_wing, y_wing)
    xac_wing = MAC_wing*0.2

    # -- total -- #
    c = np.concatenate((c_fus, c_wing))
    y = np.linspace(0, span_max/2, 20)
    MAC = (2/surface_total) * trapz(c*2, y)
    cy = c*y
    yac = (2/surface_total) * trapz(cy, y)
    xac = 0.2*MAC

    return MAC_fus, yac_fus, xac_fus, MAC_wing, yac_wing, xac_wing, MAC, yac, xac


def get_Lift_and_drag(AR, delta):

    Cl_wing, Cd_wing, Cl_max_wing, _ = wingCL()
    Cl_fuselage, Cd_fuselage, Cl_max_fus = fuselageCL()

    # --- total lift computation --- #
    Cl_tot = ((Cl_wing*surface_wing) + (Cl_fuselage*surface_fuselage))/surface_total 
    Cl_max = ((Cl_max_wing*surface_wing) + (Cl_max_fus*surface_fuselage))/surface_total 
    # --- total drag computation --- #
    Cd_induce = ((Cl_tot**2)/(np.pi* AR)) * (1+delta)
    Cd0_tot = (((Cd_wing*surface_wing) + (Cd_fuselage*surface_fuselage))/surface_total)
    Cd_tot = Cd_induce + (((Cd_wing*surface_wing) + (Cd_fuselage*surface_fuselage))/surface_total)
    
    return Cl_tot, Cd_tot, Cl_max


def wingMaxthickness():
    CL, _, _, _ = wingCL()
    b, AR_wing, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = wingGeometry()
    M_star = 1.15-(CL/(4*(np.cos(sweep_quarter)**2)))
    t_bar_over_c =  (3/(10*M)) * np.cbrt((1/(M*np.cos(sweep_quarter))) - (M*np.cos(sweep_quarter))) * (1-(((5 + (M**2)*(np.cos(sweep_quarter)**2))/(5 + (M_star**2)))**3.5))**(2/3)
    t_root = t_bar_over_c*c_root
    t_tip = t_bar_over_c*c_tip 
    
    return t_root, t_tip


def wingFuelvolume():
    b, AR_wing, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = wingGeometry()
    tau = 1
    V_fuel = 0.54*((surface_wing**2)/b) * 0.10 * ((1+(taper_ratio*np.sqrt(tau))+((taper_ratio**2)*tau))/(1+taper_ratio)**2)
    return V_fuel


def wingSurfaceWetted():
    b, AR_wing, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = wingGeometry()
    S_wet = 2*surface_wing*(1+ (1/4) * ((0.10 + (0.1*taper_ratio))/(1+taper_ratio)))
    return S_wet

def stallVelocity():
    b, AR_wing, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = wingGeometry()
    _, _, Cl_max = get_Lift_and_drag(AR, delta)

    rho_sl = 1.225

    Vs = np.sqrt((weight/surface_total) * (2/rho_sl) * (1/(1.133*Cl_max)))
    Cl_max0 = 2 * np.cos(sweep_quarter)
    W0 = 555 #landing weight
    Vs0 = np.sqrt((W0/surface_total) * (2/rho_sl) * (1/(1.133*Cl_max0)))
    return Vs, Vs0



def printFunction():

    print(f"Cl max used = {Cl_max:.2f} [-]\n")
    print(f"Total area = {surface_total:.2f} [m^2]")
    print(f"Surface of fuselage = {surface_fuselage:.2f}")
    print(f"Surface of wing = {surface_wing:.2f} \n")
    print(f"Compressibility parameter: {beta:.3f}")
    
    b, AR_wing, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = wingGeometry()
    Cl_wing, Cd_wing, Cl_max_wing, alpha_L0 = wingCL()
    print("\n-------------- wing values --------------\n")
    print(f"\nAR wing: {AR_wing:.3f} \nCL_w0 wing = {Cl_wing:.3f} \n")
    print(f"Cord at wing root: {c_root:.3f} \nCorde at wing tip: {c_tip:.3f}")
    print(f"Taper ratio: {taper_ratio:.3f}")
    print(f"sweep quater: {sweep_quarter*(180/np.pi):.3f}")
    print(f"Wing lift coefficient derivative: {a:.3f}")
    print(f"Alpha_L0: {alpha_L0*(180/np.pi):.3f}")

    b, AR_fuselage, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = fusGeometry()
    Cl_fuselage, Cd_fuselage, Cl_max_fus = fuselageCL()
    print("\n-------------- fuselage values --------------\n")
    print(f"\nAR fuselage: {AR_fuselage:.3f} \nCL_w0 fuselage = {Cl_fuselage:.3f} \n")
    print(f"Cord at fuselage root: {c_root:.3f} \nCorde at fuselage tip: {c_tip:.3f}")
    print(f"Taper ratio: {taper_ratio:.3f}")
    print(f"sweep quater: {sweep_quarter*(180/np.pi):.3f}")
    print(f"Fuselage lift coefficient derivative: {a:.3f}\n")

    MAC_fus, yac_fus, xac_fus, MAC_wing, yac_wing, xac_wing, MAC, yac, xac = getMAC()
    print(f"MAC fus: {MAC_fus:.3f} \nYac fus: {yac_fus:.3f} \nXac fus: {xac_fus:.3f} \n")
    print(f"MAC wing: {MAC_wing:.3f} \nYac wing: {yac_wing:.3f} \nXac wing: {xac_wing:.3f} \n")
    print(f"MAC: {MAC:.3f} \nYac: {yac:.3f} \nXac: {xac:.3f} \n")

    delta = 0.005 #graph slide 61 lecture 6 aerodinimics
    lift_coef, drag_coef, CL_max = get_Lift_and_drag(AR, delta)
    print(f"\n CL = {lift_coef:.3f}[-] \n CD = {drag_coef:.3f}[-] \n")
    print(f"Cl max: {CL_max:.3f}")

    t_root, t_tip = wingMaxthickness()
    print(f"Thickness root: {t_root:.3f}, thickness tip: {t_tip:.3f}")
    print(f"Mean thinckness: {(t_root+t_tip)/2:.3f}")

    V_fuel = wingFuelvolume()
    print(f"Fuel volume in wing: {V_fuel:.3f}")

    Swetted = wingSurfaceWetted()
    print(f"Surface wetted: {Swetted:.3f}")

    Vs, Vs0 = stallVelocity()
    print(f"Stall velocity: {Vs:.3f} \nStall velocity in approach config: {Vs0:.3f}")
    return

printFunction()