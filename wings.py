import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import trapz

#---Guesses---#
span_max = 29           #[m] Span  max span for airport
cabin_width = 7         #[m] 
cabin_lenght = 16.8     #[m] 
AR = 1.5                #Aspect ratio (guess)
weight = 526898.738020239# 471511.49122 #  #[N] = 106000lb (guess from weight code)
weight_empty = 253488.33 #60452.314059821154 * 9.81 #[N] 
alti = 12500            #[m]
M = 0.9                #[-] Mach number
R = 287                 #[m^2/s^2K]
gamma = 1.4
e = 0.85                #Ostxald's efficiency factor
CD0 = 0.02              # Zero lift drag coeff
delta = 0.005 #graph slide 61 lecture 6 aerodinimics
#sweep_quarter = 15     #[°] sweep angle
#Lambda = 0.6           # [-] taper ratio


#---Commande---#
polar_Cl_Cd = False
wing_plot = True
cl_plot = False
lift_and_drag_plots = False

#---Code---#

def winglet():
    #formule from Snoris
    h = 1 #[m] Choice
    delta_AR = 1.9*(h/span_max) * AR
    return delta_AR
AR = AR + winglet()

# Function to calculate air density using the ISA model
def air_density(altitude):
    # Up to 11 km (Troposphere)
    if altitude <= 11000:
        T = 288.15 - 0.0065 * altitude  # Temperature [K]
        p = 101325 * (T / 288.15) ** 5.2561  # Pressure [Pa]
    else:
        # Simplification for stratosphere, constant T [K] above 11 km
        T = 216.65  # Constant temperature [K]
        p = 22632 * np.exp(-9.81 * (altitude - 11000) / (287.05 * T))
    rho = p / (287.05 * T)  # Air density [kg/m^3]
    return rho, T
rho, T = air_density(alti)

# Function to calculate the true airspeed at a given altitude
def true_airspeed_at_altitude(altitude):
    T = air_density(altitude)[1]
    a = np.sqrt(gamma * R * T)  # [m/s] Speed of sound
    v = M * a  # [m/s] Aircraft velocity
    return v
v = true_airspeed_at_altitude(alti)

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

c_root_fus = cabin_lenght - ((cabin_width/2 +1)/np.tan(30*np.pi/180))
surface_fuselage = (cabin_lenght +(cabin_lenght - ((cabin_width/2 +1)/np.tan(30*np.pi/180))))/2 * (cabin_width+2) #(cabin_width * cabin_lenght)
surface_wing = surface_total - surface_fuselage

beta = np.sqrt(1-(M**2))


def fusGeometry():
    b = cabin_width + 2
    AR_fuselage = (b**2)/surface_fuselage
    sweep_leading =  60 #[°]
    c_root = cabin_lenght # (surface_wing/b)* (2/(1+taper_ratio))
    c_tip = cabin_lenght - ((b/2)/np.tan((90-sweep_leading)*np.pi/180)) #((2*surface_fuselage)/b) - c_root
    taper_ratio =  c_tip/c_root
    sweep_leading = sweep_leading*((np.pi)/180)
    sweep_quarter = np.arctan(np.tan(sweep_leading) + (4/AR_fuselage) * (((1-taper_ratio)/(1+taper_ratio)) * (0 - 0.25)))
    sweep_trailing = 0 # np.arctan(np.tan(sweep_quarter) + (4/AR_fuselage) * (((1-taper_ratio)/(1+taper_ratio)) * (0.25 - 1)))
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
    AoA = np.linspace(-10, 10, 51) * ((np.pi)/180)
    CL_w = np.zeros(len(AoA))
    
    k = (beta * cl_alpha)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR_fuselage)) + np.sqrt( ((1/(k * np.cos(sweep_beta)))**2) + ((2/(beta * AR_fuselage))**2) )))/beta

    for i in range(len(AoA)):    
        CL_w[i] = a*(AoA[i] - alpha_L0)
        if AoA[i] <= 0:
            if AoA[i+1] >= 0:
                CL_w0 = (CL_w[i] + a*(AoA[i+1] - alpha_L0))/2
    CL_w0 = a*(0 - alpha_L0)

    if cl_plot:
        plt.plot(AoA*(180/(np.pi)), CL_w)
        #plt.scatter(CL_max, CL_CD_max, marker="x", color="r")
        plt.xlabel('$AoA$')
        plt.ylabel('$Cl_w$')
        plt.title("Lift fuselage")
        plt.show()
    
    CL_max = np.cos(sweep_quarter) * 0.95 * ((cl_max + cl_max)/2)
    
    return CL_w, CL_w0, CD_fuselage, CL_max, a



def wingGeometry():
    b_fus, AR_fuselage, sweep_beta_fus, c_root_fus, taper_ratio_fus, sweep_quarter_fus, c_tip_fus, y_fus, leading_edge_fus, trailing_edge_fus, quarter_line_fus = fusGeometry()

    b = span_max - cabin_width - 2
    AR_wing = (b**2)/surface_wing
    sweep_leading = 40 #[°]
    """
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
    """
    
    # ---- Complex wing ---- #
    """
    h = [0.5, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, (b*0.5) - 3] #[m]
    Yposition = [0, h[0], h[0]+h[1], b/2]
    c = [c_tip_fus, 8.5, 7.8, 7.3, 6.3, 5.7, 5.3, 5.0 , 0] #[c_tip_fus, 7, 6.3, 5.8, 4.8, 4.2, 3.8, 3.5 , 0]
    """

    h = [0.5, 1.0, 1.5, 2.0, (b*0.5) - 3.5] #[m]
    Yposition = [0, h[0], h[0]+h[1], b/2]
    c = [c_tip_fus, 6.751810808525195, 5.578178429906, 4.6202933997553615, 3.7940656900653904, 0] #[c_tip_fus, 7, 6.3, 5.8, 4.8, 4.2, 3.8, 3.5 , 0]
    
    S = np.zeros(len(h))
    S_sum =0
    for i in range(len(c)-2):
        S[i] = (c[i]+ c[i+1])*0.5*h[i]
        S_sum += S[i]
    #S1 = (c[0]+ c[1])*0.5*h[0]
    #S2 = (c[1] + c[2])*0.5*h[1]
    S[-1] = (surface_wing*0.5) - S_sum
    c[-1] = (2/h[-1]) * S[-1] - c[-2] #c_tip

    sweep_leading = 40 #[°]
    sweep_leading = sweep_leading*((np.pi)/180)
    sweep_quarter = np.zeros(len(c)-1)
    sweep_trailing = np.zeros(len(c)-1)
    sweep_beta = np.zeros(len(c)-1)

    for i in range(len(c)-1):
        sweep_quarter[i] = np.arctan2(h[i], ((c[i]*0.25) - ((h[i]*np.tan(sweep_leading)) + (c[i+1]*0.25)))) - (np.pi/2)#np.arctan(np.tan(sweep_leading) + (4/AR_tmp) * (((1-taper_tmp)/(1+taper_tmp)) * (0 - 0.25)))
        sweep_trailing[i] = np.arctan2(h[i], (c[i] - ((h[i]*np.tan(sweep_leading)) + c[i+1]))) - (np.pi/2) #np.arctan(np.tan(sweep_quarter[i]) + (4/AR_tmp) * (((1-taper_tmp)/(1+taper_tmp)) * (0.25 - 1)))
        sweep_beta[i] = np.arctan2(np.tan(sweep_quarter[i]), beta)

    y = np.array([])
    quarter_line = np.array([])
    leading_edge = np.array([])
    trailing_edge = np.array([])

    quarter_base = 0
    leading_base = 0
    trailing_base = 0
    y_base = 0
    for j in range(len(c)-1):
        y_tmp = np.linspace(0, h[j], 10)
        quarter_line_tmp = np.zeros(len(y_tmp))
        leading_edge_tmp = np.zeros(len(y_tmp))
        trailing_edge_tmp = np.zeros(len(y_tmp))
        
        for i in range(len(y_tmp)):
            leading_edge_tmp[i] =  (np.tan(sweep_leading))*y_tmp[i] + leading_base
            quarter_line_tmp[i] = (np.tan(sweep_quarter[j]))*y_tmp[i] + (0.25*c[j]) + leading_base#+ quarter_base
            trailing_edge_tmp[i] = (np.tan(sweep_trailing[j]))*y_tmp[i] + c[j] + leading_base#+ trailing_base
        y_tmp = y_tmp + y_base
        quarter_base = quarter_line_tmp[-1]
        
        leading_base = leading_edge_tmp[-1]
        trailing_base = trailing_edge_tmp[-1]
        y_base = y_tmp[-1]
        y = np.concatenate((y, y_tmp))
        quarter_line = np.concatenate((quarter_line, quarter_line_tmp))
        leading_edge = np.concatenate((leading_edge, leading_edge_tmp))
        trailing_edge = np.concatenate((trailing_edge, trailing_edge_tmp))
    
    sweep_beta_tot = 0
    for i in range(len(sweep_beta)):
        sweep_beta_tot += sweep_beta[i] * h[i]
    sweep_beta_tot = sweep_beta_tot/(b/2)

    taper_ratio = c[-1]/c[0]
    #geometrical_twist = twist_angle*((taper_ratio * y_over_S)/(1-(1-taper_ratio)*y_over_S ))
    sweep_quarter_tot = 0
    for i in range(len(sweep_quarter)):
        sweep_quarter_tot += sweep_quarter[i] * h[i]
    sweep_quarter_tot = sweep_quarter_tot/(b/2)

    return b, AR_wing, sweep_beta, sweep_beta_tot, c[0], taper_ratio, sweep_quarter_tot, c[-1], y, leading_edge, trailing_edge, quarter_line

def getCalageAngle(CL):
    b_wing, AR_wing, _, sweep_beta_wing, c_root_wing, taper_ratio_wing, sweep_quarter_wing, c_tip_wing, _, _, _, _ = wingGeometry()
    _, Cl_fuselage, Cd_fuselage, Cl_max_fus, _ = fuselageCL()

    #cl_alpha wing
    cl_alpha_wing = (1.4073-0.269)/(5+5) * (180/np.pi)
    alpha_L0_wing = -7*(np.pi/180)
 

    k = (beta * cl_alpha_wing)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR_wing)) + np.sqrt((1/((k * np.cos(sweep_beta_wing))))**2 + ((2/(beta * AR_wing))**2) )))/beta

    #alpha_L0root = ((surface_fuselage*alpha_L0_fus) + (surface_wing*alpha_L0_wing))/surface_total
    #CL/a + alpha_L0_wing
    alpha_root = ((CL*surface_total - Cl_fuselage*surface_fuselage)/(surface_wing*a)) + alpha_L0_wing

    return alpha_root

def wingPlot(wing_plot):
    if wing_plot == False:
        return
    _, _, _, _, _, _, _, _, y, leading_edge, trailing_edge, quarter_line = wingGeometry() 
    plt.plot(y, leading_edge)
    plt.plot(y, trailing_edge, color='green')
    plt.plot(y, quarter_line, color='red')

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
    b, AR_wing, sweep_beta, sweep_beta_tot, c_root, taper_ratio, sweep_quarter, c_tip, _, _, _, _ = wingGeometry()
    AoA_wing = getCalageAngle(Cl_max)
    
    # --- Twist angle --- #
    twist_angle = -3 * (np.pi/180) #[°]
    alpha_01 = -0.17
    eta_a_tip = twist_angle 
    alpha_L0 = alpha_l0 + (alpha_01 * twist_angle)
    
    # --- Lift --- #
    AoA = np.linspace(-10, 10, 51) * ((np.pi)/180)
    CL_w = np.zeros(len(AoA))
    
    k = (beta * cl_alpha)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR_wing)) + np.sqrt((1/((k * np.cos(sweep_beta_tot))))**2 + ((2/(beta * AR_wing))**2) )))/beta

    for i in range(len(AoA)):    
        CL_w[i] = a*((AoA[i] + AoA_wing) - alpha_L0)

    CL_w0 = a*((AoA_wing * (np.pi/180)) - alpha_L0) # choose of AoA of the wing 

    CL_max = np.cos(sweep_quarter) * 0.95 * ((cl_max + cl_max)/2)

    if cl_plot:
        plt.plot(AoA*(180/(np.pi)), CL_w)
        #plt.scatter(CL_max, CL_CD_max, marker="x", color="r")
        plt.xlabel('$AoA$')
        plt.ylabel('$Cl_w$')
        plt.show()
    
    return CL_w, CL_w0, CD_wing, CL_max, alpha_L0, a


def totalGeometry():
    _, _, _, _, _, _, _, _, y_wing, leading_wing, trailing_wing, quarter_wing = wingGeometry() 
    _, _, _, _, _, _, _, y_fus, leading_fus, trailing_fus, quarter_fus = fusGeometry() 

    y = np.concatenate((y_fus, y_wing + y_fus[-1]))
    leading_edge = np.concatenate((leading_fus, leading_wing + leading_fus[-1]))
    trailing_edge = np.concatenate((trailing_fus, trailing_wing + leading_fus[-1]))
    quarter_chord = np.concatenate((quarter_fus, quarter_wing + leading_fus[-1]))
    return y, leading_edge, trailing_edge, quarter_chord


def getMAC():
    _, _, _, _, _, _, _, _, y_wing, leading_wing, trailing_wing, quarter_wing = wingGeometry() 
    _, _, _, _, _, _, _, y_fus, leading_fus, trailing_fus, quarter_fus = fusGeometry() 

    # -- fuselage -- #
    c_fus = trailing_fus - leading_fus
    MAC_fus = (2/surface_fuselage) * trapz(c_fus**2, y_fus) #numerical integration via method of trapez
    cy_fus = c_fus*y_fus
    yac_fus = (2/surface_fuselage) * trapz(cy_fus, y_fus)

    xac_fus = MAC_fus*0.1  # keep attention that it is an estimation the table don't give the value for this very low AR
    x_tmp = leading_fus[np.argmin(abs(y_fus - yac_fus))]
    xac_fus = x_tmp+xac_fus 

    # -- wing -- #
    c_wing = trailing_wing - leading_wing
    MAC_wing = (2/surface_wing) * trapz(c_wing**2, y_wing) #numerical integration via method of trapez
    cy_wing = c_wing*y_wing
    yac_wing = (2/surface_wing) * trapz(cy_wing, y_wing)
    xac_wing = MAC_wing*0.23
    
    x_tmp = leading_wing[np.argmin(abs(y_wing - yac_wing))]
    xac_wing = x_tmp+xac_wing + (cabin_lenght - c_wing[0])

    y_wing = y_wing + y_fus[-1]

    # -- total -- #
    c = np.concatenate((c_fus, c_wing[1:]))
    y = np.concatenate((y_fus, y_wing[1:]))
    leading = np.concatenate((leading_fus, leading_wing+leading_fus[-1]))
    MAC = (2/surface_total) * trapz(c**2, y)
    cy = c*y
    yac = (2/surface_total) * trapz(cy, y)
    xac = 0.1*MAC # keep attention that it is an estimation the table don't give the value for this very low AR

    x_tmp = leading[np.argmin(abs(y - yac))]
    xac = x_tmp+xac 

    return MAC_fus, yac_fus, xac_fus, MAC_wing, yac_wing + cabin_width/2, xac_wing, MAC, yac, xac

def plotAllWing(wing_plot):
    if wing_plot == False:
        return
    _, _, _, _, _, _, _, y_fus, leading_edge_fus, trailing_edge_fus, quarter_line_fus = fusGeometry() 
    plt.plot(y_fus, leading_edge_fus, color='blue')
    plt.plot(y_fus, trailing_edge_fus, color='green')
    #plt.plot(y_fus, quarter_line_fus, color='red')

    _, _, _, _, _, _, _, _, y_wing, leading_edge_wing, trailing_edge_wing, quarter_line_wing = wingGeometry() 
    plt.plot(y_wing + y_fus[-1], leading_edge_wing + leading_edge_fus[-1], color='blue')
    plt.plot(y_wing + y_fus[-1], trailing_edge_wing + leading_edge_fus[-1], color='green')
    #plt.plot(y_wing + y_fus[-1], quarter_line_wing + leading_edge_fus[-1], color='red')

    plt.plot([y_fus[-1], y_wing[0] + y_fus[-1]], [trailing_edge_fus[-1], trailing_edge_wing[0] + leading_edge_fus[-1]], color='green')

    # --- left part --- #
    plt.plot(-y_fus, leading_edge_fus, color='blue')
    plt.plot(-y_fus, trailing_edge_fus, color='green')
    #plt.plot(-y_fus, quarter_line_fus, color='red')

    plt.plot(-(y_wing + y_fus[-1]), (leading_edge_wing + leading_edge_fus[-1]), color='blue')
    plt.plot(-(y_wing + y_fus[-1]), (trailing_edge_wing + leading_edge_fus[-1]), color='green')
    #plt.plot(-(y_wing + y_fus[-1]), (quarter_line_wing + leading_edge_fus[-1]), color='red')

    plt.plot([-y_fus[-1], -(y_wing[0] + y_fus[-1])], [trailing_edge_fus[-1], (trailing_edge_wing[0] + leading_edge_fus[-1])], color='green')

    MAC_fus, yac_fus, xac_fus, MAC_wing, yac_wing, xac_wing, MAC, yac, xac = getMAC()

    plt.scatter(yac_fus,  xac_fus+ (yac_fus*np.tan(60*(np.pi/180))), color='red')
    plt.scatter(yac_wing,  xac_wing, color='orange')
    leading_edge_fus_x = np.interp(yac_fus, y_fus, leading_edge_fus)
    leading_edge_wing_x = np.interp(yac_wing+3.5, y_wing, leading_edge_wing)
    plt.plot((yac_fus, yac_fus), (leading_edge_fus_x, leading_edge_fus_x + MAC_fus), color='red')
    plt.plot((yac_wing, yac_wing), (leading_edge_wing_x , leading_edge_wing_x + MAC_wing), color='orange')

    plt.xlabel('$Y$')
    plt.ylabel('$X$')
    # Fixer l'échelle des axes
    plt.axis('equal')
    plt.show() 
    return
plotAllWing(wing_plot)


def get_Lift_and_drag(AR, delta):
    AoA = np.linspace(-10, 10, 51) * ((np.pi)/180)

    Cl_wing, Cl_wing_0, Cd_wing, Cl_max_wing, _, _ = wingCL()
    Cl_fuselage, CL_fus_0, Cd_fuselage, Cl_max_fus, _ = fuselageCL()

    # --- total lift computation --- #
    Cl_tot = np.zeros(len(AoA))
    Cl_tot = ((Cl_wing*surface_wing) + (Cl_fuselage*surface_fuselage))/surface_total 
    
    Cl_tot0 = float(Cl_tot[np.where(abs(AoA) <= 1e-12)])
    AoA_L0 = float(AoA[np.where(abs(Cl_tot) <= 6e-3)])*(180/np.pi)

    Cl_max = ((Cl_max_wing*surface_wing) + (Cl_max_fus*surface_fuselage))/surface_total 
    
    # --- total drag computation --- #
    Cd_induce = ((Cl_tot**2)/(np.pi* AR)) * (1+delta)
    #print(f"delta Cd = {(((Cl_tot[np.where(abs(AoA) <= 1e-12)]**2)/(np.pi* (AR-winglet()))) * (1+delta))-Cd_induce[np.where(abs(AoA) <= 1e-12)]}")
    Cd_tot = np.zeros(len(AoA))
    Cd_tot = Cd_induce + (((Cd_wing*surface_wing) + (Cd_fuselage*surface_fuselage))/surface_total)
    Cd_tot0 = float(Cd_tot[np.where(abs(AoA) <= 1e-12)])

    return Cl_tot0, Cd_tot0, Cl_max, AoA_L0, Cl_tot, Cd_tot, AoA

def plotLiftDrag(lift_and_drag_plots):
    if lift_and_drag_plots == False:
        return
    Cl_tot0, Cd_tot0, Cl_max, AoA_L0, Cl_tot, Cd_tot, AoA =  get_Lift_and_drag(AR, delta)
    plt.figure(figsize=(8,5))
    plt.plot(Cl_tot, Cl_tot/Cd_tot)
    plt.xlabel('$CL$')
    plt.ylabel('$CL/CD$')
    plt.show()

    plt.figure(figsize=(8,5))
    plt.plot(Cd_tot, Cl_tot)
    plt.xlabel('$CD$')
    plt.ylabel('$CL$')
    plt.show()

    plt.figure(figsize=(8,5))
    plt.plot(AoA, Cl_tot)
    plt.xlabel('$AoA$')
    plt.ylabel('$CL$')
    plt.show()

    plt.figure(figsize=(8,5))
    plt.plot(AoA, Cd_tot)
    plt.xlabel('$AoA$')
    plt.ylabel('$CL$')
    plt.show()
    return
plotLiftDrag(lift_and_drag_plots)

def wingMaxthickness():
    _, CL, _, _, _, _ = wingCL()
    b, AR_wing, sweep_beta, sweep_beta_tot, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = wingGeometry()
    M_star = 1.15-(CL/(4*(np.cos(sweep_quarter)**2)))
    t_bar_over_c =  (3/(10*M)) * np.cbrt((1/(M*np.cos(sweep_quarter))) - (M*np.cos(sweep_quarter))) * (1-(((5 + (M**2)*(np.cos(sweep_quarter)**2))/(5 + (M_star**2)))**3.5))**(2/3)
    t_root = t_bar_over_c*c_root
    t_tip = t_bar_over_c*c_tip 
    
    return t_root, t_tip, t_bar_over_c


def wingFuelvolume():
    b, AR_wing, sweep_beta, sweep_beta_tot, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = wingGeometry()
    t_root, t_tip, t_bar_over_c = wingMaxthickness()
    tau = 1
    V_fuel = 0.54*((surface_wing**2)/b) * t_bar_over_c * ((1+(taper_ratio*np.sqrt(tau))+((taper_ratio**2)*tau))/(1+taper_ratio)**2)
    return V_fuel


def wingSurfaceWetted():
    b, AR_wing, sweep_beta, sweep_beta_tot, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = wingGeometry()
    t_root, t_tip, t_bar_over_c = wingMaxthickness()
    S_wet_wing = 2*surface_wing*(1+ (1/4) * ((t_bar_over_c + (t_bar_over_c*taper_ratio))/(1+taper_ratio)))

    b, AR_fuselage, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = fusGeometry()
    S_wet_fus = 2*surface_fuselage*(1+ (1/4) * ((0.18 + (0.18*taper_ratio))/(1+taper_ratio)))

    S_wet_tot = S_wet_wing + S_wet_fus
    return S_wet_wing, S_wet_fus, S_wet_tot

def stallVelocity():
    b, AR_wing, sweep_beta, sweep_beta_tot, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = wingGeometry()
    _, _, Cl_max, _ , _, _, _= get_Lift_and_drag(AR, delta)
    rho_sl, T = air_density(0)

    Vs = np.sqrt((weight/surface_total) * (2/rho_sl) * (1/(1.133*Cl_max)))
    
    Cl_max0 = 2 * np.cos(sweep_quarter)
    W0 = weight_empty #landing weight
    Vs0 = np.sqrt((W0/surface_total) * (2/rho_sl) * (1/(1.133*Cl_max0)))
    return Vs, Vs0

def getReynold(altitude, c):
    rho, T = air_density(altitude)
    U_inf = true_airspeed_at_altitude(altitude)
    
    p_atmo = 99333      #[Pa]
    T = 12.0 + 273.15   #[k]
    rho = p_atmo/(287*T)
    mu = 1.716e-5 * (T/273.15)**(3/2) * ((273.15 + 110.4)/(T + 110.4)) # Sutherland's law
    Re = (rho * U_inf * c) / mu
    return Re

def printFunction():

    print(f"New AR = {AR:.3f} [-]")
    print(f"Beta = {beta:.3f} [-]\n")
    print(f"Cl max used = {Cl_max:.2f} [-]\n")
    print(f"Total area = {surface_total:.2f} [m^2]")
    print(f"Surface of fuselage = {surface_fuselage:.2f} [m^2]")
    print(f"Surface of wing = {surface_wing:.2f} [m^2] \n")
    print(f"Compressibility parameter: {beta:.3f} [-]")
    
    b, AR_wing, sweep_beta, sweep_beta_tot, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = wingGeometry()
    Cl_wing, Cl_wing_0, Cd_wing, Cl_max_wing, alpha_L0, a_wing = wingCL()
    print("\n-------------- wing values --------------\n")
    print(f"\nAR wing: {AR_wing:.3f} [-] \nCL_w0 wing = {Cl_wing_0:.3f} [-]\n")
    print(f"Cord at wing root: {c_root:.3f} [m]\nCorde at wing tip: {c_tip:.3f} [m]")
    print(f"Taper ratio: {taper_ratio:.3f} [-]")
    print(f"sweep quater: {sweep_quarter*(180/np.pi):.3f} [°]")
    print(f"Wing lift coefficient derivative: {a_wing:.3f}")
    print(f"Alpha_L0: {alpha_L0*(180/np.pi):.3f}")
    

    b, AR_fuselage, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = fusGeometry()
    Cl_fuselage, Cl_fuselage_0, Cd_fuselage, Cl_max_fus, a_fus = fuselageCL()
    print("\n-------------- fuselage values --------------\n")
    print(f"\nAR fuselage: {AR_fuselage:.3f} [-]\nCL_w0 fuselage = {Cl_fuselage_0:.3f} [-]\n")
    print(f"Cord at fuselage root: {c_root:.3f} [m]\nCorde at fuselage tip: {c_tip:.3f} [m]")
    print(f"Taper ratio: {taper_ratio:.3f} [-]")
    print(f"sweep quater: {sweep_quarter*(180/np.pi):.3f}")
    print(f"Fuselage lift coefficient derivative: {a_fus:.3f}\n")

    
    MAC_fus, yac_fus, xac_fus, MAC_wing, yac_wing, xac_wing, MAC, yac, xac = getMAC()
    print("\n-------------- MAC values --------------\n")
    print(f"MAC fus: {MAC_fus:.3f} [m]\nYac fus: {yac_fus:.3f} [m]\nXac fus: {xac_fus:.3f} [m]\n")
    print(f"MAC wing: {MAC_wing:.3f} [m]\nYac wing: {yac_wing:.3f} [m]\nXac wing: {xac_wing:.3f} [m]\n")
    print(f"MAC: {MAC:.3f} [m]\nYac: {yac:.3f} [m]\nXac: {xac:.3f} [m]\n")
    

    delta = 0.005 #graph slide 61 lecture 6 aerodinimics
    lift_coef, drag_coef, CL_max, AoA_L0, cl, _, aoa = get_Lift_and_drag(AR, delta)
    print(f"\n CL = {lift_coef:.3f}[-] \n CD = {drag_coef:.3f}[-] \n")
    print(f"Cl max: {CL_max:.3f} [-]")
    print(f"Lift coefficient derivative CL_alfa: {(cl[-1] - cl[0])/(aoa[-1] - aoa[0]):.3f} [deg^-1] \n")
    
    
    t_root, t_tip,t_bar_over_C = wingMaxthickness()
    print(f"Thickness root: {t_root:.3f}, thickness tip: {t_tip:.3f}")
    print(f"tbar over c: {t_bar_over_C:.3f} [-]\n")
    #print(f"Mean thinckness: {(t_root+t_tip)/2:.3f}\n")
    
    V_fuel = wingFuelvolume()
    print(f"Fuel volume in wing: {V_fuel:.3f} [m^3]\n")

    _, _, Swetted = wingSurfaceWetted()
    print(f"Surface wetted: {Swetted:.3f} [m^2]\n")
    
    Vs, Vs0 = stallVelocity()
    print(f"Stall velocity: {Vs:.3f} [m/s]\nStall velocity in approach config: {Vs0:.3f} [m/s]\n")
    
    AoA_root = getCalageAngle(Cl_max)
    print(f"AoA root needed: {AoA_root*(180/np.pi):.3} [°]")
    print(f"AoA zero lift: {AoA_L0:.3f} [°]")
    
    Re = getReynold(alti, MAC)
    print(f"Re_mac: {Re:.3f} [-]")

    return

printFunction()