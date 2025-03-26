import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import trapz
import pandas as pd

#---Guesses---#
span_max = 28.95           #[m] Span  max span for airport
cabin_width = 9         #[m] 
cabin_lenght = 16.8     #[m] 
#AR = 4           #Aspect ratio (guess)
#weight =   655202.751944230  #526898.7380202 #[n]          471511.49122 #  #[N] = 106000lb (guess from weight code)
weight_empty = 526898.7380202 #60452.314059821154 * 9.81 #[N]S 
alti = 12500            #[m]
M = 0.9                #[-] Mach number
R = 287                 #[m^2/s^2K]
gamma = 1.4
e = 0.85                #Ostxald's efficiency factor
delta = 0.005           #graph slide 61 lecture 6 aerodinimics
#sweep_LE_fus = 40   #[°] sweep angle fuselage
#sweep_LE_wing = 30 #[°] sweep angle wing
twist_angle = -3         #[°] twist angle
#Lambda = 0.6           # [-] taper ratio
beta = np.sqrt(1-(M**2))


#--- Figure settings ---#
plt.rcParams.update({
    "text.usetex": True,              # Use LaTeX for all text rendering
    "font.family": "serif",           # Use LaTeX's default font family
    "font.serif": ["Times"],          # Use Computer Modern for a LaTeX-like font
    "font.size": 22,                  # Global font size to match LaTeX
    "axes.titlesize": 22,             # Font size for title
    "axes.labelsize": 20,             # Font size for axis labels
    "xtick.labelsize": 20,            # Font size for x-axis ticks
    "ytick.labelsize": 20,            # Font size for y-axis ticks
    "legend.fontsize": 18             # Font size for legend
})



#---Commande---#
polar_Cl_Cd = False
wing_plot = False
lift_and_drag_plots = False
plot_airfoil = False

#---Code---#

def winglet(AR):
    #formule from Snorri
    h = 0.5 #[m] Choice
    delta_AR = 1.9*(h/span_max) * AR
    return delta_AR


#def getSweep():
#    return sweep_LE_fus, sweep_LE_wing

def getAirfoilFus():
    airfoil = 2
    if airfoil == 1:
        cl_alpha = ((1.0498+0.2062)/(5+5)) * (180/np.pi) # SC(2) 0518 M0 Re12M C_m = -0.1158
        cl_max = 1.87
        alpha_L0 = -3.5 * (np.pi/180) #[rad] 
        CD_fuselage = 0.00636 
        cm = -0.1158
    elif airfoil == 2:
        cl_alpha = (0.7712+0.2134)/(5+5) * (180/np.pi)# NACA45118 M0.85 Re12M cm = -0.0026
        cl_max = 1.7125
        alpha_L0 = -3 * (np.pi/180) #[rad]
        CD_fuselage = 0.0026
        cm = -0.0026
    elif airfoil == 3:
        cl_alpha = (0.7236+0.2617)/(5+5) * (180/np.pi)# NACA35118 M0.85 Re12M cm = -0.0022
        cl_max = 1.6876
        alpha_L0 = -2.5 * (np.pi/180) #[rad]
        CD_fuselage = 0.00624
        cm = -0.0022
    elif airfoil == 4:
        cl_alpha = (0.6332+0.3581)/(5+5) * (180/np.pi) # NACA25118 M0.85 Re12M cm = -0.0011
        cl_max = 1.6841
        alpha_L0 = -2 * (np.pi/180) #[rad]
        CD_fuselage = 0.00597
        cm = -0.0011
    elif airfoil == 5:
        cl_alpha = (0.4819+0.1866)/(5+5) * (180/np.pi) # NACA67-318 M0.85 Re12M cm = -0.0355
        cl_max = 1.4739
        alpha_L0 = -3 * (np.pi/180) #[rad]
        CD_fuselage = 0.00431
        cm = -0.0355
    elif airfoil == 6:
        cl_alpha = (0.3769+0.2766)/(5+5) * (180/np.pi) # NACA 67-118 M0.85 Re12M cm = -0.0102
        cl_max = 1.447
        alpha_L0 = -1 * (np.pi/180) #[rad]
        CD_fuselage = 0.00286
        cm = -0.0102
    elif airfoil == 7:
        cl_alpha = (0.82-0.2)/(1.5+2) * (180/np.pi) # SC(2)-0714 M0.75 Re6M
        alpha_L0 = -3.5 * (np.pi/180) #[rad]
        CD_fuselage = 0.01 
    elif airfoil == 8:
        cl_alpha = ((0.8+0.2)/(5+5)) *(180/np.pi) # Eppler 642 M0 Re1M C_m = -0.058
        alpha_L0 = -4 * (np.pi/180) #[rad] 
        CD_fuselage = 0.01 
        cm = -0.058
    return  cl_alpha, cl_max, alpha_L0, CD_fuselage, cm

def getAirfoilWing():
    airfoil = 1
    if airfoil == 1:
        cl_alpha = (1.1117+0.0543)/(5+5) * (180/np.pi) # SC(2)-0710 M0.85 Re12M cm = -0.129
        cl_max = 2.243
        alpha_l0 = -4.5*(np.pi/180)
        CD_wing = 0.00907 #at 6° aoa  and at 0° aoa = 0.006
        cm = -0.129
    elif airfoil == 2:
        cl_alpha = (1.4073-0.269)/(5+5) * (180/np.pi) # SC(2)-1010 M0.85 Re12M
        cl_max = 1.4
        alpha_l0 = -7*(np.pi/180)
        CD_wing = 0.01091 #at 6° aoa  and at 0° aoa = 0.006
        cm = -0.1518
    elif airfoil == 3:
        cl_alpha = (0.5765-0)/(5+0) * (180/np.pi) # SC(2)-0012 M0.85 Re12M
        alpha_l0 = 0
        CD_wing = 0.0059 

    elif airfoil == 4:
        cl_alpha = (0.8293+0.3558)/(5+5) * (180/np.pi) # SC(2)-0410 M0.85 Re12M
        alpha_l0 = -2*(np.pi/180)
        cl_max = 2.1599
        CD_wing = 0.006
        cm = -0.077
    elif airfoil == 5:
        cl_alpha = (1.7-0)/(4+12) * (180/np.pi) # NACA 64209 M0.85 Re6M
        alpha_l0 = -12*(np.pi/180)
        CD_wing = 0.007 
    return cl_alpha, cl_max, alpha_l0, CD_wing, cm

def getAirfoilMiddle():
    airfoil = 1
    if airfoil == 1: # NACA65-713 
        cl_alpha = (0.9712+0.1108)/(5+5) * (180/np.pi) 
        cl_max = 1.7
        alpha_l0 = -3.5*(np.pi/180)
        CD_wing = 0.0054 #at 0° aoa  
        cm = -0.064
    return cl_alpha, cl_max, alpha_l0, CD_wing, cm

import pandas as pd
import matplotlib.pyplot as plt

def plotAirfoil(plot_airfoil):
    if not plot_airfoil:
        return

    # Lire les fichiers CSV
    data_fus = pd.read_csv("NACA45118_XYZ.csv")  
    data_wing = pd.read_csv("sc20710_XYZ.csv")

    # ---- Première figure : Airfoil fuselage ----
    fig, ax = plt.subplots(figsize=(11, 3), dpi=300)  # Créer une figure et un axe
    ax.set_aspect('equal')  # Assurer une échelle identique sur les axes
    ax.plot(data_fus.iloc[:, 0], data_fus.iloc[:, 1])
    ax.set_xlabel('$x/c [-]$')  # Légende de l'axe des abscisses
    ax.set_ylabel('$y [-]$')  # Légende de l'axe des ordonnées
    plt.savefig("/Users/antoinevanhoye/Documents/M1/PI/integrated_project/Airfoils/airfoil_fus.pdf", dpi=300)  # Sauvegarder avant d'afficher
    #plt.show()

    # ---- Deuxième figure : Airfoil aile ----
    fig, ax = plt.subplots(figsize=(11, 2.7), dpi=300)
    ax.set_aspect('equal')  # Assurer une échelle identique
    ax.plot(data_wing.iloc[:, 0], data_wing.iloc[:, 1])
    ax.set_yticks([-0.08, 0, 0.08])
    ax.set_xlabel('$x/c [-]$')  # Légende de l'axe des abscisses
    ax.set_ylabel('$y [-]$')  # Légende de l'axe des ordonnées
    plt.savefig("/Users/antoinevanhoye/Documents/M1/PI/integrated_project/Airfoils/airfoil_wing.pdf", dpi=300)  
    #plt.show()
plotAirfoil(plot_airfoil)

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
def true_airspeed_at_altitude(altitude, M):
    T = air_density(altitude)[1]
    a = np.sqrt(gamma * R * T)  # [m/s] Speed of sound
    v = M * a  # [m/s] Aircraft velocity
    return v
v = true_airspeed_at_altitude(alti, M)


def guess_CL_max(AR):
    CL = np.linspace(0, 1.5, 100)
    CD0 = 0.007              # Zero lift drag coeff
    CD = CD0 + CL*CL / (np.pi * AR * 0.85) # (Zero Lift drag + Induce drag)

    CL_CD = CL / CD
    # Trouver l'indice du maximum de CL/CD
    max_index = np.argmax(CL_CD)

    # Extraire le CL correspondant au maximum de CL/CD
    CL_max = CL[max_index]
    CL_CD_max = CL_CD[max_index]

    # Afficher les résultats
    #print(f"Le maximum de CL/CD est {CL_CD_max:.2f} pour CL = {CL_max:.2f} \n")

    if polar_Cl_Cd:
        plt.plot(CL, CL/CD)
        plt.scatter(CL_max, CL_CD_max, marker="x", color="r")
        plt.xlabel('$CL$')
        plt.ylabel('$CL/CD$')
        plt.show()

    Cl_max = CL_max -  (CL_max*0.1)

    return Cl_max


def detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    surface_wing, AR, taper_wing, croot, ctip, chord_middle, sweep_leading, sweep_beta = wingGeometryIDEAL(Cl, weight, sweep_quarter_wing)

    surface_wing_ideal = span_max**2/AR #2*weight/(rho*(v**2)*Cl_max) 
    surface_fuselage = (cabin_lenght +(cabin_lenght - ((cabin_width/2)/np.tan((90-sweep_LE_fus)*np.pi/180))))/2 * (cabin_width) #(cabin_width * cabin_lenght)
    surface_wing = (chord_middle + ctip) * (span_max - cabin_width)*0.5
    surface_total = surface_wing + surface_fuselage
    #surface_fuselage = surface_total - surface_wing_ideal
    return surface_wing_ideal, surface_fuselage, surface_wing, surface_total

"""
def getCl(AR, sweep_LE_fus, weight):
    surface_total, surface_fuselage, surface_wing = detSurfac(AR, sweep_LE_fus)
    Cl = (2*weight)/(rho*(v**2)*surface_total) 
    return Cl
"""
def getSurfaceFuselage(sweep_LE_fus):
    surface_fuselage = (cabin_lenght +(cabin_lenght - ((cabin_width/2)/np.tan((90-sweep_LE_fus)*np.pi/180))))/2 * (cabin_width) #(cabin_width * cabin_lenght)
    
    return surface_fuselage 


def getSurface_And_AR(Cl, weight):
    rho, T = air_density(alti)
    v = true_airspeed_at_altitude(alti, M)
    surface_wing_ideal = weight/(0.5*rho*(v**2)*Cl)
    AR = span_max**2/surface_wing_ideal
    return surface_wing_ideal, AR


def fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    #surface_fuselage = getSurfaceFuselage(sweep_LE_fus)
    surface_wing_ideal, surface_fuselage, surface_wing, surface_total = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    b = cabin_width 
    AR_fuselage = (b**2)/surface_fuselage
    sweep_leading =  sweep_LE_fus #[°]
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

def fusPlot(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    if wing_plot == False:
        return
    _, _, _, _, _, _, _, y, leading_edge, trailing_edge, quarter_line = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight) 
    plt.plot(y, leading_edge)
    plt.plot(y, trailing_edge, color='green')
    plt.plot(y, quarter_line, color='red',)

    plt.xlabel('$Y$')
    plt.ylabel('$X$')
    # Fixer l'échelle des axes
    plt.axis('equal')
    plt.show() 
    return


def fuselageCL(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    # --- airfoil --- #
    cl_alpha, cl_max, alpha_L0, CD_fuselage, cm = getAirfoilFus()
    b, AR_fuselage, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, _, _, _, _ = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight)

    # --- Lift --- #
    AoA = np.linspace(-10, 10, 21) * ((np.pi)/180)
    CL_w = np.zeros(len(AoA))
    
    k = (beta * cl_alpha)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR_fuselage)) + np.sqrt( ((1/(k * np.cos(sweep_beta)))**2) + ((2/(beta * AR_fuselage))**2) )))/beta

    for i in range(len(AoA)):    
        CL_w[i] = a*((AoA[i]+ (0*(np.pi/180))) - alpha_L0)
        CL_w0 = a*(0 - alpha_L0)
    
    CL_w0 = a*(0*(np.pi/180) - alpha_L0)
    

    CL_max = np.cos(sweep_quarter) * 0.95 * ((cl_max + cl_max)/2)
    
    return CL_w, CL_w0, CD_fuselage, CL_max, a



def wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, weight):
    b_fus, AR_fuselage, sweep_beta_fus, c_root_fus, taper_ratio_fus, sweep_quarter_fus, c_tip_fus, y_fus, leading_edge_fus, trailing_edge_fus, quarter_line_fus = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    surface_wing_ideal, AR, taper_wing, croot, ctip, chord_middle, sweep_leading, sweep_beta = wingGeometryIDEAL(Cl, weight, sweep_quarter_wing)
    surface_wing_ideal, surface_fuselage, surface_wing, surface_total = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    b = span_max - cabin_width 
    AR_wing = AR#(b**2)/surface_wing
    
    # ---- Complex wing ---- #

    #h = [0.5, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, (b*0.5) - 3] #[m]
    #Yposition = [0, h[0], h[0]+h[1], b/2]
    #c = [c_tip_fus, 8.5, 7.8, 7.3, 6.3, 5.7, 5.3, 5.0 , 0] #[c_tip_fus, 7, 6.3, 5.8, 4.8, 4.2, 3.8, 3.5 , 0]


    #h = [0.2, 0.4, 0.6, 0.8, (b*0.5) - 1.4] #[0.5, 1.0, 1.5, 2.0, 2.5 , (b*0.5) - 3.5] #[m]
    #h = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, (b*0.5) - 1.6]
    h = [(b*0.5)]
    #Yposition = [0, h[0], h[0]+h[1], b/2]
    #c = [c_tip_fus, 7.51864968, 6.78573899, 6.18381709, 5.65259713, 0] #[c_tip_fus, 7, 6.3, 5.8, 4.8, 4.2, 3.8, 3.5 , 0]
    #c = [c_tip_fus, 7.500065654471188, 6.6515342770866, 5.919295847694265, 5.250786658155237, 4.625425101990652, 4.03259169538853, 3.4659748031966577, 2.9214890949003856, 0]
    croot = chord_middle#surface_wing/(b*0.5*(1+taper_wing))
    ctip = ctip#taper_wing*croot
    
    #c = [c_tip_fus, 4.5, 0]
    c = [croot, ctip]
    S = np.zeros(len(h))
    S_sum =0
    for i in range(len(c)-2):
        S[i] = (c[i]+ c[i+1])*0.5*h[i]
        S_sum += S[i]
    #S1 = (c[0]+ c[1])*0.5*h[0]
    #S2 = (c[1] + c[2])*0.5*h[1]
    S[-1] = (surface_wing*0.5) - S_sum
    #c[-1] = (2/h[-1]) * S[-1] - c[-2] #c_tip

    #sweep_leading = sweep_leading*((np.pi)/180)
    sweep_quarter = np.zeros(len(c)-1)
    sweep_trailing = np.zeros(len(c)-1)
    sweep_beta = np.zeros(len(c)-1)

    for i in range(len(c)-1):
        sweep_quarter[i] = np.arctan2(h[i], ((c[i]*0.25) - ((h[i]*np.tan(sweep_leading)) + (c[i+1]*0.25)))) - (np.pi/2)#np.arctan(np.tan(sweep_leading) + (4/AR_tmp) * (((1-taper_tmp)/(1+taper_tmp)) * (0 - 0.25)))
        sweep_trailing[i] = np.arctan2(h[i], (c[i] - ((h[i]*np.tan(sweep_leading)) + c[i+1]))) - (np.pi/2) #np.arctan(np.tan(sweep_quarter[i]) + (4/AR_tmp) * (((1-taper_tmp)/(1+taper_tmp)) * (0.25 - 1)))
        sweep_beta[i] = np.arctan(np.tan(sweep_quarter[i])/beta) #np.arctan2(np.tan(sweep_quarter[i]), beta)

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
    
    return b, AR_wing, sweep_beta, sweep_beta_tot, c[0], taper_ratio, sweep_quarter_tot, c[-1], y, leading_edge, trailing_edge, quarter_line, c, h


def wingGeometryIDEAL(lift_coef, weight, sweep_quarter_wing):
    surface_wing, AR = getSurface_And_AR(lift_coef, weight)
    b = span_max #- cabin_width 
    
    # ---- Complex wing ---- #
    taper_wing = 0.4

    h = [b*0.5]
    
    croot = surface_wing/(b*0.5*(1+taper_wing))
    ctip = taper_wing*croot

    c = [croot, ctip]
    
    sweep_quarter = sweep_quarter_wing*(np.pi/180)
    
    #sweep_quarter[i] = np.arctan2(h[i], ((c[i]*0.25) - ((h[i]*np.tan(sweep_leading)) + (c[i+1]*0.25)))) - (np.pi/2)

    sweep_leading = np.arctan(np.tan(sweep_quarter) + (4/AR) * (((1-taper_wing)/(1+taper_wing)) * (0.25 - 0)))
    sweep_trailing = np.arctan2(h[0], (c[0] - ((h[0]*np.tan(sweep_leading)) + c[1]))) - (np.pi/2) 
    sweep_beta = np.arctan(np.tan(sweep_quarter)/beta) 

    chord_middle = np.interp(cabin_width/2, [0.0, span_max/2], c)
    return surface_wing, AR, taper_wing, croot, ctip, chord_middle, sweep_leading, sweep_beta


def getCalageAngle(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    surface_wing_ideal, AR, taper_wing, croot, ctip, chord_middle, sweep_leading, sweep_beta = wingGeometryIDEAL(Cl, weight, sweep_quarter_wing)
    #b_wing, AR_wing, _, sweep_beta_wing, c_root_wing, taper_ratio_wing, sweep_quarter_wing, c_tip_wing, _, _, _, _, _, _= wingGeometry(AR,sweep_LE_fus, sweep_LE_wing)
    _, Cl_fuselage, Cd_fuselage, Cl_max_fus, _ = fuselageCL(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    surface_wing_ideal, surface_fuselage, surface_wing, surface_total = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, weight)

    # --- cl_alpha wing --- #
    cl_alpha_wing, cl_max, alpha_l0, CD_wing, cm = getAirfoilWing()

    # --- Twist angle --- #
    alpha_01 = -0.17
    eta_a_tip = twist_angle 
    alpha_L0_wing = alpha_l0 + (alpha_01 * twist_angle * (np.pi/180))
    
    # --- Cl of wing --- #
    
    Cl_wing = Cl - Cl_fuselage#((Cl*surface_wing_ideal) - (Cl_fuselage*surface_fuselage))/surface_wing
    
    # --- Rest --- #
    k = (beta * cl_alpha_wing)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR)) + np.sqrt((1/((k * np.cos(sweep_beta))))**2 + ((2/(beta * AR))**2) )))/beta
    
    alpha_root = (Cl_wing/a) + alpha_L0_wing
    
    #((CL*surface_total - Cl_fuselage*surface_fuselage)/(surface_wing*a)) + alpha_L0_wing     # juste all in one equation 
    return alpha_root, a

def wingPlot(wing_plot, Cl,sweep_LE_fus, sweep_quarter_wing, weight):
    if wing_plot == False:
        return
    _, _, _, _, _, _, _, _, y, leading_edge, trailing_edge, quarter_line, _, _ = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, weight) 
    plt.plot(y, leading_edge)
    plt.plot(y, trailing_edge, color='green')
    plt.plot(y, quarter_line, color='red')

    plt.xlabel('$Y$')
    plt.ylabel('$X$')
    # Fixer l'échelle des axes
    plt.axis('equal')
    plt.show() 

    return


def wingCL(Cl,sweep_LE_fus, sweep_quarter_wing, weight):
    # ----- Airfoil ----- # 
    cl_alpha, cl_max, alpha_l0, CD_wing, cm = getAirfoilWing()
    surface_wing_ideal, AR = getSurface_And_AR(Cl, weight)
    surface_wing_ideal, AR, taper_wing, croot, ctip, chord_middle, sweep_leading, sweep_beta = wingGeometryIDEAL(Cl, weight, sweep_quarter_wing)

    AoA_wing, _ = getCalageAngle(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    
    # --- Twist angle --- #
    alpha_01 = -0.17
    eta_a_tip = twist_angle 
    alpha_L0 = alpha_l0 + (alpha_01 * twist_angle * (np.pi/180))
    
    # --- Lift --- #
    AoA = np.linspace(-10, 10, 21) * ((np.pi)/180)
    CL_w = np.zeros(len(AoA))
    
    k = (beta * cl_alpha)/(2*np.pi)
    
    a = ((2*np.pi)/((2/(beta*AR)) + np.sqrt((1/((k * np.cos(sweep_beta))))**2 + ((2/(beta * AR))**2) )))/beta
    
    
    for i in range(len(AoA)):    
        CL_w[i] = a*((AoA[i] + AoA_wing) - alpha_L0)
    
    CL_w0 = a*((AoA_wing) - alpha_L0) # choose of AoA of the wing 
    
    CL_max = np.cos(sweep_quarter_wing*np.pi/180) * 0.95 * ((cl_max + cl_max)/2)
    
    return CL_w, CL_w0, CD_wing, CL_max, alpha_L0, a


def totalGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    surface_wing_ideal, AR, taper_wing, croot, ctip, chord_middle, sweep_leading, sweep_beta = wingGeometryIDEAL(Cl, weight, sweep_quarter_wing)
    _, _, _, _, _, _, _, _, y_wing, leading_wing, trailing_wing, quarter_wing, _, _ = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, weight) 
    _, _, _, _, _, _, _, y_fus, leading_fus, trailing_fus, quarter_fus = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight) 

    y = np.concatenate((y_fus, y_wing + y_fus[-1]))
    leading_edge = np.concatenate((leading_fus, leading_wing + leading_fus[-1]))
    trailing_edge = np.concatenate((trailing_fus, trailing_wing + leading_fus[-1]))
    quarter_chord = np.concatenate((quarter_fus, quarter_wing + leading_fus[-1]))
    return y, leading_edge, trailing_edge, quarter_chord


def getMAC(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    surface_wing_ideal, AR, taper_wing, croot, ctip, chord_middle, sweep_leading, sweep_beta = wingGeometryIDEAL(Cl, weight, sweep_quarter_wing)
    _, _, _, _, _, _, _, _, y_wing, leading_wing, trailing_wing, quarter_wing, _, _ = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, weight) 
    _, _, _, _, _, _, _, y_fus, leading_fus, trailing_fus, quarter_fus = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight) 
    surface_wing_ideal, surface_fuselage, surface_wing, surface_total = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    surface_fuselage = surface_total-surface_wing

    x = leading_fus[-1] - (np.tan(sweep_leading)*cabin_width/2)
    leading_wing_ideal = np.array([x, x + (np.tan(sweep_leading)*(span_max/2))]) # voir page 27 note pi ipad

    # -- fuselage -- #
    c_fus = trailing_fus - leading_fus
    #c_fus = np.array([ c_fus[0] - croot, c_fus[-1] - ctip])
    #y_fus = np.array([y_fus[0], y_fus[-1]])
    MAC_fus = (2/surface_fuselage) * trapz(c_fus**2, y_fus) #numerical integration via method of trapez
    cy_fus = c_fus*y_fus
    yac_fus = (2/surface_fuselage) * trapz(cy_fus, y_fus)
    xac_fus = MAC_fus*0.4  # keep attention that it is an estimation the table don't give the value for this very low AR
    
    #tmp_leading_fus = [leading_fus[0], leading_fus[-1]]
    x_tmp =  np.interp(yac_fus , y_fus, leading_fus)
    xac_fus = x_tmp+xac_fus 

    # -- wing -- #
    y_wing = y_wing + y_fus[-1]

    c_wing = trailing_wing - leading_wing
    MAC_wing = (2/surface_wing) * trapz(c_wing**2, y_wing) #numerical integration via method of trapez
    cy_wing = c_wing*y_wing
    yac_wing = (2/surface_wing) * trapz(cy_wing, y_wing)
    xac_wing = MAC_wing*0.29
    
    x_tmp = np.interp(yac_wing , y_wing, leading_wing)
    xac_wing = x_tmp+xac_wing + (cabin_lenght - c_fus[-1]) 

    
    """# -- wing -- #
    c_wing = np.array([chord_middle, ctip])#trailing_wing - leading_wing
    y_wing = np.array([0, (span_max - cabin_width)/2])
    MAC_wing = (2/surface_wing) * trapz(c_wing**2, y_wing) #numerical integration via method of trapez
    cy_wing = c_wing*y_wing
    yac_wing = (2/surface_wing) * trapz(cy_wing, y_wing)
    xac_wing = MAC_wing*0.28
    
    tmp_leading_wing = [leading_wing[0], leading_wing[-1]]
    x_tmp =  np.interp(yac_wing , y_wing, tmp_leading_wing)
    xac_wing = x_tmp+xac_wing + (cabin_lenght - c_fus[-1]) 

    y_wing = y_wing + y_fus[-1]
    """
    # -- total -- #
    c = np.concatenate((c_fus, c_wing[1:]))
    y = np.concatenate((y_fus, y_wing[1:]))
    leading = np.concatenate((leading_fus, leading_wing+leading_fus[-1]))
    MAC = (2/surface_total) * trapz(c**2, y)
    cy = c*y
    yac = (2/surface_total) * trapz(cy, y)

    #MAC = ((MAC_wing*surface_wing) + (MAC_fus*surface_fuselage))/surface_total
    xac = 0.275*MAC # keep attention that it is an estimation the table don't give the value for this very low AR

    x_tmp = leading[np.argmin(abs(y - yac))]
    xac = x_tmp+xac 
    xac = ((xac_fus*surface_fuselage) + (xac_wing*surface_wing))/surface_total

    return MAC_fus, yac_fus, xac_fus, MAC_wing, yac_wing, xac_wing, MAC, yac, xac

def equivalentWing(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    b_wing, AR_wing, sweep_beta_wing, sweep_beta_wing, c_root_wing, taper_ratio_wing, sweep_quarter_wing, c_tip_wing, y_wing, leading_edge_wing, trailing_edge_wing, quarter_line_wing, c, h = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, weight)
    surface_wing_ideal, AR, taper_wing, croot, ctip, chord_middle, sweep_LE_wing, sweep_beta = wingGeometryIDEAL(Cl, weight, sweep_quarter_wing)
    b_fus, AR_fuselage, sweep_beta_fus, c_root_fus, taper_ratio_fus, sweep_quarter_fus, c_tip_fus, y_fus, leading_edge_fus, trailing_edge_fus, quarter_line_fus = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    surface_wing_ideal, surface_fuselage, surface_wing, surface_total = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    c = np.concatenate(([c_root_fus], c))
    h = np.concatenate(([b_fus/2], h))
    s = np.zeros(len(c)-1)
    for i in range(len(c)-1):
        s[i] = h[i]*(c[i]+c[i+1])/2
    
    b = span_max#(b_fus+b_wing)
    Sw = 0
    for i in range(len(s)):
        Sw += (span_max/surface_total)*((s[i]*c[i]) + (s[i]*c[i+1]))
    
    Cre = 0
    for i in range(len(s)):
        Cre += (2/Sw)*(c[i]*s[i])

    Cte = 0
    for i in range(len(s)):
        Cte += (2/Sw)*(c[i+1]*s[i])
    
    lambda_E = Cte/Cre

    sweep_LE = (1/surface_total)*(sweep_LE_fus*surface_fuselage + sweep_LE_wing*surface_wing)*(np.pi/180)
    sweep_quarter = np.arctan(np.tan(sweep_LE) + ((Cre/(2*span_max) * lambda_E -1)))

    print("taper ratio equi wing", lambda_E)
    Mgc = (2*Cre/3)*(1+lambda_E+lambda_E**2)/(1+lambda_E)

    y_mgc = (span_max/6)*(1+2*lambda_E)/(1+lambda_E) 
    x_mgc = y_mgc * np.tan(sweep_LE) + ((c[0]/2) - (Cre/2))
    #addx =( c_root_fus/2) - (Cre/2) #/!\
    return Cre, Cte, lambda_E, sweep_LE, sweep_quarter, Mgc, y_mgc, x_mgc


def plotAllWing(wing_plot, sweep_LE_fus, sweep_quarter_wing, Cl, weight):
    if wing_plot == False:
        return
    Cre, Cte, lambda_E, sweep_LE, sweep_quarter, Mgc, y_mgc, x_mgc = equivalentWing(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    _, _, _, _, _, _, _, y_fus, leading_edge_fus, trailing_edge_fus, quarter_line_fus = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight) 
    plt.plot(y_fus, leading_edge_fus, color='blue')
    plt.plot(y_fus, trailing_edge_fus, color='green')
    #plt.plot(y_fus, quarter_line_fus, color='red')

    _, _, _, _, _, _, _, _, y_wing, leading_edge_wing, trailing_edge_wing, quarter_line_wing, _, _ = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, weight) 
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

    MAC_fus, yac_fus, xac_fus, MAC_wing, yac_wing, xac_wing, MAC, yac, xac = getMAC(Cl, sweep_LE_fus, sweep_quarter_wing, weight)

    plt.scatter(yac_fus,  xac_fus, color='red') #+ (yac_fus*np.tan(sweep_LE_fus*(np.pi/180)))
    plt.scatter(yac_wing,  xac_wing, color='orange')
    plt.scatter(yac, xac, color='blue')
    plt.scatter(y_mgc, x_mgc + Mgc*0.25, color='black')

    leading_edge_fus_x = np.interp(yac_fus, y_fus, leading_edge_fus)
    leading_edge_wing_x = np.interp(yac_wing , y_wing + cabin_width/2, leading_edge_wing)
    plt.plot((yac_fus, yac_fus), (leading_edge_fus_x, leading_edge_fus_x + MAC_fus), color='red')
    plt.plot((yac_wing, yac_wing), (leading_edge_wing_x + leading_edge_fus[-1], leading_edge_wing_x+ leading_edge_fus[-1] + MAC_wing), color='orange')
    
    
    #plt.plot((span_max/2, span_max/2), ((span_max/2)*np.tan(sweep_LE)+ ((trailing_edge_fus[0])/2) - (Cre/2), (span_max/2)*np.tan(sweep_LE) + Cte + ((trailing_edge_fus[0])/2) - (Cre/2)), color='black')
    #plt.plot((0, span_max/2), (((trailing_edge_fus[0])/2) - (Cre/2) , (span_max/2)*np.tan(sweep_LE)+ ((trailing_edge_fus[0])/2) - (Cre/2)), color='black')
    #plt.plot((0, span_max/2), (((trailing_edge_fus[0])/2) - (Cre/2) + Cre, (span_max/2)*np.tan(sweep_LE) + Cte + ((trailing_edge_fus[0])/2) - (Cre/2)), color='black')

    #plt.plot((-span_max/2, -span_max/2), ((span_max/2)*np.tan(sweep_LE)+ ((trailing_edge_fus[0])/2) - (Cre/2), (span_max/2)*np.tan(sweep_LE) + Cte + ((trailing_edge_fus[0])/2) - (Cre/2)), color='black')
    #plt.plot((0, -span_max/2), (((trailing_edge_fus[0])/2) - (Cre/2) , (span_max/2)*np.tan(sweep_LE)+ ((trailing_edge_fus[0])/2) - (Cre/2)), color='black')
    #plt.plot((0, -span_max/2), (((trailing_edge_fus[0])/2) - (Cre/2) + Cre, (span_max/2)*np.tan(sweep_LE) + Cte + ((trailing_edge_fus[0])/2) - (Cre/2)), color='black')

    #plt.plot((y_mgc, y_mgc), (x_mgc, x_mgc+Mgc), color='black')
    
    plt.xlabel('$Y$')
    plt.ylabel('$X$')
    # Fixer l'échelle des axes
    plt.axis('equal')
    plt.show() 
    return



def get_Lift_and_drag(Cl, delta, sweep_LE_fus, sweep_quarter_wing, weight):
    surface_wing_ideal, AR, taper_wing, croot, ctip, chord_middle, sweep_leading, sweep_beta_wing = wingGeometryIDEAL(Cl, weight, sweep_quarter_wing)
    b, AR_fuselage, sweep_beta_fus, c_root, taper_ratio, sweep_quarter, c_tip, _, _, _, _ = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    surface_wing_ideal, surface_fuselage, surface_wing, surface_total = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    AoA = np.linspace(-10, 10, 21) * (np.pi/180)
    cl_alpha_fus, _, alpha_L0_fus, _, _ = getAirfoilFus()
    cl_alpha_wing, _, _, _, _ = getAirfoilWing()
    Cl_wing, Cl_wing_0, Cd_wing, Cl_max_wing, alpha_L0_wing, a_wing = wingCL(Cl,sweep_LE_fus, sweep_quarter_wing, weight)
    Cl_fuselage, Cl_fus_0, Cd_fuselage, Cl_max_fus, a_fus = fuselageCL(Cl, sweep_LE_fus, sweep_quarter_wing, weight)

    percent_wing = Cl_wing_0/(Cl_wing_0 + Cl_fus_0)
    percent_fus = Cl_fus_0/(Cl_wing_0 + Cl_fus_0)
    CL_alfa = ((a_wing*percent_wing) + (a_fus*percent_fus))/1
    # --- total lift computation --- #
    """ Cl_tot = np.zeros(len(AoA))

    ## percentage of lift of the wing and fuselage
    CL_alfa = ((a_wing*percent_wing) + (a_fus*percent_fus))/1
    alpha_L0 = ((alpha_L0_wing*percent_wing) + (alpha_L0_fus*percent_fus))/1
    #print("CL_alfa percentage", CL_alfa)
    #print("alpha_L0 percentage", alpha_L0)
    ## ponderation of the lift of the wing and fuselage
    #CL_alfa = ((a_wing*surface_wing) + (a_fus*surface_fuselage))/surface_wing_ideal
    #alpha_L0 = ((alpha_L0_wing*surface_wing) + (alpha_L0_fus*surface_fuselage))/surface_wing_ideal
    #print("CL_alfa ponderation", CL_alfa)
    #print("alpha_L0 ponderation", alpha_L0)
    ## lift of the wing and fuselage compute with the ponderation of cl and sweep
    cl_alpha_tot = ((cl_alpha_wing*surface_wing) + (cl_alpha_fus*surface_fuselage))/surface_total
    sweep_beta_tot = (sweep_beta_wing*surface_wing + sweep_beta_fus*surface_fuselage)/surface_total
    k_tot = (beta * cl_alpha_tot)/(2*np.pi)

    #alpha_L0 = ((alpha_L0_wing*surface_wing) + (alpha_L0_fus*surface_fuselage))/surface_total
    CL_alfa = ((2*np.pi)/((2/(beta*AR)) + np.sqrt((1/((k_tot * np.cos(sweep_beta_tot))))**2 + ((2/(beta * AR))**2) )))/beta
    # print("CL_alfa ponderation composante", CL_alfa)
    # print("alpha_L0 ponderation composante", alpha_L0)
    
    for i in range(len(AoA)):    
        Cl_tot[i] = CL_alfa*(AoA[i]  - alpha_L0) """
    
    Cl_tot = Cl_wing + Cl_fuselage 
    
    AoA_L0 = np.interp(0, Cl_tot, AoA) * (180 / np.pi)

    Cl_tot0 = np.interp(0, AoA*180/np.pi, Cl_tot)
    

    cl_max = (Cl_max_wing + Cl_max_fus)/2
    

    # --- total drag computation --- #
    AR_cd = AR + winglet(AR)
    Cd_induce = ((Cl_tot**2)/(np.pi* AR_cd)) * (1+delta)
    Cd_tot = np.zeros(len(AoA))
    cd0 = 0.012 # in cruise
    Cd_tot = Cd_induce + cd0 
    Cd_tot0 = np.interp(0, AoA, Cd_tot)

    
    """
    CL_CD = Cl_tot / Cd_tot
    # Trouver l'indice du maximum de CL/CD
    max_index = np.argmax(CL_CD)
    # Extraire le CL correspondant au maximum de CL/CD
    CL_max = Cl_tot[max_index]
    CL_CD_max = CL_CD[max_index]
    print(CL_max)
    print(CL_CD_max)
    """
    
    AoA = AoA * (180 / np.pi)  
    with open("data_lift_cruise.txt", "w") as file:
        for angle, cl, cd in zip(AoA, Cl_tot, Cd_tot):  # Boucle sur les valeurs des arrays
            file.write(f"{angle:.2f} {cl:.4f} {cd:.4f}\n")  # Formatage propre"
            
    #print(((air_density(alti)[0] * (true_airspeed_at_altitude(alti, M)**2))/(air_density(0)[0] * true_airspeed_at_altitude(0, 0.45)**2)))
    
            
    return Cl_tot0, Cd_tot0, cl_max, AoA_L0, Cl_tot, Cd_tot, AoA, cd0, CL_alfa

def getClAlfa(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    surface_wing_ideal, AR, taper_wing, croot, ctip, chord_middle, sweep_LE_wing, sweep_beta = wingGeometryIDEAL(Cl, weight, sweep_quarter_wing)
    surface_wing_ideal, surface_fuselage, surface_wing, surface_total = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    sweep_LE_tot = (sweep_LE_wing*surface_wing + sweep_LE_fus*surface_fuselage)/surface_total
    sweep_LE_tot = sweep_LE_tot*(np.pi/180)
    
    cl_alpha_wing, _, _, _, _ = getAirfoilWing()
    cl_alpha_fus, _, _, _, _ = getAirfoilFus()
    cl_alpha = (cl_alpha_wing*surface_wing + cl_alpha_fus*surface_fuselage)/surface_total
    
    b_wing, AR_wing, sweep_beta_wing, sweep_beta_tot_wing, c_root_wing, taper_ratio_wing, sweep_quarter_wing, c_tip_wing, y_wing, leading_edge_wing, trailing_edge_wing, quarter_line_wing, c_wing, h_wing = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, weight)
    b_fus, AR_fuselage, sweep_beta_fus, c_root_fus, taper_ratio_fus, sweep_quarter_fus, c_tip_fus, y_fus, leading_edge_fus, trailing_edge_fus, quarter_line_fus = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    taper_ratio_tot = (taper_ratio_wing*surface_wing + taper_ratio_fus*surface_fuselage)/surface_total

    sweep_quarter = np.arctan(np.tan(sweep_LE_tot) + (4/AR) * (((1-taper_ratio_tot)/(1+taper_ratio_tot)) * (0 - 0.25)))
    sweep_beta_tot = np.arctan2(np.tan(sweep_quarter), beta) 
    
    k = (beta * cl_alpha)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR)) + np.sqrt((1/((k * np.cos(sweep_beta_tot))))**2 + ((2/(beta * AR))**2) )))/beta
    return a


def plotLiftDrag(lift_and_drag_plots, Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    if lift_and_drag_plots == False:
        return
    Cl_tot0, Cd_tot0, Cl_max, AoA_L0, Cl_tot, Cd_tot, AoA, Cd_tot0, CL_alfa =  get_Lift_and_drag(Cl, delta, sweep_LE_fus, sweep_quarter_wing, weight)
    plt.figure(figsize=(8,5))
    plt.plot(Cl_tot, Cl_tot/Cd_tot)
    plt.xlabel('$CL$')
    plt.ylabel('$CL/CD$')
    plt.show()

    plt.figure(figsize=(10, 6),dpi=300)
    plt.plot(Cd_tot, Cl_tot)
    plt.xlabel('$C_D$ [-]')
    plt.ylabel('$C_L$ [-]')
    plt.savefig("/Users/antoinevanhoye/Documents/M1/PI/integrated_project/Airfoils/drag_polar.pdf", dpi=300)
    plt.show()

    plt.figure(figsize=(8,5))
    plt.plot(AoA, Cl_tot)
    plt.xlabel('$AoA$ [rad]')
    plt.ylabel('$CL$')
    plt.show()

    plt.figure(figsize=(8,5))
    plt.plot(AoA, Cd_tot)
    plt.xlabel('$AoA$ [rad]')
    plt.ylabel('$CD$')
    plt.show()
    return


def wingMaxthickness(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    _, CL, _, _, _, _ = wingCL(Cl,sweep_LE_fus, sweep_quarter_wing, weight)
    b, AR_wing, sweep_beta, sweep_beta_tot, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line, c, h = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, weight)
    M_star = 1.15-(CL/(4*(np.cos(sweep_quarter)**2)))
    t_bar_over_c =  (3/(10*M)) * np.cbrt((1/(M*np.cos(sweep_quarter))) - (M*np.cos(sweep_quarter))) * (1-(((5 + (M**2)*(np.cos(sweep_quarter)**2))/(5 + (M_star**2)))**3.5))**(2/3)
    t_root = t_bar_over_c*c_root
    t_tip = t_bar_over_c*c_tip 
    
    return t_root, t_tip, t_bar_over_c


def wingFuelvolume(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    surface_wing_ideal, surface_fuselage, surface_wing, surface_total = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    b, AR_wing, sweep_beta, sweep_beta_tot, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line, c, h = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, weight)
    t_root, t_tip, t_bar_over_c = wingMaxthickness(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    tau = 1
    V_fuel = 0.54*((surface_wing**2)/b) * t_bar_over_c * ((1+(taper_ratio*np.sqrt(tau))+((taper_ratio**2)*tau))/(1+taper_ratio)**2)
    return V_fuel


def wingSurfaceWetted(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    surface_wing_ideal, surface_fuselage, surface_wing, surface_total = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    b, AR_wing, sweep_beta, sweep_beta_tot, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line, c, h = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, weight)
    t_root, t_tip, t_bar_over_c = wingMaxthickness(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    S_wet_wing = 2*surface_wing*(1+ (1/4) * ((t_bar_over_c + (t_bar_over_c*taper_ratio))/(1+taper_ratio)))

    b, AR_fuselage, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    S_wet_fus = 2*surface_fuselage*(1+ (1/4) * ((0.18 + (0.18*taper_ratio))/(1+taper_ratio)))

    S_wet_tot = S_wet_wing + S_wet_fus
    return S_wet_wing, S_wet_fus, S_wet_tot

def stallVelocity(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    surface_wing_ideal, surface_fuselage, surface_wing, surface_total = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    b, AR_wing, sweep_beta, sweep_beta_tot, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line, c, h = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, weight)
    _, _, Cl_max, _ , _, _, _, _, _= get_Lift_and_drag(Cl, delta, sweep_LE_fus, sweep_quarter_wing, weight)
    rho_sl, T = air_density(0)

    Vs = np.sqrt((weight/surface_total) * (2/rho_sl) * (1/(1.133*Cl_max)))
    
    Cl_max0 = 2.1 #Amos' value   #2 * np.cos(sweep_quarter)
    W0 = weight_empty #landing weight
    Vs0 = np.sqrt((W0/surface_total) * (2/rho_sl) * (1/(1.133*Cl_max0)))
    return Vs, Vs0

def getReynold(altitude, c):
    rho, T = air_density(altitude)
    U_inf = true_airspeed_at_altitude(altitude, M)
    
    p_atmo = 99333      #[Pa]
    T = 12.0 + 273.15   #[k]
    rho = p_atmo/(287*T)
    mu = 1.716e-5 * (T/273.15)**(3/2) * ((273.15 + 110.4)/(T + 110.4)) # Sutherland's law
    Re = (rho * U_inf * c) / mu
    return Re

def getHighLiftDevice(Cl, sweep_LE_fus, sweep_quarter_wing, weight):
    surface_wing_ideal, surface_fuselage, surface_wing, surface_total = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    surface_wing_ideal, AR, taper_wing, croot, c_tip_wing, chord_middle, sweep_LE_wing, sweep_beta = wingGeometryIDEAL(Cl, weight, sweep_quarter_wing)
    _, _, _, c_root_fus, _, _, _, _, _, _, _ = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    taper_ratio = c_tip_wing/c_root_fus
    CL_max_tot = get_Lift_and_drag(Cl, delta, sweep_LE_fus, sweep_quarter_wing, weight)[2]
    sweep_LE_tot = ((sweep_LE_wing*surface_wing + sweep_LE_fus*surface_fuselage)/surface_total)*(np.pi/180)
    sweep_HL = np.arctan(np.tan(sweep_LE_tot) + (4/AR) * (((1-taper_ratio)/(1+taper_ratio)) * (0 - 0.70)))
    #print(sweep_HL*(180/np.pi))
    delta_CL_max = 2.1 - CL_max_tot  #2.1 Amos' value And 1.9 for landing configuration
    delta_cl_max = delta_CL_max* (1/(0.8* np.cos(sweep_HL)))
    return delta_cl_max

def main():
    Cl, sweep_LE_fus, sweep_quarter_wing, weight = 0.45, 50.0, 29.0,  469794.4881175733
    
    surface_wing_ideal, surface_fuselage, surface_wing, surface_total = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    _, _, _, _, _, _, _, c_tip_wing, _, _, _, _, _, _ = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, weight)
    surface_wing_ideal, AR, taper_wing, croot, c_tip_wing, chord_middle, sweep_LE_wing, sweep_beta = wingGeometryIDEAL(Cl, weight, sweep_quarter_wing)
    surface_wing_ideal, AR = getSurface_And_AR(Cl, weight)
    _, _, _, c_root_fus, _, _, _, _, _, _, _ = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    delta_cl_max = getHighLiftDevice(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    print("\n-------------- Total values --------------\n")
    print(f"New AR = {AR:.3f} [-]")
    print(f"Cl used = {Cl:.2f} [-]\n")
    print(f"Total area = {surface_total:.2f} [m^2]")
    print(f"Surface of fuselage = {surface_fuselage:.2f} [m^2]")
    print(f"Surface of wing = {surface_wing:.2f} [m^2]")
    print(f"Surface of wing ideal = {surface_wing_ideal:.2f} [m^2]")
    print(f"Compressibility parameter beta: {beta:.3f} [-]")
    print(f"Taper ratio: {c_tip_wing/c_root_fus:.3f} [-]")
    print(f"High lift device delta clmax: {delta_cl_max:.3f} [-]\n")
    print(f"Cl alfa: {getClAlfa(Cl, sweep_LE_fus, sweep_quarter_wing, weight):.3f} [rad^-1]\n")
    
    b, AR_wing, sweep_beta, sweep_beta_tot, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line, c, h = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, weight)
    Cl_wing, Cl_wing_0, Cd_wing, Cl_max_wing, alpha_L0, a_wing = wingCL(Cl,sweep_LE_fus, sweep_quarter_wing, weight)
    cl_alpha_wing, cl_max_wing, alpha_l0_wing, CD_wing, cm_wing = getAirfoilWing()
    print("\n-------------- wing values --------------\n")
    print(f"\nAR wing: {AR_wing:.3f} [-] \nCL_w0 wing = {Cl_wing_0:.3f} [-]\n")
    print(f"Wing lift percentage, {float((Cl_wing_0*(surface_wing/surface_total))/Cl):.3f} ")
    print(f"Chord at wing root: {croot:.3f} [m]\nChord at wing tip: {c_tip_wing:.3f} [m]\nChord at middle wing: {chord_middle:.3f} [m]")
    print(f"Taper ratio: {taper_ratio:.3f} [-]")
    print(f"sweep quater: {sweep_quarter*(180/np.pi):.3f} [°]")
    print(f"Wing lift coefficient derivative: {a_wing:.3f}")
    print(f"Alpha_L0: {alpha_L0*(180/np.pi):.3f}")
    print(f"Cl max wing: {Cl_max_wing:.3f} [-]")
    print(f"cl alpha wing airfoil: {cl_alpha_wing:.3f} [rad^-1]\n")
    
    

    b, AR_fuselage, sweep_beta, c_root, taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    Cl_fuselage, Cl_fuselage_0, Cd_fuselage, Cl_max_fus, a_fus = fuselageCL(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    cl_alpha_fus, cl_max_fus, alpha_l0_fus, CD_fus, cm_fus = getAirfoilFus()
    print("\n-------------- fuselage values --------------\n")
    print(f"\nAR fuselage: {AR_fuselage:.3f} [-]\nCL_w0 fuselage = {Cl_fuselage_0*surface_fuselage/surface_wing_ideal:.3f} [-]\n")
    print(f"Cord at fuselage root: {c_root:.3f} [m]\nCorde at fuselage tip: {c_tip:.3f} [m]")
    print(f"Taper ratio: {taper_ratio:.3f} [-]")
    print(f"sweep quater: {sweep_quarter*(180/np.pi):.3f}")
    print(f"Fuselage lift coefficient derivative: {a_fus:.3f}")
    print(f"Cl max fuselage: {Cl_max_fus:.3f} [-]")
    print(f"cl alpha fuselage airfoil: {cl_alpha_fus:.3f} [rad^-1]\n")
    
    MAC_fus, yac_fus, xac_fus, MAC_wing, yac_wing, xac_wing, MAC, yac, xac = getMAC(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    #Mgc, y_mgc, Cwre = getMAC2(AR, sweep_LE_fus, sweep_LE_wing)
    print("\n-------------- MAC values --------------\n")
    print(f"MAC fus: {MAC_fus:.3f} [m]\nYac fus: {yac_fus:.3f} [m]\nXac fus: {xac_fus:.3f} [m]\n")
    print(f"MAC wing: {MAC_wing:.3f} [m]\nYac wing: {yac_wing:.3f} [m]\nXac wing: {xac_wing:.3f} [m]\n")
    print(f"MAC: {MAC:.3f} [m]\nYac: {yac:.3f} [m]\nXac: {xac:.3f} [m]\n")
    #print(f"MGC: {Mgc:.3f} [m]\nY_MGC: {y_mgc:.3f} [m]\nCwre: {Cwre:.3f} [m]\n")
    

    delta = 0.005 #graph slide 61 lecture 6 aerodynimics
    lift_coef, drag_coef, CL_max, AoA_L0, cl, _, aoa, Cd_tot0, CL_alfa = get_Lift_and_drag(Cl, delta, sweep_LE_fus, sweep_quarter_wing, weight)
    print("\n-------------- Lift and drag --------------\n")
    print(f"\n CL = {lift_coef:.3f}[-] \n CD = {drag_coef:.5f}[-] \n")
    print(f"Lift to drag ratio: {(lift_coef/drag_coef):.3f} [-]")
    print(f"Cl max: {CL_max:.3f} [-]")
    print(f"Lift coefficient derivative CL_alfa: {CL_alfa:.3f} [rad^-1]")
    print(f"CD0: {Cd_tot0:.3f} [-]\n")
    
    
    t_root, t_tip,t_bar_over_C = wingMaxthickness(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    print("\n-------------- Other value --------------\n")
    print(f"Thickness root: {t_root:.3f}, thickness tip: {t_tip:.3f}")
    print(f"tbar over c: {t_bar_over_C:.3f} [-]\n")
    #print(f"Mean thinckness: {(t_root+t_tip)/2:.3f}\n")
    
    V_fuel = wingFuelvolume(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    print(f"Fuel volume in wing: {V_fuel:.3f} [m^3]\n")

    _, _, Swetted = wingSurfaceWetted(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    print(f"Surface wetted: {Swetted:.3f} [m^2]\n")
    
    Vs, Vs0 = stallVelocity(Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    print(f"Stall velocity: {Vs:.3f} [m/s]\nStall velocity in approach config: {Vs0:.3f} [m/s]\n")
    


    AoA_root,_ = getCalageAngle(Cl, sweep_LE_fus, sweep_quarter_wing, weight) #right values ????
    print(f"Setting angle: {AoA_root*(180/np.pi):.3} [°]")
    print(f"AoA zero lift: {AoA_L0:.3f} [°]")
    
    Re = getReynold(alti, MAC)
    print(f"Re_mac: {Re:.3f} [-]")

    plotAllWing(wing_plot, sweep_LE_fus, sweep_quarter_wing, Cl, weight)
    plotLiftDrag(lift_and_drag_plots, Cl, sweep_LE_fus, sweep_quarter_wing, weight)
    return


if __name__ == '__main__':
    main()
