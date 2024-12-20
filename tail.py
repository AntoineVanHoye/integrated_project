# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 16:14:27 2024

@author: diego
"""


"""
Group member:
- Leo Di Puma           - Louis Belboom
- Emile de Lamalle      - Diego Schyns
- Amos David            - Antoine Van Hoye
- Guillaume Vaillot
"""

import numpy as np
import matplotlib as plt
from scipy.integrate import trapz
###Variables a changer car pas sur des données#############


S = 182.56 #m² area of total surface Wing+ Fuselage
MAC = 3.171 #m mean chord of the line 
b = 29 #m span
M=0.9 # nombre de mach attention!!!!!!!!!!!!!!!!!! pas le même que les ailes A CHANGER
R=287 
M_C=0.9

###############CHOIX#####################

####NACA
AoA_tail=0
S_tot_tail=30# [m^2] 30% of the wing area
gamma_h= 30*np.pi/180 # [rad] dihedral angle

AR_tail=2 #   Aspect ratio of the V-tail
lambda_tail=0.95 #Taper ratio of the V-tail
sweep_angle_tail_c4=30 #[°] sweep angle at chord quarter
c_root_tail=1.5 #[m] chord à la racine de la tail 
sweep_leading_tail=45 #[°]
K=0.66
#########PLOT###############
cl_plot=True

  #[m]bras de levier du centre aerodynamique de la tail 
 #hypothèses 

#NACA 0012 for vertical parts of the tail
def geomtail():
    # tail doit être dans les 10% avec la technique V_mc  
    #calcul des surfaces d une tail conventionnel 
    
    ####ROSKAM
    # S_h=V_h*S*MAC/(x_h)
    # S_v=V_v*S*b/x_v
    # S_tot_tail=np.sqrt(S_v**2+S_h**2)
    #calcul du butterfly angle
    # gamma_h=np.arctan(S_v/S_h)
    
    ####NACA
    
    S_h=S_tot_tail*(np.cos(gamma_h))**2
    
    S_v=S_tot_tail*(np.sin(gamma_h))**2
    
    #Calcul de l AR de la tail

    b_tail=np.sqrt(AR_tail*S_tot_tail) #m longueur de la queue
    c_root_tail=2*S_tot_tail/(b_tail*(1+lambda_tail*np.cos(gamma_h)))
    c_tip_tail=lambda_tail*c_root_tail
    chord_tail=2/3*c_root_tail*((lambda_tail**2+lambda_tail+1)/(lambda_tail+1))
    
    return  c_tip_tail, chord_tail,gamma_h,S_h,S_v,c_root_tail ,b_tail,S_tot_tail


c_tip_tail,chord_tail,gamma,S_h,S_v,c_root_tail,b_tail,S_tot_tail_real=geomtail()


def getMACTail():
    beta=np.sqrt(1-M**2)
    sweep_leading_tail1 = sweep_leading_tail*((np.pi)/180)
    sweep_quarter_tail = np.arctan(np.tan(sweep_leading_tail1) + (4/AR_tail) * (((1-lambda_tail)/(1+lambda_tail)) * (0 - 0.25)))
    sweep_trailing_tail = np.arctan(np.tan(sweep_quarter_tail) + (4/AR_tail) * (((1-lambda_tail)/(1+lambda_tail)) * (0.25 - 1)))
    sweep_beta_tail = np.arctan2(np.tan(sweep_quarter_tail), beta)
    y_tail = np.linspace(0, b_tail/2, 10)
    quarter_line_tail = np.zeros(len(y_tail))
    leading_edge_tail = np.zeros(len(y_tail))
    trailing_edge_tail = np.zeros(len(y_tail))
    for i in range(len(y_tail)):
        quarter_line_tail[i] = (np.tan(sweep_quarter_tail))*y_tail[i] + (0.25*c_root_tail)
        leading_edge_tail[i] = (np.tan(sweep_leading_tail1))*y_tail[i] 
        trailing_edge_tail[i] = (np.tan(sweep_trailing_tail))*y_tail[i] + c_root_tail
# -- tail -- #
    c_tail = trailing_edge_tail - leading_edge_tail
    MAC_tail = (2/S_tot_tail) * trapz(c_tail**2, y_tail) #numerical integration via method of trapez
    cy_tail = c_tail*y_tail
    yac_wing = (2/S_tot_tail) * trapz(cy_tail, y_tail)
    xac_wing = MAC_tail*0.2
    return xac_wing, yac_wing,MAC_tail,sweep_beta_tail, sweep_quarter_tail


xac, yac, MACtail, sweep_beta_tail, sweep_quarter_tail=getMACTail()


def NACA():
    beta=np.sqrt(1-M**2)
    cl_alpha=(1.5461-1.0909)/(15-10)* (180/np.pi)
    alpha_L0=0
    cl_max= 1.4551
    C_d=0.011781 #at Aoa=0°
    c_tip_tail, chord_tail,gamma_h,S_h,S_v,c_root_tail ,b_tail,S_tot_tail= geomtail()
    xac_wing, yac_wing,MAC_tail,sweep_beta_tail,sweep_quarter_tail=getMACTail()
    k = (beta * cl_alpha)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR_tail)) + np.sqrt((1/((k * np.cos(sweep_beta_tail))))**2 + ((2/(beta * AR_tail))**2) )))/beta
    AoA = np.linspace(-10, 10, 200) * ((np.pi)/180)
    CL_t = np.zeros(len(AoA))
    CL_tN = np.zeros(len(AoA))
    for i in range(len(AoA)):    
        CL_t[i] = a*(AoA[i] - alpha_L0)
        if AoA[i] <= 0:
            if AoA[i+1] >= 0:
                CL_t0 = (CL_t[i] + a*(AoA[i+1] - alpha_L0))/2
    CL_t0 = a*((AoA_tail * (np.pi/180)) - alpha_L0) # choose of AoA of the wing 
    

    CL_t_max = np.cos(sweep_quarter_tail) * 0.95 * ((cl_max + cl_max)/2)
    
    return CL_t,CL_t0, CL_t_max

CL_t,CL_t0,CL_t_max=NACA()



def changergamma_h(CL):
    beta=np.sqrt(1-M**2)
    cl_alpha=(1.5461-1.0909)/(15-10)* (180/np.pi)
    alpha_L0=0
    cl_max= 1.4551
    C_d=0.011781 #at Aoa=0°
    k = (beta * cl_alpha)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR_tail)) + np.sqrt((1/((k * np.cos(sweep_beta_tail))))**2 + ((2/(beta * AR_tail))**2) )))/beta
    
    alpha_root =CL/a
    return alpha_root
    

    
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

# Function to calculate the true airspeed at a given altitude
def true_airspeed_at_altitude(altitude):
    T = air_density(altitude)[1]
    a = np.sqrt(gamma1 * R * T)  # [m/s] Speed of sound
    v = M_C * a  # [m/s] Aircraft velocity
    return v

def getReynold(altitude):
    rho, T = air_density(altitude)
    U_inf = true_airspeed_at_altitude(altitude)
    p_atmo = 99333      #[Pa]
    T = 12.0 + 273.15   #[k]
    rho = p_atmo/(287*T)
    mu = 1.716e-5 * (T/273.15)**(3/2)*((273.15 + 110.4)/(T + 110.4)) # Sutherland's law
    Re = (rho * U_inf * c_root_tail) / mu
    return Re
gamma1=1.4
altitude=12500
rho, T=   air_density(altitude)
Re=getReynold(altitude)
speed=true_airspeed_at_altitude(altitude)
    
Need_CL=-55387.2468/(1/2*rho*speed**2*S_h)






print("CL requis",Need_CL)
    
alpha_root=changergamma_h(Need_CL)

print("angle d attaque",alpha_root*180/np.pi)
    

    
    
    
    
    
    