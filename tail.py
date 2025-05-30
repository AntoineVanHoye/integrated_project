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
"""
import numpy as np
import matplotlib as plt
from scipy.integrate import trapz
import matplotlib.pyplot as plt

###### Ligne 33 et 34 bras de levier, Ligne 204 force à changer dans la formule. Verifier l'angle d'attaque dans le terminal

M=0.9 # nombre de mach attention!!!!!!!!!!!!!!!!!! pas le même que les ailes A CHANGER
R=287
M_C=0.9

###############CHOIX#####################
AR_tail=3.5 # choix chelou mais ancien projet prenait ça
sweep_leading_tail=3 #[°] roskam
incidence_angle=0 # comme dans roskam
lambda_tail=1 #roskam
L_VT=4.35047143507#[ft] distance frome the tail quarter chord (25% of the mean chord length measure back ffrom the leading of the mean chord) to the wing quarter chord
L_HT=L_VT#[ft] distance frome the tail quarter chord (25% of the mean chord length measure back ffrom the leading of the mean chord) to the wing quarter chord
c_HT=0.720625 # Horizontal tail volume coefficient for jet transport Daniel P.Raymer
c_VT=0.072875 # Vertical tail volume coefficient for jet transport Daniel P.Raymer
b_w= 20#[ft] wing span
S_w=61.75# [ft^2]wing area####NACA
##PLOT###############
cl_plot=True

print("bras de levier vertical de la tail :", L_VT)
print("bras de levier horizontal de la tail :", L_HT)

  #[m]bras de levier du centre aerodynamique de la tail
 #hypothèses

#NACA 0012 for vertical parts of the tail
# def geomtail():
#     # tail doit être dans les 10% avec la technique V_mc  
#     #calcul des surfaces d une tail conventionnel
   
#     ####ROSKAM
#     # S_h=V_h*S*MAC/(x_h)
#     # S_v=V_v*S*b/x_v
#     # S_tot_tail=np.sqrt(S_v**2+S_h**2)
#     #calcul du butterfly angle
#     # gamma_h=np.arctan(S_v/S_h)
   
#     ####NACA
   
#     S_h=S_tot_tail*(np.cos(gamma_h))**2
   
#     S_v=S_tot_tail*(np.sin(gamma_h))**2
   
#     #Calcul de l AR de la tail

#     b_tail=np.sqrt(AR_tail*S_tot_tail) #m longueur de la queue
#     c_root_tail=2*S_tot_tail/(b_tail*(1+lambda_tail*np.cos(gamma_h)))
#     c_tip_tail=lambda_tail*c_root_tail
#     chord_tail=2/3*c_root_tail*((lambda_tail**2+lambda_tail+1)/(lambda_tail+1))
   
#     return  c_tip_tail, chord_tail,gamma_h,S_h,S_v,c_root_tail ,b_tail,S_tot_tail

   
def geomtail():
    ###TOUS LES CHOIX DE GEOMETRIE ON ETE FAIT POUR RESPECTER LES TABLES 8.13 ET 8.14 DE ROSKAM PART 2 EN HORIZONTAL ET VERTICAL
   
    C_bar_w=S_w/b_w
    S_VT=c_VT*b_w*S_w/L_VT
    S_HT=c_HT*C_bar_w*S_w/L_HT
    S_tot_tail=np.sqrt(S_VT**2+S_HT**2)
    gamma_h=np.arctan(S_VT/S_HT)
    b_tail=np.sqrt(AR_tail*S_tot_tail) #m longueur de la queue
    c_root_tail=2*S_tot_tail/(b_tail*(1+lambda_tail*np.cos(gamma_h)))
    c_tip_tail=lambda_tail*c_root_tail
    C_bar_tail=2/3*c_root_tail*((lambda_tail**2+lambda_tail+1)/(lambda_tail+1))
   

    return  c_tip_tail, C_bar_tail,gamma_h,S_HT,S_VT,c_root_tail ,b_tail,S_tot_tail
b_tail = _,_,_,_,_,_,b_tail,_ = geomtail()
print("span =",b_tail)

c_tip_tail,chord_tail,gamma,S_h,S_v,c_root_tail,b_tail,S_tot_tail_real=geomtail()
print("Surface horizontale de la tail:",S_h)
print("L'angle de la tail:",gamma*180/np.pi)
print("La surface totale de la tail:",S_tot_tail_real)
print("Corde à la racine de la tail et dans notre cas aussi la MAC:",c_root_tail)
gamma_deg=gamma*180/np.pi  

def getMACTail():
    c_tip_tail, C_bar_tail,gamma_h,S_HT,S_VT,c_root_tail ,b_tail,S_tot_tail=geomtail()
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
    c_tail = C_bar_tail
    MAC_tail = c_tail
   
    yac_wing = (2/S_tot_tail) * c_tail * 1/2 * (b_tail/2)**2
    xac_wing = MAC_tail*0.2
    return xac_wing, yac_wing,MAC_tail,sweep_beta_tail, sweep_quarter_tail


xac, yac, MACtail, sweep_beta_tail, sweep_quarter_tail=getMACTail()


def NACA():
    AoA_tail=0
    beta=np.sqrt(1-M**2)
    cl_alpha=(0.5413-0.6572)/(5-6)/np.pi*180
    alpha_L0=0
    cl_max= 1.423
    C_d=0.014 #at Aoa=0°
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
   
   
    return a,CL_t,CL_t0, CL_t_max,cl_alpha

a,CL_t,CL_t0,CL_t_max,cl_alpha=NACA()

print("CLapha is", a)

def changergamma_h(CL):
    beta=np.sqrt(1-M**2)
    cl_alpha=(1.5461-1.0909)/(15-10)* (180/np.pi)
    alpha_L0=0
    cl_max= 1.423
    C_d=0.014#at Aoa=0°
    k = (beta * cl_alpha)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR_tail)) + np.sqrt((1/((k * np.cos(sweep_beta_tail))))**2 + ((2/(beta * AR_tail))**2) )))/beta
   
    alpha_root =CL/a
    return alpha_root,a
   

   
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
S_h_en_mettrecarré=S_h#/10.7639
Need_CL=-305990.760081231/(1/2*rho*speed**2*S_h_en_mettrecarré)



print(Need_CL)
   
alpha_root,a=changergamma_h(Need_CL)

print("Angle d'attaque de la tail:",alpha_root*180/np.pi)
   
def stab_directional():  
    S_fs = ...  # Surface du fuselage
    length_fus = ...  # Longueur du fuselage
    Sb = ...  # Surface de l'aile
    hf1 = ...  # Hauteur fuselage à l'avant
    hf2 = ...  # Hauteur fuselage à l'arrière
    bf1 = ...  # Largeur fuselage à l'avant
    bf2 = ...  # Largeur fuselage à l'arrière
    hf_max = ...  # Hauteur maximale du fuselage
    lcg = ...  # Longueur jusqu'au centre de gravité
    L_f=... # Distance entre centre de gravité et centre aerodynamique de la Fin (les autres le font sans faire varier ce point mais nous appparemment on devrait parce les voyageurs se déplacent selon l'axe y, selon Leo )
    S=... #surface des ailes??????????
    b=... #span des ailes ???????????
   
    #J'initialise la CN_beta_tot à zéro, il est important de noter que pour la stabilité directionnel, les ailes le fuselage et l'empennage vertical on chacun une contribution séparée.
    CN_beta_tot=0
   

   
    # Calcul de CN_beta_f, formule donnée dans le cours qui requiert divers données géométriques de l'avion que nous n'avons pas encore il me semble
    K_beta = 0.3 * (lcg / length_fus) + 0.75 * (hf_max / length_fus) - 0.105
    CN_beta_fuselage = -K_beta * (S_fs * length_fus / Sb) * ((hf1 / hf2) ** 0.5) * ((bf2 / bf1) ** (1/3))
   
    # La contribution des ailes semble être une correlation empirique en fonction de comment sont positionee les ailes sur l'avion
    CN_beta_w=0.012 #{High, mid, low}-mounted wing effect = {-0.017,0.012,0.024}
       
    #Je reprend mes données calculée en haut, notamment la pente de CL de mon airfoil
    a,CL_t,CL_t0, CL_t_max=NACA()
    c_tip_tail, chord_tail,gamma_h,S_h,S_v,c_root_tail ,b_tail,S_tot_tail=geomtail()
   
    #C' est ici que les romains s'empoignèrent, c'est la partie la plus tactique du code, slide 56 du cours de conceptual
    #design de Noels, on remarque que le moment C_N dépend de C_LF de la fin. Le vent battant qui arrive sur l'empennage
    #vertical n'est autre que l'analogie à l'angle d'attaque qui arrive sur l'empennage horizontal, fais un petit dessin
    #pour mieux te le représenter. C'est pourquoi lorsqu'on dérive CN en fonction de beta on sors de la dérivée les
    #constantes et on ne dérive alors plus que CLF en fonction de beta. La on se demande qu'est ce que c'est CLF en fonction
    # de beta ?? C'est la que l'analogie à l'empennage horizontal nous dit que c'est comme si l'airfoil de la fin possedait
    #un angle d'attaque beta, c'est pourquoi la dérivée de CLF n'est autre que la pente de notre CL de l'airfoil qui est "a"
    CN_beta_fin=a*S_v*L_f/(S*b)
   
    CN_beta_tot=CN_beta_fin+CN_beta_w+CN_beta_fuselage
   
    if (CN_beta_tot>0):
        print ("je ne sais pas encore si C_N_beta doit être plus grand ou plus petit que zero, pas clair dans le cours")
    else:
        print ("je ne sais pas encore si C_N_beta doit être plus grand ou plus petit que zero, pas clair dans le cours")

"""



import numpy as np
import matplotlib as plt
from scipy.integrate import trapz
import matplotlib.pyplot as plt
from wings import detSurfac
from wings import air_density
from wings import true_airspeed_at_altitude


M=0.9 
beta=np.sqrt(1-M**2)
R=287
M_C=0.9

###############PLOT OF THE HORIZONTAL TAIL#####################
tail_plot = False
####################################


speed = true_airspeed_at_altitude(12500,0.9)
rho = air_density(12500)[0]

l_cabin = 10.1
l_cockpit = 2.01

sweep_leading_tail=37*np.pi/180 #choice
incidence_angle=0 #choice
taper_ratio=0.6 #choice
tail_pos = l_cabin + l_cockpit + 1 #choice


def surf_tail():
    surf_hor_tail = 35
    surf_vert_tail = 19
    return surf_hor_tail, surf_vert_tail

def geomtail():
    surf_hor_tail = surf_tail()[0]
    surf_vert_tail = surf_tail()[1]
    surf_tot_tail = np.sqrt(surf_vert_tail**2+surf_hor_tail**2)
    c_root_tail = 3.5 #choice
    
    span_hor = 2*surf_hor_tail/(taper_ratio*c_root_tail+c_root_tail)
    span_vert = 2*surf_vert_tail/(taper_ratio*c_root_tail+c_root_tail)
    gamma_h = np.arctan(surf_vert_tail/surf_hor_tail)
    span_vert = np.tan(gamma_h) * span_hor
    span_tot = np.sqrt(span_hor**2+span_vert**2)
    
    """
    gamma_h = np.arctan(np.sqrt(surf_vert_tail/surf_hor_tail))
    surf_tot_tail = surf_hor_tail/np.cos(gamma_h)**2
    span_tot = 2*surf_tot_tail/(c_root_tail*(1+taper_ratio))
    span_hor = np.cos(gamma_h)*span_tot
    span_vert = np.sin(gamma_h)*span_tot
    """
    AR_h = span_hor**2/surf_hor_tail
    AR = span_tot**2/surf_tot_tail
    
    c_tip_tail = taper_ratio*c_root_tail
    
    y = np.linspace(0,span_hor/2,100)
    y1 = np.linspace(0,span_vert/2,100)
    c = np.linspace(c_root_tail, c_tip_tail, len(y))
    quarter_line_tail = np.zeros(len(y))
    leading_edge_tail = np.zeros(len(y))
    trailing_edge_tail = np.zeros(len(y))
    sweep_quarter_tail = np.arctan(np.tan(sweep_leading_tail) + (4/AR) * (((1-taper_ratio)/(1+taper_ratio)) * (0 - 0.25)))
    sweep_trailing_tail = np.arctan(np.tan(sweep_quarter_tail) + (4/AR) * (((1-taper_ratio)/(1+taper_ratio)) * (0.25 - 1)))
    for i in range(len(y)):
        quarter_line_tail[i] = (np.tan(sweep_quarter_tail))*y[i] + (0.25*c_root_tail)
        leading_edge_tail[i] = (np.tan(sweep_leading_tail))*y[i]
        trailing_edge_tail[i] = (np.tan(sweep_trailing_tail))*y[i] + c_root_tail
    
    MAC_tail = 2/surf_tot_tail*trapz(c**2, y)
    cy = c*y
    yac = (2/surf_vert_tail) *trapz(cy,y1)
    xac = MAC_tail*0.2 + yac * np.tan(sweep_leading_tail)
    return  c_root_tail, span_hor, span_vert,AR_h, AR,gamma_h, surf_tot_tail, MAC_tail,yac,xac

c_root_tail,span_hor,span_vert,AR_h, AR,gamma_h, surf_tot_tail, MAC_tail,yac,xac = geomtail()


def plotTail():
    if tail_plot == False:
        return
    
    y = np.linspace(0,span_hor/2,100)

    quarter_line_tail = np.zeros(len(y))
    leading_edge_tail = np.zeros(len(y))
    trailing_edge_tail = np.zeros(len(y))
    sweep_quarter_tail = np.arctan(np.tan(sweep_leading_tail) + (4/AR) * (((1-taper_ratio)/(1+taper_ratio)) * (0 - 0.25)))
    sweep_trailing_tail = np.arctan(np.tan(sweep_quarter_tail) + (4/AR) * (((1-taper_ratio)/(1+taper_ratio)) * (0.25 - 1)))
    for i in range(len(y)):
        quarter_line_tail[i] = (np.tan(sweep_quarter_tail))*y[i] + (0.25*c_root_tail)
        leading_edge_tail[i] = (np.tan(sweep_leading_tail))*y[i]
        trailing_edge_tail[i] = (np.tan(sweep_trailing_tail))*y[i] + c_root_tail
    plt.plot(y, leading_edge_tail)
    plt.plot(y, trailing_edge_tail)
    plt.show()

    return

plotTail()


def getSweepTail():
    
    c_root_tail,span_hor,span_vert,AR_h, AR,gamma_h, surf_tot_tail, MAC_tail,yac,xac = geomtail()
    sweep_quarter_tail = np.arctan(np.tan(sweep_leading_tail) + (4/AR_h) * (((1-taper_ratio)/(1+taper_ratio)) * (0 - 0.25)))
    sweep_trailing_tail = np.arctan(np.tan(sweep_quarter_tail) + (4/AR_h) * (((1-taper_ratio)/(1+taper_ratio)) * (0.25 - 1)))
    sweep_beta_tail = np.arctan2(np.tan(sweep_quarter_tail), beta)

    return sweep_beta_tail, sweep_quarter_tail

def LiftCurveSlope():

    c_root_tail,span_hor,span_vert,AR_h, AR,gamma_h, surf_tot_tail, MAC_tail,yac,xac = geomtail()
    cl_alpha=(0.5413-0.6572)/(5-6)/np.pi*180 # NACA 0012
    #cl_alpha = (-0.0799-1.0135)/(-5-5)/np.pi*180 # NACA 4412
    sweep_beta_tail,sweep_quarter_tail=getSweepTail()
    k = (beta * cl_alpha)/(2*np.pi)
    a = ((2*np.pi)/((2/(beta*AR_h)) + np.sqrt((1/((k * np.cos(sweep_beta_tail))))**2 + ((2/(beta * AR_h))**2) )))/beta

    return a

a=LiftCurveSlope()


def need_CL(force):
    surf_hor_tail = surf_tail()[0]
    Need_CL=-force/(1/2*rho*speed**2*surf_hor_tail)
    return Need_CL

def setting_angle(force):
    a = LiftCurveSlope()
    CL = need_CL(force)
    #alpha_L0 = -4.35*np.pi/180
    #alpha_root =CL/a - alpha_L0/a
    alpha_root =CL/a
    return alpha_root

force =-70000

def main():
    c_root_tail,span_hor,span_vert,AR_h, AR,gamma_h, surf_tot_tail, MAC_tail,yac,xac = geomtail()
    sweep_beta_tail,sweep_quarter_tail=getSweepTail()
    a = LiftCurveSlope()
    alpha_root = setting_angle(force)
    print("True dihedral angle is", np.arctan(surf_tail()[1]/surf_tail()[0])*180/np.pi)
    print("--------------------------ROOT CHORD--------------------------------------------")
    print("The root chord of the tail is",c_root_tail,"m or",c_root_tail*3.28084,"ft")
    print("--------------------------LIFT COEFFICIENT--------------------------------------------")
    print("The lift coefficient of the tail is",need_CL(force))
    print("--------------------------SPANS--------------------------------------------")
    print("The horizontal span of the tail is",span_hor,"m or",span_hor*3.28084,"ft")
    print("The vertical span of the tail is",span_vert,"m or",span_vert*3.28084,"ft")
    print("The total span of the tail is",np.sqrt(span_hor**2+span_vert**2),"m or",np.sqrt(span_hor**2+span_vert**2)*3.28084,"ft")
    print("--------------------------AR--------------------------------------------")
    print("The aspect ratio of the horizontal tail is",AR_h)
    print("The aspect ratio of the vertical tail is",span_vert**2/surf_tail()[1])
    print("The total aspect ratio of the tail is",AR)
    print("--------------------------SWEEP QUARTER--------------------------------------------")
    print("The sweep quarter chord of the tail is",sweep_quarter_tail*180/np.pi,"degrees")
    print("--------------------------DIHEDRAL ANGLE--------------------------------------------")
    print("The dihedral angle gamma of the tail is",gamma_h*180/np.pi,"degrees")
    print("--------------------------AREAS--------------------------------------------")
    print("The total area of the tail is",surf_tot_tail,"m^2 (",surf_tot_tail*3.28084**2,"ft^2")# and it represents",surf_tot_tail/surf_tot*100,"% of the total lifting surface")
    print("The horizontal area of the tail is",surf_tail()[0],"m^2 (",surf_tail()[0]*3.28084**2,"ft^2")# and it represents",surf_tail()[0]/surf_tot*100,"% of the total lifting surface")
    print("The vertical area of the tail is",surf_tail()[1],"m^2 (",surf_tail()[1]*3.28084**2,"ft^2")# and it represents",surf_tail()[1]/surf_tot*100,"% of the total lifting surface")
    print("--------------------------MAC--------------------------------------------")
    print("The MAC of the tail is",MAC_tail,"m or",MAC_tail*3.28084,"ft")
    print("--------------------------AC--------------------------------------------")
    print("The y position of the aerodynamic center of the tail is",yac,"m or",yac*3.28084,"ft.")
    print("The x position of the aerodynamic center of the tail is",xac,"m or",xac*3.28084,"ft in the local axis and",tail_pos+xac,"m and",(tail_pos+xac)*3.28084,"ft with respect to the nose of the aircraft.")
    print("--------------------------LIFT CURVE SLOPE--------------------------------------------")
    print("The lift curve slope of the tail is", a)
    print("--------------------------SETTING ANGLE--------------------------------------------")
    print("The needed setting angle of the tail is", alpha_root*180/np.pi,"degrees.")
    return

if __name__ == "__main__":
    main()






