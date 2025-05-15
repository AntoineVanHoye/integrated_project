import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import trapz
from wings import fusGeometry as fusgeom
from wings import wingGeometry as winggeom
from sympy import symbols, Eq, solve
from wings import air_density
from wings import true_airspeed_at_altitude
from weight import get_weight
from wings import wingFuelvolume
from wings import getMAC
from wings import fusGeometry
from wings import detSurfac
from wings import wingCL
from wings import fuselageCL
from wings import getAirfoilFus
from wings import getAirfoilWing
from wings import getClAlfa
from wings import wingCL
from wings import wingGeometryIDEAL
from tail import LiftCurveSlope
from tail import geomtail
from tail import surf_tail
from tail import LiftCurveSlope
from wings import getSurface_And_AR
from wings import wingGeometry


##################################################################
######GENERAL QUANTTITES
##################################################################

rho = air_density(12500)[0]
speed = true_airspeed_at_altitude(12500,0.9)

delta = 0.005
b = 28.95
l_fus = 16.8
width_cabin = 9
l_cabin = 10.1
l_cockpit = 2.01 
span_wings = b-width_cabin 
l_aft = l_fus - l_cabin - l_cockpit


##################################################################
######COMPUTATION OF THE CM0_TOT
##################################################################

Cm0_wing = getAirfoilWing()[4]
Cm0_fus = getAirfoilFus()[4]

def Cm0(Cm0_airfoil_fus,Cm0_airfoil_wing, Cl ,sweep_LE_fus, sweep_quarter_wing, force):
    surface_wing_ideal = 119.05
    surf_fus  = 121.07
    surf_wing = 71.11
    MAC_fus = 14.293
    MAC_wing = 3.706
    MAC_tot = 11.098
    y_wing = np.array([0, 1.10833333, 2.21666667, 3.325,      4.43333333, 5.54166667,6.65,       7.75833333, 8.86666667, 9.975     ] )
    leading_wing = np.array([0,         0.68183035, 1.36366069, 2.04549104, 2.72732139, 3.40915173,4.09098208, 4.77281242, 5.45464277, 6.13647312])
    trailing_wing =np.array([4.77878354, 5.1907293 , 5.60267505, 6.01462081, 6.42656656, 6.83851232,7.25045808, 7.66240383 ,8.07434959, 8.48629534])
    y_fus = np.array([0, 0.5, 1,  1.5 ,2,  2.5 ,3, 3.5 ,4, 4.5])
    leading_fus = np.array([0,        0.5958768,  1.19175359, 1.78763039, 2.38350719 ,2.97938398,3.57526078, 4.17113757 ,4.76701437 ,5.36289117])
    trailing_fus = np.array([16.8 ,16.8, 16.8, 16.8 ,16.8 ,16.8 ,16.8, 16.8, 16.8 ,16.8])
    c_fus = trailing_fus - leading_fus
    c_wing = trailing_wing - leading_wing
    Cm0_wing = (2/(surf_wing*MAC_wing)) * trapz(Cm0_airfoil_wing*c_wing**2, y_wing)
    Cm0_fus = (2/(surf_fus*MAC_fus)) * trapz(Cm0_airfoil_fus*c_fus**2, y_fus)
    M0_wing = Cm0_wing*(1/2*rho*speed**2*MAC_wing*surf_wing)
    M0_fus = Cm0_fus*(1/2*rho*speed**2*MAC_fus*surf_fus)
    M0_tot = M0_wing + M0_fus
    Cm0_tot = M0_tot/(1/2*rho*speed**2*MAC_tot*surface_wing_ideal)
    
    #Cm0_wing = trapz(Cm0_airfoil_wing*c_wing**2, y_wing)/(1/2*MAC_tot*surface_wing_ideal)
    #Cm0_fus =trapz(Cm0_airfoil_fus*c_fus**2, y_fus)/(1/2*MAC_tot*surface_wing_ideal)
    #Cm0_tot = Cm0_wing + Cm0_fus
    #Cm0_tot = 2/(surface_wing_ideal*MAC_tot)*( trapz(Cm0_airfoil_wing*c_wing**2, y_wing) + trapz(Cm0_airfoil_fus*c_fus**2, y_fus))
    #Cm0_tot = (Cm0_wing*surf_wing*MAC_wing + Cm0_fus*surf_fus*MAC_fus)/(surface_wing_ideal*MAC_tot)
    #Cm0_tot = (Cm0_wing*surf_wing + Cm0_fus*surf_fus)/surf_tot
    
    return Cm0_tot,Cm0_fus,Cm0_wing


##################################################################
######GEOMETRICAL POSITIONS OF IMPORTANT POINTS
##################################################################

z_AC_tot = 0
z_CG_tot = 0
z_CG_motors = 1.7

##################################################################
######TAIL QUANTITIES
##################################################################

c_root_tail,span_hor_tail,span_vert_tail,AR_h_tail, AR_tail,surf_vert_tail, surf_tot_tail, MAC_tail,y_AC_tail,x_AC_tail_local = geomtail()

hor_tail_surf = surf_tail()[0]
a1 = LiftCurveSlope()
x_AC_tail = l_cabin + l_cockpit + x_AC_tail_local + 1
l_tail = MAC_tail

z_AC_tail = y_AC_tail

##################################################################
######CONFIGURATION SETTING
##################################################################

config =2
fuel = 2

##################################################################
######CG POSITION
##################################################################

def passengers(i): 

    passengers_weight = 8*215
    if i == 1 : #if no passengers 
       passengers_pos = 1 #arbitrary
       passengers_weight = 0
    if i == 2 : #if passengers are as close as possible to the nose
       passengers_pos = l_cockpit
    if i == 3 : #if passengers are as far as possible from the nose 
       passengers_pos = l_cockpit + l_cabin
    
    return passengers_weight, passengers_pos

def CG_position(i,d, Cl, sweep_LE_fus, sweep_quarter_wing, force): 
    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(Cl, sweep_LE_fus, sweep_quarter_wing, force)
    fus_weight, aft_weight, wing_weight, land_gear_weight,motors_weight,nacelle_weight,APU_weight,enginst_weight,instr_weight,hydr_syst_weight,furn_weight,air_cond_weight,payload_weight, ops_weight,elec_syst_weight,surf_cont_weight,_,_,_= get_weight()
    _,_,_,_,_,_,chord_tip_fus,_,_,_,_ = fusGeometry(Cl, sweep_LE_fus, sweep_quarter_wing, force)
    surface_wing_ideal, AR, taper_wing, croot, ctip, chord_middle, sweep_LE_wing, sweep_beta = wingGeometryIDEAL(Cl, force, sweep_quarter_wing)

    #force = CL(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing, Cl, sweep_LE_fus, sweep_quarter_wing, force)[0] ## Boucle infinie en utilisant ceci

    sweep_angle_wing = sweep_LE_wing*np.pi/180
    sweep_angle_fus = sweep_LE_fus*np.pi/180

    wing_pos = (l_fus - chord_tip_fus) + MAC_wing*0.3 + y_AC_wing*np.tan(sweep_angle_wing)
    
    fus_pos = (l_cabin + l_cockpit)*0.5 #0.2*MAC_fus + y_AC_fus*np.tan(sweep_angle_fus)

    aft_pos = l_cockpit + l_cabin + l_aft*0.3

    APU_pos = l_cockpit + l_cabin

    hydr_pos = 0.75*wing_pos + 0.25 *x_AC_tail

    payload_pos = 0.5*(l_cockpit + l_cabin) + 0.5*l_cockpit

    ops_pos = l_cockpit 

    land_gear_pos = 9.703*0.7 + 0.3*0.673
    surf_cont_pos = (l_fus - chord_tip_fus) + MAC_wing*0.35 + y_AC_wing*np.tan(sweep_angle_wing)

    instr_pos = 0.35*l_cockpit

    furn_pos = 0.4*l_fus

    air_cond_pos = l_cockpit + l_cabin
    
    motors_pos = l_cockpit + l_cabin + 1

    nacelle_pos = motors_pos

    enginst_pos = motors_pos

    elec_syst_pos = 0.75*l_fus/2 + 0.25*motors_pos
  
    available_fuel_vol = 15 * 1000 #wingFuelvolume(AR, sweep_LE_fus, sweep_LE_wing, force)*1000 
    if d == 1 : #no fuel
        fuel_weight = 0
        fuel_pos = 1
        pourc_wings = 0
    if d == 2 : #full of fuel
        vol_fuel = 25032.322462
        fuel_weight = vol_fuel*0.8*2.20462
        #pourc_wings = available_fuel_vol/vol_fuel
        pourc_wings = 1
        #fuel_pos = wing_pos*pourc_wings + (1-pourc_wings)*l_fus*0.44
        fuel_pos = wing_pos*pourc_wings 

    passengers_weight = passengers(i)[0]
    passengers_pos = passengers(i)[1] 

    """
    print("fuel:",fuel_pos,"m and",fuel_pos*3.28084,"ft ->",fuel_pos/l_fus*100)
    print("wing:",wing_pos, "m and",wing_pos*3.28084,"ft ->", wing_pos/l_fus *100)
    print("fus:",fus_pos, "m and",fus_pos*3.28084,"ft ->", fus_pos/l_fus *100)
    print("aft:",aft_pos, "m and",aft_pos*3.28084,"ft ->", aft_pos/l_fus *100)
    print("APU:",APU_pos, "m and", APU_pos*3.28084,"ft ->", APU_pos/l_fus *100)
    print("hydr:",hydr_pos, "m and", hydr_pos*3.28084,"ft ->", hydr_pos/l_fus *100)
    print("payload:",payload_pos, "m and",payload_pos*3.28084,"ft ->", payload_pos/l_fus *100) 
    print("ops:",ops_pos, "m and",ops_pos*3.28084,"ft ->", ops_pos/l_fus *100)
    print("land_gear:",land_gear_pos, "m and", land_gear_pos*3.28084,"ft ->", land_gear_pos/l_fus *100)
    print("surf_cont:",surf_cont_pos, "m and", surf_cont_pos*3.28084,"ft ->", surf_cont_pos/l_fus *100)
    print("instr:",instr_pos, "m and", instr_pos*3.28084,"ft ->", instr_pos/l_fus *100)
    print("furn:",furn_pos, "m and",furn_pos*3.28084,"ft ->", furn_pos/l_fus *100)
    print("air_cond:",air_cond_pos, "m and",air_cond_pos*3.28084,"ft ->", air_cond_pos/l_fus *100)
    print("motors:",motors_pos, "m and",motors_pos*3.28084,"ft ->", motors_pos/l_fus *100)
    print("nacelle:",nacelle_pos, "m and", nacelle_pos*3.28084,"ft ->", nacelle_pos/l_fus *100)
    print("enginst:",enginst_pos, "m and",enginst_pos*3.28084,"ft ->", enginst_pos/l_fus *100)
    print("elec_syst:",elec_syst_pos, "m and",elec_syst_pos*3.28084,"ft ->", elec_syst_pos/l_fus *100)
    """
    
    total_mom = (wing_weight*wing_pos) + (fus_weight*fus_pos) + (land_gear_weight*land_gear_pos) + (surf_cont_weight*surf_cont_pos) + (instr_weight*instr_pos) + (elec_syst_weight*elec_syst_pos) + (furn_weight*furn_pos) + (air_cond_weight*air_cond_pos) + (passengers_weight*passengers_pos) + (motors_weight*motors_pos) + (fuel_pos*fuel_weight) + (aft_pos*aft_weight) + (nacelle_pos*nacelle_weight) + (APU_pos*APU_weight) + (enginst_pos*enginst_weight) + (hydr_pos*hydr_syst_weight) + (payload_pos*payload_weight) + (ops_pos*ops_weight) 
    total_weight = wing_weight + fus_weight + land_gear_weight + surf_cont_weight + instr_weight + elec_syst_weight + furn_weight + air_cond_weight + passengers_weight + motors_weight + fuel_weight + aft_weight + nacelle_weight + APU_weight + enginst_weight + hydr_syst_weight + payload_weight + ops_weight 
    position = total_mom/total_weight

    position = 10
    return position, pourc_wings, motors_pos/MAC_tot, total_weight 


##################################################################
######EQUILIBRIUM IN PITCH
##################################################################


#tail volume ratio effectivness 
def tail_eff(i,d, Cl, sweep_LE_fus, sweep_quarter_wing, force):
    surface_wing_ideal, surf_fus, surf_wing, surf_tot = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, force)
    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(Cl, sweep_LE_fus, sweep_quarter_wing, force)

    x_CG_tot = CG_position(i,d, Cl, sweep_LE_fus, sweep_quarter_wing, force)[0]
    V_T = hor_tail_surf * (x_AC_tail - x_CG_tot)/(surf_tot* MAC_tot)
    return V_T

def prop_force(): 
    A = 2*1.32**2 *np.pi
    m_dot = rho*A*speed*0.453592/32.2
    Fp = m_dot*speed
    return Fp



#i : configuration numerotation
def CL(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing, Cl, sweep_LE_fus, sweep_quarter_wing, force): 

    surface_wing_ideal, surf_fus, surf_wing, surf_tot = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, force)
    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(Cl, sweep_LE_fus, sweep_quarter_wing, force)
    x_AC_tot = 8.38
    L_tot, L_T = symbols('L_tot L_T')
    V_T = tail_eff(i,d, Cl, sweep_LE_fus, sweep_quarter_wing, force)
    T =38778.1#134926.83#take-off : 134926.83 32374.06 
    x_CG_tot = CG_position(i,d, Cl, sweep_LE_fus, sweep_quarter_wing, force)[0]
    weight = CG_position(i,d, Cl, sweep_LE_fus, sweep_quarter_wing, force)[3]*9.81*0.453592 + passengers(i)[0]*9.81*0.453592
    
    Cm0_tot = Cm0(Cm0_airfoil_fus,Cm0_airfoil_wing, Cl ,sweep_LE_fus, sweep_quarter_wing, force)[0]
    M0 = Cm0_tot*((1/2)*rho*(speed**2)*MAC_tot*surface_wing_ideal)
    x_CG_motors = CG_position(i,d, Cl, sweep_LE_fus, sweep_quarter_wing, force)[2]*MAC_tot
    Fp = prop_force()

    #translation equilibrium 
    eq1 = Eq(L_tot + L_T - weight + Fp, 0)
    
    #rotation equilibrium 
    eq2 = Eq(M0 + L_tot*(x_CG_tot-x_AC_tot)-L_T*(x_AC_tail - x_CG_tot)-T*(z_CG_motors - z_CG_tot) + Fp*(x_CG_tot-x_CG_motors),0) 

    solution = solve((eq1, eq2),(L_tot, L_T))
    L_tot = float(solution[L_tot])
    L_T = float(solution[L_T])
    
    #if abs(L_tot/force) > 0.1 or abs(force/L_tot) > 0.1:
    #    L_tot, L_T = CL(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing, Cl, sweep_LE_fus, sweep_quarter_wing, float(L_tot))
    return L_tot,L_T

def boucleForce(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing, Cl, sweep_LE_fus, sweep_quarter_wing):
    force = 200000.0
    L_tot = 1.0
    while abs(L_tot - force)/L_tot > 0.05:
        L_tot, L_T = CL(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing, Cl, sweep_LE_fus, sweep_quarter_wing, force)
        force = L_tot
    return L_tot, L_T


##################################################################
######LONGITUDINAL STABILITY
##################################################################


def long_stat_stab_cruise(i,d, Cl, sweep_LE_fus, sweep_quarter_wing, force): #in the pitching plane
    surface_wing_ideal, surf_fus, surf_wing, surf_tot = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, force)
    a = 3
    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(Cl, sweep_LE_fus, sweep_quarter_wing, force)
    #check the stability
    x_AC_tot = 8.56
    #neutral point : position of the cg in order to have the derivative equals 0
    x_CG_tot = CG_position(i,d, Cl, sweep_LE_fus, sweep_quarter_wing, force)[0]
    engines_pos = CG_position(i,d, Cl, sweep_LE_fus, sweep_quarter_wing, force)[2]
    deps = 0#2*3.17/(np.pi*7.04)
    Fp = prop_force()
    eta = 0.9
    #hn = x_AC_tot/MAC_tot + V_T*a1_over_a*(1- deps) + dalpha_prop * Fp * (engines_pos-x_CG_tot)/MAC_tot#- 0.5*fus_width**2 * fus_length/(S_wing*a*MAC_wing)  #position of the neutral point  

    #page 416 Raymer
    q = 1/2*rho*speed**2
    hn = (a * x_AC_tot/MAC_tot +eta*a1*hor_tail_surf/surf_tot*(1-deps)*x_AC_tail/MAC_tot+Fp/(q*surf_tot)*(1-deps)*engines_pos)/(a + eta*hor_tail_surf/surf_tot*a1*(1-deps)+Fp/(q*surf_tot))

    Kn = hn - x_CG_tot/MAC_tot 
    return Kn, hn


##################################################################
######CG POSITIONS RANGE
##################################################################


def get_CG(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing,Kn, Cl, sweep_LE_fus, sweep_quarter_wing, force): 
    surface_wing_ideal, surf_fus, surf_wing, surf_tot = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, force)
    a = 3
    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(Cl, sweep_LE_fus, sweep_quarter_wing, force)
    engines_pos = CG_position(i,d, Cl, sweep_LE_fus, sweep_quarter_wing, force)[2]
    deps = 0#2*3.17/(np.pi*7.04)
    x_AC_tot = 8.3241
    V_T = tail_eff(i,d, Cl, sweep_LE_fus, sweep_quarter_wing, force)
    dalpha_prop = 1 - deps
    Fp = prop_force()
    eta = 0.9
    q = 1/2*rho*speed**2
    hn = (a * x_AC_tot/MAC_tot + eta*a1*hor_tail_surf/surf_tot*(1-deps)*x_AC_tail/MAC_tot+Fp/(q*surf_tot)*(1-deps)*engines_pos)/(a + eta*hor_tail_surf/surf_tot*a1*(1-deps)+Fp/(q*surf_tot))

    #x_CG = (hn*MAC_tot) - (Kn*MAC_tot)
    x_CG1 = (hn*MAC_tot) - (0.05*MAC_tot)
    x_CG2 = (hn*MAC_tot) - (0.15*MAC_tot)
    return x_CG1, x_CG2


##################################################################
######DIRECTIONAL STABILITY
##################################################################


def interpolation(x1, y1, x2, y2, x3):
    t=(x3-x1)/(x2-x1) 
    y3 = y1 + t * (y2 - y1)
    
    return y3

def dir_stat_stab_cruise(CG_position, Cl, sweep_LE_fus, sweep_quarter_wing, force):  
    surface_wing_ideal, surf_fus, surf_wing, surf_tot = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, force)
    
    """
    hf1 = (interpolation(0.2481847, 0.129909, 0.2632646, 0.1303305, 0.25)+interpolation(0.2491559,0.04703807 , 0.2613637,0.0475382 , 0.25))*l_fus  # forward fuselage height
    hf2 = (interpolation(0.7477641, 0.04391509, 0.7610847, 0.04077155, 0.75)+interpolation(0.7403715, 0.05064632, 0.7543387, 0.04952601, 0.75))*l_fus # rear fuselage height
    bf1 = 5.881743321*3.28084# forward fuselage width
    bf2 = 9  # rear fuselage width
    K_beta = 0.3 * (x_CG / l_fus) + 0.75 * (hf_max / l_fus) - 0.105
    CN_beta_fuselage = -K_beta * (surf_fus*l_fus/(surf_wing*span_wings))*((hf1/hf2)**0.5)*((bf2/bf1)**(1/3))
    CN_beta_fin=a*surf_vert_tail*L_f/(surf_wing*span_wings)
    CN_beta_w=0.012 #{High, mid, low}-mounted wing effect = {-0.017,0.012,0.024}
    """

    hf_max = 0.179*16.8*3.28084  # maximum fuselage height
    x_CG = CG_position *3.28084
    b = 28.95 *3.28084
    L_f= x_AC_tail*3.28084 - x_CG # distance between center of gravity and aerodynamic center of the tail
    surf_fus = surf_fus*3.28084**2
    fus_vol = 224.98*3.28084**3
    h = 1.762*3.28084 #mean fuselage depth
    w = 7.564*3.28084 #mean fuselage width
    surface_wing_ideal = surface_wing_ideal*3.28084**2
    surf_vert_tail = surf_tail()[1]*3.28084**2
    V_T = surf_vert_tail* L_f/(surface_wing_ideal* b)
    _, AR = getSurface_And_AR(Cl, force)
    mean_chord = AR/b
    sweep_quarter_wing = sweep_quarter_wing*np.pi/180
    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(Cl, sweep_LE_fus, sweep_quarter_wing, force)
    x_AC_tot = 8.37895
    CN_beta_fuselage = -1.3*fus_vol/(surface_wing_ideal*b)*(h/w)

    CN_beta_w = Cl**2 *(1/(4*np.pi*AR) - np.tan(sweep_quarter_wing)/(np.pi*AR*(AR+4*np.cos(sweep_quarter_wing)))*(np.cos(sweep_quarter_wing) - AR/2 - AR**2/(8*np.cos(sweep_quarter_wing)) + 6*np.abs(x_CG - x_AC_tot*3.28084)/(mean_chord*3.28084) * np.sin(sweep_quarter_wing/AR)))

    z_w = 0
    a = LiftCurveSlope()
    term = 0.724 + ((3.06*(surf_vert_tail/surface_wing_ideal))/(1 + np.cos(sweep_quarter_wing))) + 0.4 *(z_w/hf_max) + 0.009*AR
    CN_beta_fin = a*term*V_T
    
    CN_beta_tot=CN_beta_fin+CN_beta_w+CN_beta_fuselage
    
    return CN_beta_fin, CN_beta_fuselage, CN_beta_w, CN_beta_tot


##################################################################
######LATERAL STABILITY
##################################################################

def lat_stat_stab_cruise(CG_position,dihedral_angle, Cl,sweep_LE_fus, sweep_quarter_wing, force): 

    surface_wing_ideal, surf_fus, surf_wing, surf_tot = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, force)
    b, AR_wing, sweep_beta, sweep_beta_tot, c_root, wings_taper_ratio, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line, c, h = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, force)
    gamma = dihedral_angle*np.pi/180
    CL_wings, CL_w0, CD_wing, CL_max, alpha_L0, CL_alpha_wings = wingCL(Cl,sweep_LE_fus, sweep_quarter_wing, force)

    CL_beta_dihedral = -0.25 *CL_alpha_wings * gamma* (2 * (1 + 2*wings_taper_ratio)/(3*(1+wings_taper_ratio)))

    CL_beta_wings_fus = 0 #because middle mounted wing

    hf_max = 0.179*16.8*3.28084  # maximum fuselage height
    surface_wing_ideal = surface_wing_ideal*3.28084**2
    surf_vert_tail = surf_tail()[1]*3.28084**2
    _, AR = getSurface_And_AR(Cl, force)
    z_w = 0
    term = 0.724 + ((3.06*(surf_vert_tail/surface_wing_ideal))/(1 + np.cos(sweep_quarter_wing))) + 0.4 *(z_w/hf_max) + 0.009*AR
    a = LiftCurveSlope()
    x_CG = CG_position *3.28084
    b = 28.95 *3.28084
    CL_beta_tail = -a*term*surf_vert_tail/surface_wing_ideal*z_AC_tail/b
    
    #x = wings_taper_ratio
    graph_value  = np.interp(wings_taper_ratio, [0.5, 1], [-0.21, -0.19]) #Nicolai's book, page 590

    CL_beta_wing_sweep = graph_value * CL_w0
    CL_beta_wings = CL_beta_wing_sweep + CL_beta_dihedral

    CL_beta_tot = CL_beta_wings_fus + CL_beta_wings + CL_beta_tail

    return CL_beta_tot, CL_beta_wings, CL_beta_wings_fus, CL_beta_tail

##################################################################
######PRINT OF STABILITY
##################################################################


def main():
    Cl, sweep_LE_fus, sweep_quarter_wing, dihedral_angle = 0.45, 50.0, 29.0, 0

    force, force_tail = boucleForce(config, fuel,Cm0_fus,Cm0_wing, Cl, sweep_LE_fus, sweep_quarter_wing) #CL(config, fuel,Cm0_fus,Cm0_wing, Cl, sweep_LE_fus, sweep_quarter_wing)

    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(Cl, sweep_LE_fus, sweep_quarter_wing, force)
    CN_beta_fin, CN_beta_fuselage, CN_beta_w, CN_beta_tot = dir_stat_stab_cruise(CG_position(config,fuel, Cl, sweep_LE_fus, sweep_quarter_wing, force)[0], Cl, sweep_LE_fus, sweep_quarter_wing, force)
    CL_beta_tot, CL_beta_wings, CL_beta_wings_fus, CL_beta_tail = lat_stat_stab_cruise(CG_position(config,fuel, Cl, sweep_LE_fus, sweep_quarter_wing, force)[0],dihedral_angle, Cl,sweep_LE_fus, sweep_quarter_wing, force)

    print("--------------------------AERODYNAMIC CENTER--------------------------------------------")
    print("The aerodynamic center is positioned at",x_AC_tot,"from the nose of the airplane, which represents",x_AC_tot*100/l_fus,"% of the total length.")
    print("----------------------------------------------------------------------")

    print("--------------------------CM0--------------------------------------------")
    print("The blended wing body has a Cm0_tot of",Cm0(Cm0_fus,Cm0_wing, Cl ,sweep_LE_fus, sweep_quarter_wing, force)[0])
    print("The wings have a pitching moment (up) of",Cm0(Cm0_fus,Cm0_wing, Cl ,sweep_LE_fus, sweep_quarter_wing, force)[2],". It represents",Cm0(Cm0_fus,Cm0_wing, Cl ,sweep_LE_fus, sweep_quarter_wing, force)[2]/Cm0(Cm0_fus,Cm0_wing, Cl ,sweep_LE_fus, sweep_quarter_wing, force)[0]*100,"% of the pitching moment of the airplane.")
    print("The fuselage has a pitching moment (up) of",Cm0(Cm0_fus,Cm0_wing, Cl ,sweep_LE_fus, sweep_quarter_wing, force)[1],". It represents",Cm0(Cm0_fus,Cm0_wing, Cl ,sweep_LE_fus, sweep_quarter_wing, force)[1]/Cm0(Cm0_fus,Cm0_wing, Cl ,sweep_LE_fus, sweep_quarter_wing, force)[0]*100,"% of the pitching moment of the airplane.")
    print("----------------------------------------------------------------------")

    print("--------------------------TAIL--------------------------------------------")
    print("The position of the aerodynamic center of the tail is at",x_AC_tail," m and the end of the tail is at",x_AC_tail- x_AC_tail_local + c_root_tail,"m from the nose and the distance between x_AC_tail and x_CG_tot is",x_AC_tail - CG_position(config, fuel, Cl, sweep_LE_fus, sweep_quarter_wing, force)[0],"m.")
    print("----------------------------------------------------------------------")

    print("--------------------------FUEL STORAGE--------------------------------------------")
    print(CG_position(config, fuel, Cl, sweep_LE_fus, sweep_quarter_wing, force)[1]*100,"% of the total fuel volume is stored in the wings and we need",(1-CG_position(config, fuel, Cl, sweep_LE_fus, sweep_quarter_wing, force)[1])*26854.56/1000,"m³ in the fuselage.")
    print("----------------------------------------------------------------------")

    print("--------------------------CENTER OF GRAVITY--------------------------------------------")
    print("The center of gravity is positioned at",CG_position(config, fuel, Cl, sweep_LE_fus, sweep_quarter_wing, force)[0],"m (in feet :",CG_position(config, fuel, Cl, sweep_LE_fus, sweep_quarter_wing, force)[0]*3.28," (",CG_position(config, fuel, Cl, sweep_LE_fus, sweep_quarter_wing, force)[0]*100/l_fus,"%) from the nose of the airplane.")
    print("----------------------------------------------------------------------")

    print("--------------------------CRUISE--------------------------------------------")
    print("The new lift force generated by the body is",force,"[N]")
    print("The new lift force generated by the tail is",force_tail,"[N]")
    print("----------------------------------------------------------------------")

    print("--------------------------STATIC MARGIN AND NEUTRAL POINT--------------------------------------------")
    Kn, hn = long_stat_stab_cruise(config, fuel, Cl, sweep_LE_fus, sweep_quarter_wing, force)
    if Kn >= 0.05 and Kn < 0.3 : 
        print("The static margin has a correct value and is equal to : ", (Kn*100), "%  and the neutral point is positioned at",hn*MAC_tot,"from the nose, which represents",hn*MAC_tot*100/l_fus,"% of the total length.")
    else : 
        print("The static margin has to be changed and is equal to : ", (Kn*100), "%.  and the neutral point is positioned at",hn*MAC_tot,"m from the nose, which represents",hn*MAC_tot*100/l_fus,"% of the total length.")
    print("----------------------------------------------------------------------")

    print("--------------------------ACCEPTABLE POSITIONS FOR THE CENTER OF GRAVITY--------------------------------------------")
    print("The center of gravity is positioned at",get_CG(config, fuel,Cm0_fus,Cm0_wing,Kn, Cl, sweep_LE_fus, sweep_quarter_wing, force)[0],"m from the nose when the static margin is equal to",0.05*100,"%. It is the maximal value of the range.")
    print("The center of gravity is positioned at",get_CG(config, fuel,Cm0_fus,Cm0_wing,Kn, Cl, sweep_LE_fus, sweep_quarter_wing, force)[1],"m from the nose when the static margin is equal to",0.15*100,"%. It is the minimal value of the range.")
    print("----------------------------------------------------------------------")

    print("--------------------------DIRECTIONAL STABILITY--------------------------------------------")
    if (CN_beta_tot<0): #better to be > 0.1 and < 0.25
        print ("The aircraft is not directionally stable because CN_beta_tot is negative and equals",CN_beta_tot)
        print("The contribution of the fuselage is",CN_beta_fuselage)
        print("The contribution of the wings is",CN_beta_w)
        print("The contribution of the fin is",CN_beta_fin)
    else:
        print ("The aircraft is directionally stable because CN_beta_tot is positive and equals",CN_beta_tot)
        print("The contribution of the fuselage is",CN_beta_fuselage)
        print("The contribution of the wings is",CN_beta_w)
        print("The contribution of the fin is",CN_beta_fin)
    print("----------------------------------------------------------------------")

    print("--------------------------LATERAL STABILITY--------------------------------------------")
    print("The CL_beta coefficient is equal to",CL_beta_tot)
    print("The contribution of the wings is",CL_beta_wings)
    print("The contribution of the tail is",CL_beta_tail)
    print("----------------------------------------------------------------------")
    return

if __name__ == "__main__":
    main()
